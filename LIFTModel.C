#include "LIFTModel.H"
#include "DimensionValidator.H"
#include "mathematicalConstants.H"
#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "fvc.H"
#include "fvm.H"
#include "findRefCell.H"
#include "adjustPhi.H"

namespace Foam 
{

LIFTModel::LIFTModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    mixture_(const_cast<immiscibleIncompressibleTwoPhaseMixture&>
        (mesh.lookupObject<immiscibleIncompressibleTwoPhaseMixture>("mixture"))),
    reporter_(std::make_unique<ValidationReporter>(mesh.time(), dict)),
    
    // Initialize fields
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.225)
    ),
    U_
    (
        IOobject("U", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0, 0, 0))
    ),
    p_
    (
        IOobject("p", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("p", dimPressure, 1e5)
    ),
    alpha1_
    (
        IOobject("alpha.titanium", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("alpha", dimless, 0)
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    Te_
    (
        IOobject("Te", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("Te", dimTemperature, 300)
    ),
    Tl_
    (
        IOobject("Tl", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("Tl", dimTemperature, 300)
    ),
    laserSource_
    (
        IOobject("laserSource", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("laserSource", dimPower/dimVolume, 0)
    ),
    phaseChangeRate_
    (
        IOobject("phaseChangeRate", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("phaseChangeRate", dimMass/dimVolume/dimTime, 0)
    ),
    phaseIndicator_
    (
        IOobject("phaseIndicator", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("phaseIndicator", dimless, 0)
    ),
    interfaceEnergy_
    (
        IOobject("interfaceEnergy", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("interfaceEnergy", dimEnergy/dimArea, 0)
    ),

    // Initialize models with make_unique
    ttm_(std::make_unique<twoTemperatureModel>(mesh, dict)),
    laser_(std::make_unique<femtosecondLaserModel>(mesh, dict)),
    shockWave_(std::make_unique<ultraFastShockWaveModel>(mesh, dict)),
    droplet_(std::make_unique<dropletModel>(mesh, dict, rho_, U_)),
    fsDroplet_(std::make_unique<femtosecondDropletModel>(mesh, dict, rho_, U_)),
    material_(std::make_unique<extremeConditionMaterialProperties>(mesh, dict)),
    phaseChange_(std::make_unique<nonEquilibriumPhaseChangeModel>(mesh, Tl_, alpha1_, dict)),
    analytics_(std::make_unique<timeResolvedAnalytics>(mesh, dict)),
    visualization_(std::make_unique<multiScaleVisualization>(mesh, dict)),
    timeStep_(std::make_unique<adaptiveTimeStep>(mesh, dict)),
    stressModel_(std::make_unique<residualStressModel>(Tl_, alpha1_, dict)),
    mixtureProps_(std::make_unique<twoPhaseMixtureProperties>(mesh, dict)),
    interfaceCapturing_(std::make_unique<advancedInterfaceCapturing>(mesh, alpha1_, phi_, mixture_)),
    curvature_(std::make_unique<curvatureModel>(mesh, alpha1_)),
    dropletMetrics_(std::make_unique<dropletMetrics>(mesh, alpha1_)),
    radiation_(std::make_unique<radiationModel>(mesh, dict)),
    meshRefinement_(nullptr),
    vacuumBC_(nullptr),

    // Initialize physical parameters
    Tm_(dimensionedScalar("Tm", dimTemperature, dict.get<scalar>("meltingTemperature"))),
    Tv_(dimensionedScalar("Tv", dimTemperature, dict.get<scalar>("vaporizationTemperature"))),
    Lm_(dimensionedScalar("Lm", dimEnergy/dimMass, dict.get<scalar>("latentHeatMelting"))),
    Lv_(dimensionedScalar("Lv", dimEnergy/dimMass, dict.get<scalar>("latentHeatVaporization"))),
    sigma_(dimensionedScalar("sigma", dimForce/dimLength, dict.get<scalar>("surfaceTension"))),
    interfaceWidth_(dimensionedScalar("interfaceWidth", dimLength, dict.get<scalar>("interfaceWidth"))),
    
    // Initialize energy tracking
    initialEnergy_(calculateTotalEnergy()),
    lastEnergy_(initialEnergy_),
    maxEnergyChange_(dict.lookupOrDefault<scalar>("maxEnergyChange", 0.01)),

    // Initialize solution controls
    solveEnergy_(dict.lookupOrDefault("solveEnergy", true)),
    solveMomentum_(dict.lookupOrDefault("solveMomentum", true)),
    adaptiveMesh_(dict.lookupOrDefault("adaptiveMesh", false)),
    writeFields_(dict.lookupOrDefault("writeFields", true))
{
    // Validate mixture
    if (!validateMixture())
    {
        reportError("Invalid mixture initialization");
        FatalError << "Model initialization failed" << abort(FatalError);
    }

    // Initialize fields and validate dimensions
    initializeFields();
    validateDimensions();

    // Initialize laser source
    initializeLaserSource();

    // Initialize material regions
    initializeMaterialRegions();

    // Set up mesh refinement if requested
    if (adaptiveMesh_)
    {
        if (dynamicRefineFvMesh* refineMesh = dynamic_cast<dynamicRefineFvMesh*>(const_cast<fvMesh*>(&mesh_)))
        {
            meshRefinement_ = std::make_unique<adaptiveMeshRefinement>(*refineMesh, alpha1_);
        }
        else
        {
            reportWarning("Adaptive mesh requested but mesh is not dynamic. Disabling adaptive mesh.");
            adaptiveMesh_ = false;
        }
    }

    // Set up vacuum boundary condition if patch exists
    label patchIndex = mesh.boundaryMesh().findPatchID("vacuum");
    if (patchIndex >= 0)
    {
        const fvPatch& patch = mesh.boundary()[patchIndex];
        vacuumBC_ = std::make_unique<vacuumBoundaryCondition>
        (
            patch,
            mesh.lookupObject<volScalarField>("p"),
            dict
        );
    }
    else
    {
        Info<< "Warning: 'vacuum' patch not found. Using default boundary conditions." << endl;
    }

    // Final validation
    if (!valid())
    {
        reportError("Model validation failed during initialization");
        FatalError << "Model initialization failed" << abort(FatalError);
    }

    // Report successful initialization
    reporter_->addResult
    (
        "initialization",
        true,
        0.0,
        0.0,
        "Model initialized successfully",
        "all"
    );

    // Write initial state if requested
    if (writeFields_)
    {
        write();
    }
}
// ... (continued from previous part)

void LIFTModel::solve()
{
    // Pre-solve validation
    if (!validateFields())
    {
        reportError("Pre-solve field validation failed");
        return;
    }

    // Store initial energy for conservation check
    dimensionedScalar energyBefore = calculateTotalEnergy();

    try
    {
        // Energy equations
        if (solveEnergy_)
        {
            // Laser-material interaction
            laser_->update();
            ttm_->solve(laser_->source());
            
            // Phase change and material properties
            phaseChange_->correct(Tl_);
            material_->update(Tl_);
            
            // Update radiation model
            radiation_->correct(Tl_);
        }

        // Momentum equations
        if (solveMomentum_)
        {
            // Calculate interface energy
            calculateInterfaceEnergy();
            
            // Interface tracking
            interfaceCapturing_->correct();
            
            // Update physics models
            shockWave_->update(p_, U_);
            droplet_->update(Tl_, p_, U_);
            fsDroplet_->update(Tl_, p_, U_);
            mixtureProps_->correct(alpha1_);
        }
        
        // Mesh refinement if enabled
        if (adaptiveMesh_ && meshRefinement_)
        {
            meshRefinement_->refine();
            if (!validateMeshQuality())
            {
                reportWarning("Mesh quality metrics exceeded after refinement");
            }
        }

        // Update phase indicator
        updatePhaseIndicator();

        // Analysis and visualization
        stressModel_->calculate(Tl_, alpha1_);
        analytics_->calculate(Tl_, p_);
        visualization_->update(Tl_, alpha1_);

        // Post-solve energy conservation check
        dimensionedScalar energyAfter = calculateTotalEnergy();
        if (!checkEnergyConservation(energyBefore, energyAfter))
        {
            reportWarning("Energy conservation violated");
        }

        // Update energy tracking
        updateEnergyTracking();

        // Write fields if requested
        if (writeFields_ && mesh_.time().writeTime())
        {
            write();
            writeSummary();
        }
    }
    catch (const Foam::error& e)
    {
        reportError("Error in solve step: " + std::string(e.message()));
        throw;
    }
}

void LIFTModel::writeSummary() const
{
    Info<< "Time = " << mesh_.time().timeName() << nl
        << "  Phase fractions:" << nl
        << "    Min alpha = " << min(alpha1_).value() << nl
        << "    Max alpha = " << max(alpha1_).value() << nl
        << "  Temperatures:" << nl
        << "    Min Te = " << min(Te_).value() << nl
        << "    Max Te = " << max(Te_).value() << nl
        << "    Min Tl = " << min(Tl_).value() << nl
        << "    Max Tl = " << max(Tl_).value() << nl
        << "  Energy:" << nl
        << "    Total energy = " << lastEnergy_.value() << nl
        << "    Energy change = " << 
           (lastEnergy_.value() - initialEnergy_.value())/initialEnergy_.value() * 100
        << "%" << endl;

    if (droplet_->isDropletFormed())
    {
        Info<< "  Droplet metrics:" << nl
            << "    Aspect ratio = " << dropletMetrics_->aspectRatio() << nl
            << "    Circularity = " << dropletMetrics_->circularity() << endl;
    }
}

void LIFTModel::updateFields
(
    volScalarField& alpha1,
    volVectorField& U,
    volScalarField& p,
    volScalarField& p_rgh,
    volScalarField& rho
)
{
    alpha1 = alpha1_;
    U = U_;
    p = p_;
    rho = rho_;
    
    // Update p_rgh based on hydrostatic pressure
    p_rgh = p_ - rho_*mesh_.C().component(vector::Z)*dimensionedScalar("g", dimAcceleration, 9.81);
}

void LIFTModel::calculateInterfaceEnergy()
{
    volVectorField gradAlpha = fvc::grad(alpha1_);
    interfaceEnergy_ = sigma_ * mag(gradAlpha);
}
// Add these implementations to the LIFTModel.C file

dimensionedScalar LIFTModel::calculateTotalEnergy() const
{
    // Thermal energy (electron and lattice subsystems)
    dimensionedScalar thermalEnergy = fvc::domainIntegrate
    (
        (ttm_->Ce() * Te_ + ttm_->Cl() * Tl_) * rho_
    );

    // Kinetic energy
    dimensionedScalar kineticEnergy = fvc::domainIntegrate
    (
        0.5 * rho_ * magSqr(U_)
    );

    // Surface energy from interface
    dimensionedScalar surfaceEnergy = fvc::domainIntegrate(interfaceEnergy_);

    // Phase change energy
    dimensionedScalar phaseEnergy = fvc::domainIntegrate
    (
        alpha1_ * rho_ * (
            Lm_ * pos(Tl_ - Tm_) +  // Melting contribution
            Lv_ * pos(Tl_ - Tv_)    // Vaporization contribution
        )
    );

    return thermalEnergy + kineticEnergy + surfaceEnergy + phaseEnergy;
}

void LIFTModel::updatePhaseIndicator()
{
    forAll(mesh_.C(), cellI)
    {
        if (Tl_[cellI] < Tm_.value())
        {
            phaseIndicator_[cellI] = 0;  // Solid
        }
        else if (Tl_[cellI] < Tv_.value())
        {
            phaseIndicator_[cellI] = 1;  // Liquid
            
            // Add melting transition zone
            scalar meltFraction = (Tl_[cellI] - Tm_.value())/(Tv_.value() - Tm_.value());
            phaseIndicator_[cellI] = meltFraction;
        }
        else
        {
            phaseIndicator_[cellI] = 2;  // Vapor
        }
    }
}

void LIFTModel::initializeMaterialRegions()
{
    // Get geometry parameters from dictionary
    const scalar donorThickness = dict_.get<scalar>("donorThickness");
    const scalar substrateThickness = dict_.get<scalar>("substrateThickness");
    const scalar gapThickness = dict_.lookupOrDefault<scalar>("gapThickness", 10e-6);
    //const vector donorCenter = dict_.get<vector>("donorCenter");

    forAll(mesh_.C(), cellI)
    {
        const vector& pos = mesh_.C()[cellI];
        const scalar z = pos.z();
        //const scalar r = mag(vector(pos.x() - donorCenter.x(), pos.y() - donorCenter.y(), 0));
        
        // Donor layer (titanium)
        if (z < donorThickness)
        {
            alpha1_[cellI] = 1.0;
            Te_[cellI] = dict_.lookupOrDefault<scalar>("donorInitialTemp", 300.0);
            Tl_[cellI] = Te_[cellI];
        }
        // Gap region
        else if (z < donorThickness + gapThickness)
        {
            alpha1_[cellI] = 0.0;
            Te_[cellI] = dict_.lookupOrDefault<scalar>("ambientTemp", 300.0);
            Tl_[cellI] = Te_[cellI];
        }
        // Substrate (glass)
        else if (z < donorThickness + gapThickness + substrateThickness)
        {
            alpha1_[cellI] = 0.0;
            Te_[cellI] = dict_.lookupOrDefault<scalar>("substrateInitialTemp", 300.0);
            Tl_[cellI] = Te_[cellI];
        }
        // Surrounding medium
        else
        {
            alpha1_[cellI] = 0.0;
            Te_[cellI] = dict_.lookupOrDefault<scalar>("ambientTemp", 300.0);
            Tl_[cellI] = Te_[cellI];
        }
        
        // Initialize phase indicator based on initial temperature
        if (Tl_[cellI] < Tm_.value())
        {
            phaseIndicator_[cellI] = 0;
        }
        else if (Tl_[cellI] < Tv_.value())
        {
            phaseIndicator_[cellI] = 1;
        }
        else
        {
            phaseIndicator_[cellI] = 2;
        }
    }
}

bool LIFTModel::checkEnergyConservation
(
    const dimensionedScalar& energyBefore,
    const dimensionedScalar& energyAfter
) const
{
    // Calculate relative energy change
    scalar relativeChange = mag
    (
        (energyAfter.value() - energyBefore.value())/
        (energyBefore.value() + SMALL)
    );

    // Get energy conservation parameters
    const scalar maxRelativeChange = 
        dict_.lookupOrDefault<scalar>("maxEnergyChange", 0.01);
    const scalar absoluteTolerance = 
        dict_.lookupOrDefault<scalar>("energyAbsTolerance", 1e-6);

    // Check both relative and absolute changes
    bool conserved = 
        (relativeChange < maxRelativeChange) || 
        (mag(energyAfter.value() - energyBefore.value()) < absoluteTolerance);

    if (!conserved)
    {
        Info<< "Energy conservation warning:" << nl
            << "  Before: " << energyBefore.value() << nl
            << "  After:  " << energyAfter.value() << nl
            << "  Change: " << relativeChange * 100 << "%" << nl
            << "  Limit:  " << maxRelativeChange * 100 << "%" << endl;
    }

    return conserved;
}

bool LIFTModel::valid() const
{
    bool isValid = true;

    // Check basic requirements
    isValid &= validateMixture();
    isValid &= validateFields();
    isValid &= validateMeshQuality();
    
    // Check models
    if (!ttm_ || !ttm_->valid())
    {
        reportError("Invalid two-temperature model");
        isValid = false;
    }

    if (!laser_ || !laser_->valid())
    {
        reportError("Invalid laser model");
        isValid = false;
    }

    if (!phaseChange_ || !phaseChange_->valid())
    {
        reportError("Invalid phase change model");
        isValid = false;
    }

    if (!material_ || !material_->valid())
    {
        reportError("Invalid material properties model");
        isValid = false;
    }

    if (!mixtureProps_ || !mixtureProps_->valid())
    {
        reportError("Invalid mixture properties model");
        isValid = false;
    }

    // Check physical parameters
    if (Tm_.value() <= 0 || Tv_.value() <= Tm_.value())
    {
        reportError("Invalid temperature limits");
        isValid = false;
    }

    if (Lm_.value() <= 0 || Lv_.value() <= 0)
    {
        reportError("Invalid latent heat values");
        isValid = false;
    }

    if (sigma_.value() <= 0)
    {
        reportError("Invalid surface tension");
        isValid = false;
    }

    if (interfaceWidth_.value() <= 0)
    {
        reportError("Invalid interface width");
        isValid = false;
    }

    return isValid;
}
bool LIFTModel::validateMixture() const
{
    bool valid = true;

    // Check mixture properties - update to use more specific checks
    if (!(mixture_.alpha1().size() == mesh_.nCells()))
    {
        reportError("Invalid mixture alpha1 field size");
        valid = false;
    }

    // Check viscosity limits
    dimensionedScalar maxViscosity("maxViscosity", dimViscosity, 1e3);
    dimensionedScalar minViscosity("minViscosity", dimViscosity, 1e-6);

    if (max(mixture_.nu()).value() > maxViscosity.value() || 
        min(mixture_.nu()).value() < minViscosity.value())
    {
        reportWarning("Mixture viscosity out of normal range");
    }

    return valid;
}

bool LIFTModel::validateFields() const
{
    bool valid = true;

    // Check field sizes
    if (rho_.size() != mesh_.nCells() || U_.size() != mesh_.nCells() ||
        p_.size() != mesh_.nCells() || alpha1_.size() != mesh_.nCells() ||
        Te_.size() != mesh_.nCells() || Tl_.size() != mesh_.nCells())
    {
        reportError("Field size mismatch");
        valid = false;
    }

    // Check phase fraction bounds
    if (min(alpha1_).value() < -SMALL || max(alpha1_).value() > 1.0 + SMALL)
    {
        reportError("Phase fraction out of bounds");
        valid = false;
    }

    // Check temperature bounds
    dimensionedScalar minTemp("minTemp", dimTemperature, 200);  // Realistic minimum
    dimensionedScalar maxTemp("maxTemp", dimTemperature, 1e5);  // Upper limit

    if (min(Te_).value() < minTemp.value() || max(Te_).value() > maxTemp.value() ||
        min(Tl_).value() < minTemp.value() || max(Tl_).value() > maxTemp.value())
    {
        reportError("Temperature out of physical bounds");
        valid = false;
    }

    // Check pressure bounds
    dimensionedScalar minPressure("minPressure", dimPressure, 10);  // Near vacuum
    if (min(p_).value() < minPressure.value())
    {
        reportError("Pressure below minimum threshold");
        valid = false;
    }

    // Check velocity bounds
    dimensionedScalar maxVelocity("maxVelocity", dimVelocity, 1e4);  // Reasonable maximum
    if (max(mag(U_)).value() > maxVelocity.value())
    {
        reportError("Velocity magnitude exceeds maximum threshold");
        valid = false;
    }

    return valid;
}

void LIFTModel::initializeFields()
{
    // Initialize phi based on velocity
    phi_ = fvc::interpolate(U_) & mesh_.Sf();

    // Initialize phase fractions
    const scalar donorThickness = dict_.get<scalar>("donorThickness");
    const scalar substrateThickness = dict_.get<scalar>("substrateThickness");

    forAll(mesh_.C(), cellI)
    {
        const scalar z = mesh_.C()[cellI].z();
        
        if (z < donorThickness)
        {
            alpha1_[cellI] = 1.0;  // Donor material (titanium)
            Te_[cellI] = 300.0;    // Room temperature
            Tl_[cellI] = 300.0;
        }
        else if (z < (donorThickness + substrateThickness))
        {
            alpha1_[cellI] = 0.0;  // Substrate (glass)
            Te_[cellI] = 300.0;
            Tl_[cellI] = 300.0;
        }
        else
        {
            alpha1_[cellI] = 0.0;  // Air/vacuum region
            Te_[cellI] = 300.0;
            Tl_[cellI] = 300.0;
        }
    }

    // Update density based on phase fractions
    rho_ = alpha1_*mixture_.rho1() + (1.0 - alpha1_)*mixture_.rho2();

    // Initialize interface energy
    calculateInterfaceEnergy();
}

void LIFTModel::initializeLaserSource()
{
    if (!laser_)
    {
        reportError("Laser model not initialized");
        return;
    }

    // Get laser parameters from dictionary
    const scalar pulseEnergy = dict_.get<scalar>("pulseEnergy");
    const scalar pulseWidth = dict_.get<scalar>("pulseWidth");
    const scalar spotSize = dict_.get<scalar>("spotSize");
    const vector laserPosition = dict_.get<vector>("laserPosition");

    // Calculate initial laser intensity distribution
    forAll(mesh_.C(), cellI)
    {
        const point& cellCenter = mesh_.C()[cellI];
        const scalar r = mag(cellCenter - laserPosition);
        
        // Gaussian spatial profile
        const scalar spatialFactor = exp(-2.0*sqr(r/spotSize));
        
        // Initial temporal profile (t = 0)
        const scalar temporalFactor = exp(-4.0*log(2.0)*sqr(0.0/pulseWidth));
        
        laserSource_[cellI] = pulseEnergy * spatialFactor * temporalFactor;
    }
}

void LIFTModel::validateDimensions() const
{
    // Check field dimensions (using the template version for GeometricFields)
    DimensionValidator::checkFieldDimensions(Te_, dimTemperature, "electron temperature");
    DimensionValidator::checkFieldDimensions(Tl_, dimTemperature, "lattice temperature");
    DimensionValidator::checkFieldDimensions(p_, dimPressure, "pressure");
    DimensionValidator::checkFieldDimensions(U_, dimVelocity, "velocity");
    DimensionValidator::checkFieldDimensions(alpha1_, dimless, "phase fraction");
    DimensionValidator::checkFieldDimensions(laserSource_, dimPower/dimVolume, "laser source");
    DimensionValidator::checkFieldDimensions(phaseChangeRate_, dimMass/dimVolume/dimTime, "phase change rate");

    // Check TTM model dimensions
    if (ttm_)
    {
        // For dimensioned quantities
        const word ceFieldName = "electron heat capacity";
        const word clFieldName = "lattice heat capacity";
        const word gFieldName = "electron-phonon coupling";

        DimensionValidator::checkFieldDimensions
        (
            ttm_->Ce().dimensions(),
            dimEnergy/dimVolume/dimTemperature,
            ceFieldName
        );

        DimensionValidator::checkFieldDimensions
        (
            ttm_->Cl().dimensions(),
            dimEnergy/dimVolume/dimTemperature,
            clFieldName
        );

        // For the coupling term, get dimensions from the tmp field safely
        tmp<volScalarField> GPtr = ttm_->electronPhononCoupling();
        DimensionValidator::checkFieldDimensions
        (
            GPtr().dimensions(),
            dimEnergy/dimVolume/dimTime/dimTemperature,
            gFieldName
        );
    }

    // Check physical parameters (using the dimensionSet version)
    DimensionValidator::checkFieldDimensions(Tm_.dimensions(), dimTemperature, "melting temperature");
    DimensionValidator::checkFieldDimensions(Tv_.dimensions(), dimTemperature, "vaporization temperature");
    DimensionValidator::checkFieldDimensions(Lm_.dimensions(), dimEnergy/dimMass, "latent heat of melting");
    DimensionValidator::checkFieldDimensions(Lv_.dimensions(), dimEnergy/dimMass, "latent heat of vaporization");
    DimensionValidator::checkFieldDimensions(sigma_.dimensions(), dimForce/dimLength, "surface tension");
}
scalar LIFTModel::calculateDeltaT() const
{
    if (!timeStep_)
    {
        return mesh_.time().deltaTValue();
    }

    scalar deltaT = timeStep_->calculateDeltaT(Te_, U_);

    // Additional constraints
    if (laser_)
    {
        // Resolve laser pulse
        const scalar pulseWidth = dict_.get<scalar>("pulseWidth");
        deltaT = min(deltaT, pulseWidth/10.0);
    }

    if (phaseChange_)
    {
        // Resolve phase change
        const scalar phaseChangeTime = 
            dict_.lookupOrDefault<scalar>("phaseChangeTimeScale", 1e-9);
        deltaT = min(deltaT, phaseChangeTime/5.0);
    }

    return deltaT;
}

void LIFTModel::write()
{
    if (!writeFields_)
    {
        return;
    }

    // Write primary fields
    rho_.write();
    U_.write();
    p_.write();
    alpha1_.write();
    phi_.write();

    // Write temperature fields
    Te_.write();
    Tl_.write();

    // Write process fields
    laserSource_.write();
    phaseChangeRate_.write();
    phaseIndicator_.write();
    interfaceEnergy_.write();

    // Write model data
    if (analytics_)
    {
        analytics_->write();
    }
    if (visualization_)
    {
        visualization_->write();
    }

    // Write summary
    writeSummary();
}

LIFTModel::~LIFTModel()
{
    // Write final state if configured
    if (writeFields_ && mesh_.time().writeTime())
    {
        write();
    }

    // Report final statistics
    Info<< "\nFinal state:" << nl
        << "  Total simulation time: " << mesh_.time().elapsedCpuTime() << " s" << nl
        << "  Total energy change: " << 
           (lastEnergy_.value() - initialEnergy_.value())/initialEnergy_.value() * 100
        << "%" << endl;
}

// Add these methods just before the final namespace closing brace:

void LIFTModel::reportError(const std::string& msg) const
{
    reporter_->addResult("error", false, 0.0, 0.0, msg, "none");
    FatalError << msg << endl;
}

void LIFTModel::reportWarning(const std::string& msg) const
{
    reporter_->addResult("warning", true, 0.0, 0.0, msg, "none");
    Warning << msg << endl;
}

bool LIFTModel::validateMeshQuality() const
{
    scalar maxSkewness = 0.0;
    scalar maxNonOrthogonality = 0.0;
    scalar minVolume = GREAT;
    
    // Get quality metrics from dictionary
    const scalar maxAllowedSkewness = 
        dict_.subDict("meshQuality").get<scalar>("maxSkewness");
    const scalar maxAllowedNonOrthogonality = 
        dict_.subDict("meshQuality").get<scalar>("maxNonOrthogonality");
    const scalar minAllowedVolume = 
        dict_.subDict("meshQuality").get<scalar>("minVolume");

    // Check cell metrics
    forAll(mesh_.cells(), cellI)
    {
        minVolume = min(minVolume, mesh_.V()[cellI]);
        
        // Calculate skewness and non-orthogonality
        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            label faceID = c[faceI];
            if (mesh_.isInternalFace(faceID))
            {
                const point& fc = mesh_.faceCentres()[faceID];
                const point& cc = mesh_.cellCentres()[cellI];
                const vector d = fc - cc;
                const vector s = mesh_.faceAreas()[faceID];
                
                // Calculate skewness
                scalar sk = mag(d & s)/(mag(d)*mag(s) + SMALL);
                maxSkewness = max(maxSkewness, sk);
                
                // Calculate non-orthogonality
                scalar no = acos(mag(d & s)/(mag(d)*mag(s) + SMALL))*180.0/constant::mathematical::pi;
                maxNonOrthogonality = max(maxNonOrthogonality, no);
            }
        }
    }

    bool valid = true;
    valid &= (maxSkewness < maxAllowedSkewness);
    valid &= (maxNonOrthogonality < maxAllowedNonOrthogonality);
    valid &= (minVolume > minAllowedVolume);

    if (!valid)
    {
        Info<< "Mesh quality metrics:" << nl
            << "  Non-orthogonality: " << maxNonOrthogonality << "/" << maxAllowedNonOrthogonality << nl
            << "  Skewness: " << maxSkewness << "/" << maxAllowedSkewness << nl
            << "  Min volume: " << minVolume << "/" << minAllowedVolume << endl;
    }

    return valid;
}


void LIFTModel::updateEnergyTracking()
{
    lastEnergy_ = calculateTotalEnergy();
}

} // End namespace Foam
