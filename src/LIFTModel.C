/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of main solver class for LIFT process simulation.
    
    Handles:
    - Field initialization and evolution
    - Model coupling and synchronization
    - Energy conservation tracking
    - Phase change and interface dynamics
    - Solution of coupled equation systems
    
    The implementation follows a segregated solution approach:
    1. Phase equation for interface tracking
    2. Temperature equations for thermal evolution
    3. Momentum equation for fluid dynamics
    4. Pressure equation for incompressibility
    
    Special attention is paid to:
    - Energy conservation
    - Phase change accuracy
    - Interface sharpness
    - Numerical stability
\*---------------------------------------------------------------------------*/

#include "LIFTModel.H"
#include "twoTemperatureModel.H"
#include "femtosecondLaserModel.H"
#include "phaseChangeModel.H"
#include "dropletModel.H"
#include "fvOptions.H"
#include "MULES.H"
#include "dynamicFvMesh.H"


namespace Foam
{
    defineTypeNameAndDebug(LIFTModel, 0);
    
    // Define thermal conductivity dimensions
    const dimensionSet dimThermalConductivity(1, 1, -3, -1, 0, 0, 0);
}
/*---------------------------------------------------------------------------*\
                        LIFTModel Implementation
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LIFTModel::LIFTModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    //dynamicMesh_(new dynamicRefineFvMesh(runTime_)),  // Add this line
    dict_(dict),
    runTime_(const_cast<Time&>(mesh.time())),

    // Get field references
    U_(mesh.lookupObjectRef<volVectorField>("U")),
    phi_(mesh.lookupObjectRef<surfaceScalarField>("phi")),
    p_(mesh.lookupObjectRef<volScalarField>("p")),
    alpha1_(mesh.lookupObjectRef<volScalarField>("alpha.titanium")),
    rho_(mesh.lookupObjectRef<volScalarField>("rho")),
    Te_(mesh.lookupObjectRef<volScalarField>("Te")),
    Tl_(mesh.lookupObjectRef<volScalarField>("Tl")),
    
    // Create laser and phase change fields
    laserSource_
    (
        IOobject("laserSource", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimPower/dimVolume, 0.0)
    ),
    phaseChangeRate_
    (
        IOobject("phaseChangeRate", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    ),
    phaseIndicator_
    (
        IOobject("phaseIndicator", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    // Create interface tracking fields
    interfaceEnergy_
    (
        IOobject("interfaceEnergy", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
    ),
    interfaceNormal_
    (
        IOobject("interfaceNormal", runTime_.timeName(), mesh),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    curvature_
    (
        IOobject("curvature", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0.0)
    ),

    // Create pressure fields
    p_rgh_
    (
        IOobject("p_rgh", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    rhoPhi_
    (
        IOobject("rhoPhi", runTime_.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),
    gh_("gh", mesh.C() & g_),
    ghf_("ghf", mesh.Cf() & g_),

    // Initialize models
    mixturePtr_(new immiscibleIncompressibleTwoPhaseMixture(U_, phi_)),
    mixture_(*mixturePtr_),
    ttm_(new twoTemperatureModel(mesh_, dict_)),
    laser_(new femtosecondLaserModel(mesh_, dict_)),
    phaseChange_(new phaseChangeModel(mesh_, Te_, alpha1_, dict_, true)),
    droplet_(new dropletModel(mesh_, dict_, rho_, U_, true)),
    turbulence_(nullptr),
       MRF_(mesh),
    fvOptions_(fv::options::New(mesh)),
    UEqn_(U_, dimForce/dimVolume),
// In LIFTModel constructor
interfaceCapturing_
(
    new advancedInterfaceCapturing
    (
        mesh,
        mesh_.lookupObjectRef<volScalarField>("alpha.titanium"),  // Use lookupObjectRef instead
        phi_,
        mixture_,
        Tl_
    )
),
    radiation_(new radiationModel(mesh, dict)),
    thermo_(basicThermo::New(mesh)),

    // Physical parameters from dictionary
    Tm_("Tm", dimTemperature, dict.get<scalar>("meltingTemperature")),
    Tv_("Tv", dimTemperature, dict.get<scalar>("vaporizationTemperature")),
    Lm_("Lm", dimEnergy/dimMass, dict.get<scalar>("latentHeatMelting")),
    Lv_("Lv", dimEnergy/dimMass, dict.get<scalar>("latentHeatVaporization")),
    sigma_("sigma", dimForce/dimLength, dict.get<scalar>("surfaceTension")),
    interfaceWidth_("interfaceWidth", dimLength, dict.get<scalar>("interfaceWidth")),
    criticalPressure_("pc", dimPressure, dict.get<scalar>("criticalPressure")),
    g_("g", dimAcceleration, dict.getOrDefault<vector>("g", vector(0, -9.81, 0))),

    // Material properties
    rho1_("rho1", dimDensity, dict.subDict("donor").get<scalar>("density")),
    rho2_("rho2", dimDensity, dict.subDict("ambient").get<scalar>("density")),
    mu1_("mu1", dimDynamicViscosity, dict.subDict("donor").get<scalar>("viscosity")),
    mu2_("mu2", dimDynamicViscosity, dict.subDict("ambient").get<scalar>("viscosity")),
    Cp1_("Cp1", dimSpecificHeatCapacity, dict.subDict("donor").get<scalar>("Cp")),
    Cp2_("Cp2", dimSpecificHeatCapacity, dict.subDict("ambient").get<scalar>("Cp")),
    k1_("k1", dimThermalConductivity, dict.subDict("donor").get<scalar>("thermalConductivity")),
    k2_("k2", dimThermalConductivity, dict.subDict("ambient").get<scalar>("thermalConductivity")),

    // Solution controls
    pimple_(const_cast<pimpleControl&>(mesh.lookupObject<pimpleControl>("PIMPLE"))),
    maxCo_(dict.getOrDefault<scalar>("maxCo", 0.5)),
    alphaCoNum_(0.0),
    pRefCell_(0),
    pRefValue_(0.0),
    cumulativeContErr_(0.0)
{
    // Initialize turbulence model
    turbulence_ = incompressible::turbulenceModel::New(U_, phi_, mixture_);
// In constructor
/*if (!mesh_.foundObject<volScalarField>("alpha1"))
{
    mesh_.objectRegistry::store
    (
        new volScalarField("alpha1", alpha1_)
    );
}
    // Set phase field
    dynamicMesh_->setPhaseField(alpha1_);
*/
    if (debug)
    {
        Info<< "Created alpha1 alias for alpha.titanium" << endl;
    }
    Info<< "\nInitialized LIFT model with:" << nl
        << "  Melting temperature: " << Tm_.value() << " K" << nl
        << "  Vaporization temperature: " << Tv_.value() << " K" << nl
        << "  Latent heat of melting: " << Lm_.value() << " J/kg" << nl
        << "  Latent heat of vaporization: " << Lv_.value() << " J/kg" << nl
        << "  Surface tension: " << sigma_.value() << " N/m" << nl
        << "  Interface width: " << interfaceWidth_.value() << " m" << endl;
}

Foam::autoPtr<Foam::LIFTModel> 
Foam::LIFTModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    return autoPtr<LIFTModel>(new LIFTModel(mesh, dict));
}

Foam::dimensionedScalar 
Foam::LIFTModel::calculateTotalEnergy() const
{
    // Calculate thermal energy
    const dimensionedScalar thermalEnergy = fvc::domainIntegrate
    (
        (ttm_->Ce() * Te_ + ttm_->Cl() * Tl_) * rho_
    );

    // Calculate kinetic energy
    const dimensionedScalar kineticEnergy = fvc::domainIntegrate
    (
        0.5 * rho_ * magSqr(U_)
    );

    // Calculate interface energy
    const dimensionedScalar surfaceEnergy = fvc::domainIntegrate(interfaceEnergy_);

    // Calculate phase change energy
const dimensionedScalar phaseEnergy = fvc::domainIntegrate
(
    mesh_.lookupObject<volScalarField>("alpha.titanium") * rho_ * 
    (
        Lm_ * pos(Tl_ - Tm_) + 
        Lv_ * pos(Tl_ - Tv_)
    )
);

    // Get laser contribution if any
    const dimensionedScalar laserEnergy = fvc::domainIntegrate
    (
        laserSource_ * mesh_.time().deltaT()
    );

    // Return total energy
    return thermalEnergy + kineticEnergy + surfaceEnergy + 
           phaseEnergy + laserEnergy;
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void LIFTModel::solve()
{
    // Solve phase equation
    solvePhaseEquation();
    
    // Solve temperature equations
    solveTemperatureEquation();
    
    // Solve momentum equation
    solveMomentumEquation(); 
    
    // Solve pressure equation  
    solvePressureEquation();
    
    // Update models  
    ttm_->solve(laser_->source());
    phaseChange_->correct(Tl_);
    
    if (droplet_)
    {
        droplet_->update(Tl_, p_, U_);
    }
}

void LIFTModel::solvePhaseEquation()
{
    // Track initial energy
    dimensionedScalar energyBefore = calculateTotalEnergy();

    // Solve alpha equation with phase change
    surfaceScalarField alphaPhi
    (
        "alphaPhi",
        phi_*fvc::interpolate(alpha1_)
      + fvc::interpolate(phaseChangeRate_/mixture_.rho1())
    );

    // Apply MULES limiter with explicit bounds
    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1_,
        phi_,
        alphaPhi,
        zeroField(),
        zeroField(),
        scalarField(alpha1_.size(), 1.0),  // Replace alphaMax
        scalarField(alpha1_.size(), 0.0)   // Replace alphaMin
    );

    // Update mixture properties
    mixture_.correct();
    rho_ = alpha1_*mixture_.rho1() + (1.0 - alpha1_)*mixture_.rho2();

    // Check energy conservation
    checkEnergyConservation(energyBefore, calculateTotalEnergy());
}

void LIFTModel::solveTemperatureEquation()
{
    // Initial energy tracking
    dimensionedScalar energyBefore = calculateTotalEnergy();

    // Get temperature fields and properties from two-temp model
    tmp<volScalarField> tke = ttm_->ke();
    tmp<volScalarField> tkl = ttm_->kl();
    
    // Update laser source contribution
    laser_->correct();
    
    // Electron temperature equation with electron-phonon coupling
    fvScalarMatrix TeEqn
    (
        fvm::ddt(ttm_->Ce()*rho_, Te_)
      - fvm::laplacian(tke(), Te_)
      + fvm::Sp(ttm_->G()*rho_, Te_)
      - ttm_->G()*rho_*Tl_
     ==
        laser_->source()
      + phaseChange_->electronSource()
    );

    TeEqn.relax();
    TeEqn.solve();

    // Lattice temperature equation with phase change contribution
    fvScalarMatrix TlEqn
    (
        fvm::ddt(ttm_->Cl()*rho_, Tl_)
      + fvm::div(phi_, Tl_)
      - fvm::laplacian(tkl(), Tl_)
      - ttm_->G()*rho_*Te_
      + fvm::Sp(ttm_->G()*rho_, Tl_)
     ==
        phaseChange_->latticeSource()
    );

    TlEqn.relax();
    TlEqn.solve();

    // Energy conservation check
    checkEnergyConservation(energyBefore, calculateTotalEnergy());
}

void LIFTModel::solveMomentumEquation() 
{
    // Surface tension with temperature effects
    volScalarField localSigma(sigma_*(1.0 - 0.26e-3*(Te_ - Tm_)));

    // Calculate effective viscosity
    volScalarField muEff(mixture_.mu());
    if (turbulence_)
    {
        muEff += rho_*turbulence_->nut();
    }

    // Momentum equation
    UEqn_ = 
    (
        fvm::ddt(rho_, U_)
      + fvm::div(rhoPhi_, U_)    
      - fvm::laplacian(muEff, U_)
     ==
        fvc::reconstruct
        (
            mixture_.surfaceTensionForce()*mesh_.magSf() +
            (
              - ghf_*fvc::snGrad(rho_) 
              - fvc::snGrad(p_rgh_)
            )*mesh_.magSf()
        )
      + phaseChange_->momentumSource()
    );

    UEqn_.relax();

    if (pimple_.momentumPredictor())
    {
        // Use fvSolution::solve() instead of member solve()
        fvVectorMatrix UEqnf = UEqn_;
        UEqnf.source() -= fvc::grad(p_rgh_);
        UEqnf.solve();
    }
}

void LIFTModel::solvePressureEquation()
{
    volScalarField rAU(1.0/UEqn_.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U_);
    HbyA = rAU*UEqn_.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U_, phi_)
    );

    adjustPhi(phiHbyA, U_, p_rgh_);

    // Non-orthogonal pressure corrector loop
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh_) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell_, pRefValue_);
        p_rghEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA - p_rghEqn.flux();
            U_ = HbyA - rAU*fvc::reconstruct(p_rghEqn.flux()/rAUf);
            p_ = p_rgh_ + rho_*gh_;
        }
    }
}
bool 
Foam::LIFTModel::valid() const 
{
    return ttm_.valid() && laser_.valid() && mixturePtr_.valid();
}
bool LIFTModel::checkEnergyConservation
(
    const dimensionedScalar& energyBefore,
    const dimensionedScalar& energyAfter
) const
{
    scalar energyError = mag
    (
        (energyAfter.value() - energyBefore.value())/
        (mag(energyBefore.value()) + SMALL)
    );
    
    if (energyError > 1e-6)
    {
        Info<< "Energy conservation error: " << energyError*100 << "%" << endl;
        return false;
    }
    return true;
}

void 
Foam::LIFTModel::write()
{
    if (mesh_.time().writeTime())
    {
        U_.write();
        phi_.write();
        p_.write();
        Te_.write();
        Tl_.write();
                mesh_.lookupObject<volScalarField>("alpha.titanium").write();  // Changed from alpha1_
        rho_.write();
        laserSource_.write();
        phaseChangeRate_.write();
    }
}
