/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of phase change model for LIFT process simulation.
    
    Handles:
    - Phase transition calculations
    - Temperature-dependent properties
    - Residual stress computation
    - Energy conservation tracking
    - Source term generation
    
    The model supports both equilibrium and non-equilibrium phase changes,
    with temperature and pressure-dependent melting conditions.
    
\*---------------------------------------------------------------------------*/

#include "phaseChangeModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(phaseChangeModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
    Constructor
    Initializes phase change model with mesh, temperature field, phase fraction,
    dictionary settings and non-equilibrium flag.
\*---------------------------------------------------------------------------*/

phaseChangeModel::phaseChangeModel
(
    const fvMesh& mesh,
    const volScalarField& T,
    volScalarField& alpha1,
    const dictionary& dict,
    bool isNonEquilibrium
)
:
    mesh_(mesh),
    T_(T),
    alpha1_(alpha1), 
    dict_(dict),
    mixture_(nullptr),
    
    // Material properties
    Tm_("Tm", dimTemperature, dict.get<scalar>("meltingTemperature")),
    Tv_("Tv", dimTemperature, dict.get<scalar>("vaporizationTemperature")),
    Lm_("Lm", dimEnergy/dimMass, dict.get<scalar>("latentHeatMelting")),
    Lv_("Lv", dimEnergy/dimMass, dict.get<scalar>("latentHeatVaporization")), 
    interfaceWidth_("interfaceWidth", dimLength, dict.lookupOrDefault<scalar>("interfaceWidth", 5e-6)),
    dTmdp_("dTmdp", dimTemperature/dimPressure, dict.lookupOrDefault<scalar>("dTmdp", 0)),
    pRef_("pRef", dimPressure, dict.lookupOrDefault<scalar>("pRef", 1e5)),

    // Residual stress properties
    E_(dict.lookupOrDefault<scalar>("youngModulus", 2.0e11)),
    nu_(dict.lookupOrDefault<scalar>("poissonRatio", 0.3)),
    alpha_(dict.lookupOrDefault<scalar>("thermalExpansionCoeff", 1.2e-5)),
    Tref_(dict.lookupOrDefault<scalar>("referenceTemperature", 300.0)),

    // Kinetic parameters
    gradualMeltingRate_("meltRate", dimless/dimTime, dict.lookupOrDefault<scalar>("gradualMeltingRate", 1.0)),
    gradualSolidificationRate_("solidRate", dimless/dimTime, dict.lookupOrDefault<scalar>("gradualSolidificationRate", 1.0)),
    undercoolingCoeff_("undercooling", dimTemperature/dimLength, dict.lookupOrDefault<scalar>("undercoolingCoeff", 0)),

    // Initialize fields
    phaseIndicator_
    (
        IOobject("phaseIndicator", mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh,
                IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    ),
    phaseChangeRate_(),
    lastTotalEnergy_("energyTotal", dimEnergy, 0),
    energyInitialized_(false),
    isNonEquilibrium_(isNonEquilibrium)
{
    if (mesh_.foundObject<immiscibleIncompressibleTwoPhaseMixture>("mixture"))
    {
        mixture_ = &mesh_.lookupObject<immiscibleIncompressibleTwoPhaseMixture>("mixture");
        updateDensity();
    }

    if (!validate())
    {
        FatalErrorIn("phaseChangeModel::phaseChangeModel")
            << "Invalid phase change model configuration"
            << abort(FatalError);
    }
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
autoPtr<phaseChangeModel> phaseChangeModel::New
(
    const word& modelType,
    const fvMesh& mesh,
    const volScalarField& T,
    volScalarField& alpha1,
    const dictionary& dict
)
{
    // Default to non-equilibrium mode for femtosecond model
    bool isNonEquilibrium = (modelType == "femtosecondPhaseChange");
    
    // Create and return the model
    return autoPtr<phaseChangeModel>
    (
        new phaseChangeModel(mesh, T, alpha1, dict, isNonEquilibrium)
    );
}
/*---------------------------------------------------------------------------*\
    Calculate phase change rate based on temperature and phase state
\*---------------------------------------------------------------------------*/
// Add these helper functions at the start of the class
void Foam::phaseChangeModel::checkAndLimitTemperature(scalar& T) const
{
    if (T < 0)
    {
        T = 0;
        WarningInFunction
            << "Negative temperature " << T  
            << " detected and corrected to 0"
            << endl;
    }
    
    scalar maxTemp = dict_.lookupOrDefault<scalar>("maxTemperature", 5000.0);
    if (T > maxTemp)
    {
        T = maxTemp;
        WarningInFunction
            << "Temperature " << T 
            << " exceeded maximum limit and was capped at " << maxTemp
            << endl;
    }
}


bool Foam::phaseChangeModel::checkPropertyBounds() const
{
    if (Tm_.value() <= 0 || Tv_.value() <= Tm_.value())
    {
        FatalErrorInFunction
            << "Invalid temperature bounds: Tm = " << Tm_.value()
            << ", Tv = " << Tv_.value()
            << abort(FatalError);
        return false;
    }

    if (Lm_.value() <= 0 || Lv_.value() <= 0)
    {
        FatalErrorInFunction
            << "Invalid latent heat values: Lm = " << Lm_.value()
            << ", Lv = " << Lv_.value()
            << abort(FatalError);
        return false;
    }

    return true;
}

// Modify calculatePhaseChangeRate() to improve stability
void phaseChangeModel::calculatePhaseChangeRate() const
{
    if (!phaseChangeRate_.valid())
    {
        phaseChangeRate_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject("phaseChangeRate", mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
            )
        );
    }

    volScalarField& pcr = phaseChangeRate_.ref();
    
    // Calculate temperature gradient safely
    tmp<volVectorField> gradT = fvc::grad(T_);
    const volVectorField& gradTRef = gradT();
    
    // Maximum allowed phase change rate
    const scalar maxRate = dict_.lookupOrDefault<scalar>("maxPhaseChangeRate", 1e6);
    
    forAll(mesh_.C(), cellI)
    {
        // Check and limit temperature
        scalar T = T_[cellI];
        checkAndLimitTemperature(T);
        
        // Calculate effective melting temperature with safeguards
        scalar gradMag = mag(gradTRef[cellI]);
        if (gradMag > 1e6)
        {
            gradMag = 1e6;  // Limit extreme gradients
        }
        
        scalar Tm_effective = Tm_.value() - undercoolingCoeff_.value() * gradMag;
        
        // Ensure Tm_effective stays in physical range
        Tm_effective = max(Tm_effective, 0.5*Tm_.value());
        
        // Calculate local density safely
        scalar localDensity = 0.0;
        if (mixture_)
        {
            scalar alpha = max(min(alpha1_[cellI], 1.0), 0.0);
            localDensity = alpha*mixture_->rho1().value() + 
                          (1.0 - alpha)*mixture_->rho2().value();
        }
        else
        {
            localDensity = max(rho_[cellI], SMALL);
        }

        // Calculate phase change rate with limiting
        scalar rate = 0.0;
        
        if (isNonEquilibrium_)
        {
            if (T > Tv_.value())
            {
                rate = min(Lv_.value() * localDensity * gradualMeltingRate_.value(), maxRate);
            }
            else if (T > Tm_effective)
            {
                rate = min(Lm_.value() * localDensity * gradualMeltingRate_.value(), maxRate);
            }
        }
        else
        {
            if (T > Tm_effective)
            {
                rate = min(Lm_.value() * localDensity * gradualMeltingRate_.value(), maxRate);
            }
            else
            {
                rate = max(-Lm_.value() * localDensity * gradualSolidificationRate_.value(), -maxRate);
            }
        }

        // Apply rate with smooth transition
        const scalar transitionWidth = 1.0;  // Temperature range for smooth transition
        scalar blend = 0.5*(1.0 + tanh((T - Tm_effective)/transitionWidth));
        pcr[cellI] = rate * blend;
    }

    // Ensure field is bounded
    pcr.max(-maxRate);
    pcr.min(maxRate);
}

// Modify residualStress() calculation for better stability
tmp<volSymmTensorField> phaseChangeModel::residualStress() const
{
    tmp<volSymmTensorField> tStress
    (
        new volSymmTensorField
        (
            IOobject
            (
                "residualStress",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            symmTensor::zero
        )
    );

    volSymmTensorField& stress = tStress.ref();
    
    // Create dimensioned temperature difference with bounds
    dimensionedScalar minTemp("minTemp", dimTemperature, -1000.0);
    dimensionedScalar maxTemp("maxTemp", dimTemperature, 1000.0);
    dimensionedScalar Tref("Tref", dimTemperature, Tref_);

    volScalarField limitedTempDiff = 
        min
        (
            max
            (
                T_ - Tref,
                minTemp
            ),
            maxTemp
        );
        
    volScalarField thermalStrain = alpha_ * limitedTempDiff;

    // Calculate elastic constants with checks
    scalar lambda = 0.0;
    scalar mu = 0.0;
    
    if (nu_ > -1.0 && nu_ < 0.5)  // Valid Poisson's ratio range
    {
        lambda = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
        mu = E_ / (2.0 * (1.0 + nu_));
    }
    else
    {
        WarningInFunction
            << "Invalid Poisson's ratio. Using default values."
            << endl;
        lambda = E_ * 0.3;  // Default assuming nu = 0.3
        mu = E_ / 2.6;
    }

    // Maximum allowed stress with proper dimensions
    dimensionedScalar maxStress
    (
        "maxStress",
        dimPressure,
        dict_.lookupOrDefault<scalar>("maxStress", 1e9)
    );

    // Negative maximum stress
    dimensionedScalar minStress = -maxStress;

    forAll(mesh_.C(), cellI)
    {
        if (alpha1_[cellI] > 0.5)  // Only calculate for solid phase
        {
            scalar strain = thermalStrain[cellI];
            scalar sigmaxx = 
                min
                (
                    max
                    (
                        (2.0 * mu + lambda) * strain,
                        minStress.value()
                    ),
                    maxStress.value()
                );
            stress[cellI] = symmTensor(sigmaxx, 0, 0, sigmaxx, 0, sigmaxx);
        }
    }

    return tStress;
}
void phaseChangeModel::updateDensity()
{
    if (mixture_)
    {
        rho_ = alpha1_*mixture_->rho1() + (1.0 - alpha1_)*mixture_->rho2();
    }
}

void phaseChangeModel::correct(const volScalarField& T)
{
    calculatePhaseChangeRate();
    updateDensity();
    
    if (!checkEnergyConservation())
    {
        WarningIn("phaseChangeModel::correct")
            << "Energy conservation violation detected"
            << endl;
    }
    updateEnergyTracking();
}

void phaseChangeModel::correctPressure(const volScalarField& p)
{
    volScalarField Tm_p = Tm_ + dTmdp_*(p - pRef_);
    
    forAll(mesh_.C(), cellI)
    {
        if (T_[cellI] > Tm_p[cellI])
        {
            phaseIndicator_[cellI] = 1.0;
        }
        else
        {
            phaseIndicator_[cellI] = 0.0;
        }
    }
}

tmp<fvVectorMatrix> phaseChangeModel::momentumSource() const
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    
    tmp<fvVectorMatrix> tMomentumSource
    (
        new fvVectorMatrix
        (
            U,
            dimForce/dimVolume
        )
    );
    
    if (phaseChangeRate_.valid())
    {
        tMomentumSource.ref() += fvm::Sp(phaseChangeRate_(), U);
        
        if (isNonEquilibrium_)
        {
            tMomentumSource.ref() -= fvc::grad(recoilPressure());
        }
    }
    
    return tMomentumSource;
}

tmp<volScalarField> phaseChangeModel::volumetricForce() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "volumetricForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimForce/dimVolume, 0.0)
        )
    );
}

tmp<volScalarField> phaseChangeModel::recoilPressure() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "recoilPressure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    );
}

tmp<volScalarField> phaseChangeModel::source() const
{
    tmp<volScalarField> tSource
    (
        new volScalarField
        (
            IOobject
            (
                "phaseChangeSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& source = tSource.ref();
    
    if (phaseChangeRate_.valid())
    {
        const volScalarField& pcr = phaseChangeRate_();
        
        forAll(mesh_.C(), cellI)
        {
            scalar localDensity = mixture_ ? 
                (alpha1_[cellI]*mixture_->rho1().value() + 
                 (1.0 - alpha1_[cellI])*mixture_->rho2().value()) :
                rho_[cellI];
                
            source[cellI] = pcr[cellI] * localDensity;
        }
    }

    return tSource;
}

tmp<volScalarField> phaseChangeModel::massSource() const
{
    return phaseChangeRate_/Lm_;
}

tmp<fvScalarMatrix> phaseChangeModel::energySource(const volScalarField& h) const
{
    tmp<fvScalarMatrix> tSh
    (
        new fvScalarMatrix
        (
            h,
            dimEnergy/dimVolume/dimTime
        )
    );

    if (phaseChangeRate_.valid())
    {
        const volScalarField& pcr = phaseChangeRate_();
        tSh.ref() -= fvm::Sp(pcr/mixture_->rho1(), h);
        
        scalar L = isNonEquilibrium_ ? (Lm_.value() + Lv_.value()) : Lm_.value();
        tSh.ref() -= pcr * L;
    }

    return tSh;
}

tmp<fvScalarMatrix> phaseChangeModel::electronSource() const
{
    const volScalarField& Te = mesh_.lookupObject<volScalarField>("Te");
    
    tmp<fvScalarMatrix> tSource
    (
        new fvScalarMatrix
        (
            Te,
            dimEnergy/dimVolume/dimTime
        )
    );
    
    if (isNonEquilibrium_ && phaseChangeRate_.valid())
    {
        tSource.ref() -= fvm::Sp(phaseChangeRate_(), Te);
    }
    
    return tSource;
}

tmp<fvScalarMatrix> phaseChangeModel::latticeSource() const
{
    const volScalarField& Tl = mesh_.lookupObject<volScalarField>("Tl");
    
    tmp<fvScalarMatrix> tSource
    (
        new fvScalarMatrix
        (
            Tl,
            dimEnergy/dimVolume/dimTime
        )
    );
    
    if (phaseChangeRate_.valid())
    {
        tSource.ref() -= fvm::Sp(phaseChangeRate_(), Tl);
    }
    
    return tSource;
}

scalar phaseChangeModel::liquidFraction() const
{
    return fvc::domainIntegrate(pos(phaseIndicator_ - 0.5)).value() / 
           gSum(mesh_.V());
}

scalar phaseChangeModel::solidFraction() const
{
    return fvc::domainIntegrate(neg(phaseIndicator_ - 0.5)).value() / 
           gSum(mesh_.V());
}

scalar phaseChangeModel::vaporFraction() const
{
    if (!isNonEquilibrium_) 
        return 0.0;
    
    return fvc::domainIntegrate(pos(T_ - Tv_)).value() / 
           gSum(mesh_.V());
}

scalar phaseChangeModel::interfaceArea() const
{
    return fvc::domainIntegrate(mag(fvc::grad(phaseIndicator_))).value();
}

tmp<volScalarField> phaseChangeModel::subcooling() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject("subcooling", mesh_.time().timeName(), mesh_),
            Tm_ - T_
        )
    );
}

bool phaseChangeModel::valid() const
{
    // Check base parameters
    if (Tm_.value() <= 0 || Lm_.value() <= 0 || interfaceWidth_.value() <= 0)
    {
        return false;
    }

    // Check non-equilibrium parameters if applicable
    if (isNonEquilibrium_)
    {
        if (Tv_.value() <= Tm_.value() || Lv_.value() <= 0)
        {
            return false;
        }
    }

    // Check field existence in registry
    if (!mesh_.foundObject<volScalarField>("T") ||
        !mesh_.foundObject<volScalarField>("alpha.titanium"))
    {
        return false;
    }

 // Check residual stress parameters
    if (E_ <= 0 || nu_ <= -1 || nu_ >= 0.5 || alpha_ <= 0)
    {
        return false;
    }

    // Check for negative or invalid values in fields
    forAll(mesh_.C(), cellI)
    {
        if (T_[cellI] < 0 || alpha1_[cellI] < -SMALL || alpha1_[cellI] > 1.0 + SMALL)
        {
            return false;
        }
    }

    // Check mixture
    if (mixture_)
    {
        // Check for valid density values
        if (mixture_->rho1().value() <= 0 || mixture_->rho2().value() <= 0)
        {
            return false;
        }
    }

    // Check phase change rate field if it exists
    if (phaseChangeRate_.valid())
    {
        const volScalarField& pcr = phaseChangeRate_();
        forAll(mesh_.C(), cellI)
        {
            if (!std::isfinite(pcr[cellI]))
            {
                return false;
            }
        }
    }

    return true;
}

// Add energy conservation check with tolerance
bool phaseChangeModel::checkEnergyConservation() const
{
    if (!energyInitialized_) 
        return true;

    dimensionedScalar currentEnergy = totalPhaseChangeEnergy();
    scalar energyError = mag((currentEnergy.value() - lastTotalEnergy_.value())/
                           (mag(lastTotalEnergy_.value()) + SMALL));

    // Get tolerance from dictionary with safe default
    scalar tolerance = dict_.lookupOrDefault<scalar>("energyTolerance", 1e-6);
    
    if (energyError > tolerance)
    {
        WarningIn("phaseChangeModel::checkEnergyConservation")
            << "Energy conservation error: " << energyError * 100 << "%"
            << endl;
        return false;
    }

    return true;
}
void phaseChangeModel::updateEnergyTracking() const
{
    lastTotalEnergy_ = totalPhaseChangeEnergy();
    energyInitialized_ = true;
}
bool phaseChangeModel::validate() const
{
    // Base validation handled by valid() method
    return valid();
}
dimensionedScalar phaseChangeModel::totalPhaseChangeEnergy() const
{
    return fvc::domainIntegrate
    (
        rho_ * (isNonEquilibrium_ ? 
            (Lm_*pos(phaseIndicator_ - 0.5) + Lv_*pos(T_ - Tv_)) :
            (Lm_*phaseIndicator_))
    );
}

void phaseChangeModel::write() const
{
    Info<< "Phase Change Model:" << nl
        << "  Melting temperature: " << Tm_.value() << " K" << nl
        << "  Vaporization temperature: " << Tv_.value() << " K" << nl
        << "  Phase fractions:" << nl
        << "    Solid: " << solidFraction() << nl
        << "    Liquid: " << liquidFraction();
    
    if (isNonEquilibrium_)
    {
        Info<< nl << "    Vapor: " << vaporFraction();
    }
    
    Info<< nl << "  Energy error: " 
        << mag((totalPhaseChangeEnergy().value() - lastTotalEnergy_.value())/
              (mag(lastTotalEnergy_.value()) + SMALL))
        << nl << "Residual Stress Parameters:" << nl
        << "  Young's modulus: " << E_ << nl
        << "  Poisson's ratio: " << nu_ << nl
        << "  Thermal expansion: " << alpha_ << nl
        << "  Reference temperature: " << Tref_ << endl;

    phaseIndicator_.write();
    rho_.write();
    if (phaseChangeRate_.valid())
    {
        phaseChangeRate_().write();
    }

    // Write residual stress field
    tmp<volSymmTensorField> stress = residualStress();
    stress().write();
}

} // End namespace Foam
