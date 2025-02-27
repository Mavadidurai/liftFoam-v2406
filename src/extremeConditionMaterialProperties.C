#include "extremeConditionMaterialProperties.H"

namespace Foam
{

extremeConditionMaterialProperties::extremeConditionMaterialProperties
(
    const fvMesh& mesh,
    const dictionary& dict  
)
:
    mesh_(mesh),
    dict_(dict),
    
    // Initialize base properties for titanium
    rho0_("rho0", dimDensity, dict.get<scalar>("rho0")),
    Cp0_("Cp0", dimEnergy/dimMass/dimTemperature, dict.get<scalar>("Cp0")), 
    k0_("k0", dimPower/dimLength/dimTemperature, dict.get<scalar>("k0")),
    Tm_("Tm", dimTemperature, 1941.0),  // Titanium melting point
    Tv_("Tv", dimTemperature, 3560.0),  // Titanium boiling point
    Lm_("Lm", dimEnergy/dimMass, 3.65e5),  // Titanium latent heat of melting
    Lv_("Lv", dimEnergy/dimMass, 8.85e6),  // Titanium latent heat of vaporization
    sigma_("sigma", dimForce/dimLength, 1.64),  // Titanium surface tension
    
    // Temperature coefficients
    alphaT_(dict.getOrDefault<scalar>("alphaT", 8.6e-6)),  // Thermal expansion
    betaT_(dict.getOrDefault<scalar>("betaT", 0.02)),      // Conductivity temp coeff
    gammaT_(dict.getOrDefault<scalar>("gammaT", 0.01)),    // Heat capacity temp coeff
    
    // Electron-phonon coupling
    G0_(dict.getOrDefault<scalar>("G0", 2.6e17)),    // Coupling coefficient
    Ce0_(dict.getOrDefault<scalar>("Ce0", 215.0)),   // Electron heat capacity coeff
    Cl0_(dict.getOrDefault<scalar>("Cl0", 2.4e6))    // Lattice heat capacity coeff
{}

bool extremeConditionMaterialProperties::valid() const
{
    bool valid = true;

    // Check base properties
    valid &= rho0_.value() > 0;
    valid &= Cp0_.value() > 0;
    valid &= k0_.value() > 0;
    valid &= Tm_.value() > 0;

    // Check phase change properties
    valid &= Tv_.value() > Tm_.value();
    valid &= Lm_.value() > 0;
    valid &= Lv_.value() > 0;
    valid &= sigma_.value() > 0;

    // Check coefficients
    valid &= alphaT_ >= 0;
    valid &= betaT_ >= 0;
    valid &= gammaT_ >= 0;
    valid &= G0_ > 0;
    valid &= Ce0_ > 0;
    valid &= Cl0_ > 0;

    return valid;
}

void extremeConditionMaterialProperties::update(const volScalarField& T)
{
    forAll(mesh_.C(), cellI)
    {
        scalar dT = T[cellI] - 300.0;  // Temperature difference from reference
        k0_.value() *= (1.0 + betaT_*dT);
        Cp0_.value() *= (1.0 + gammaT_*dT);
    }
}

tmp<volScalarField> extremeConditionMaterialProperties::rho
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho0_ * (1.0 - alphaT_*(T - 300.0))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::Cp
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Cp0_ * (1.0 + gammaT_*(T - 300.0))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::k
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            k0_ * (1.0 + betaT_*(T - 300.0))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::electronThermalConductivity
(
    const volScalarField& Te
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ke",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            k0_ * pow(Te/300.0, 1.5) * (1.0 + betaT_*sqr(Te/300.0))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::electronHeatCapacity
(
    const volScalarField& Te
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Ce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Ce0_ * Te/300.0
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::electronPhononCoupling
(
    const volScalarField& Te,
    const volScalarField& Tl
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "G",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            G0_ * (Te/300.0) * (1.0 - exp(-(Te - Tl)/300.0))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::surfaceTension
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sigma_ * (1.0 - (T - Tm_)/(Tv_ - Tm_))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::phaseChangeRate
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tphaseChangeRate
    (
        new volScalarField
        (
            IOobject
            (
                "phaseChangeRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& pcr = tphaseChangeRate.ref();

    // Calculate phase change rate
    forAll(mesh_.C(), cellI)
    {
        if (T[cellI] > Tm_.value()) 
        {
            pcr[cellI] = Lm_.value() * rho0_.value() * 
                        (T[cellI] - Tm_.value())/(mesh_.time().deltaTValue());
        }
        if (T[cellI] > Tv_.value())
        {
            pcr[cellI] += Lv_.value() * rho0_.value() * 
                         (T[cellI] - Tv_.value())/(mesh_.time().deltaTValue());
        }
    }

    return tphaseChangeRate;
}

} // End namespace Foam
