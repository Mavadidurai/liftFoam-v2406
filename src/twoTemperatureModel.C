/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of the two-temperature model for femtosecond laser-material
    interaction in LIFT process. Handles:
    - Temperature field initialization and evolution
    - Material property calculations
    - Energy conservation tracking
    - Coupled electron-lattice temperature solution
    
    The model uses temperature-dependent material properties and includes
    electron-phonon coupling for accurate simulation of ultrafast laser heating.

\*---------------------------------------------------------------------------*/

#include "twoTemperatureModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    twoTemperatureModel Implementation
\*---------------------------------------------------------------------------*/

twoTemperatureModel::twoTemperatureModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    Te_
    (
        IOobject
        (
            "Te",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Te", dimTemperature, 300.0)
    ),
    Tl_
    (
        IOobject
        (
            "Tl",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tl", dimTemperature, 300.0)
    ),
    Ce_
    (
        "Ce",
        dimEnergy/dimVolume/dimTemperature,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "Ce",
            dimensionedScalar
            (
                "Ce_default",
                dimEnergy/dimVolume/dimTemperature,
                210.0
            )
        ).value()
    ),
    Cl_
    (
        "Cl",
        dimEnergy/dimVolume/dimTemperature,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "Cl",
            dimensionedScalar
            (
                "Cl_default",
                dimEnergy/dimVolume/dimTemperature,
                2.3e6
            )
        ).value()
    ),
    G_
    (
        "G",
        dimEnergy/dimVolume/dimTime/dimTemperature,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "G",
            dimensionedScalar
            (
                "G_default",
                dimEnergy/dimVolume/dimTime/dimTemperature,
                2.6e17
            )
        ).value()
    ),
    lastTotalEnergy_
    (
        "lastTotalEnergy",
        dimEnergy,
        0.0
    ),
    energyInitialized_(false)
{
    if (!validateParameters())
    {
        FatalErrorInFunction
            << "Invalid model parameters"
            << abort(FatalError);
    }

    // Register fields if not already present
    if (!mesh_.foundObject<volScalarField>("Te"))
    {
        mesh_.objectRegistry::store(new volScalarField(Te_));
    }
    if (!mesh_.foundObject<volScalarField>("Tl"))
    {
        mesh_.objectRegistry::store(new volScalarField(Tl_));
    }

    // Initialize energy tracking
    updateEnergyTracking();
}

twoTemperatureModel::~twoTemperatureModel()
{
    if (mesh_.foundObject<volScalarField>("Te"))
    {
        mesh_.objectRegistry::checkOut("Te");
    }
    if (mesh_.foundObject<volScalarField>("Tl"))
    {
        mesh_.objectRegistry::checkOut("Tl");
    }
}

bool twoTemperatureModel::validateParameters() const
{
    bool valid = true;

    // Check dimensions
    if (Te_.dimensions() != dimTemperature || 
        Tl_.dimensions() != dimTemperature)
    {
        FatalErrorInFunction
            << "Invalid temperature dimensions"
            << abort(FatalError);
        valid = false;
    }

    if (Ce_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        Cl_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        G_.dimensions() != dimEnergy/dimVolume/dimTime/dimTemperature)
    {
        FatalErrorInFunction
            << "Invalid material property dimensions"
            << abort(FatalError);
        valid = false;
    }

    // Check property values
    if (Ce_.value() <= 0 || Cl_.value() <= 0 || G_.value() <= 0)
    {
        FatalErrorInFunction
            << "Non-positive material properties detected"
            << abort(FatalError);
        valid = false;
    }

    // Check temperature values
    forAll(Te_, cellI)
    {
        if (Te_[cellI] < 0 || Tl_[cellI] < 0)
        {
            FatalErrorInFunction
                << "Negative temperature detected at cell " << cellI
                << abort(FatalError);
            valid = false;
            break;
        }
    }

    return valid;
}

bool twoTemperatureModel::validateFields() const
{
    bool valid = true;

    forAll(Te_, cellI)
    {
        if (Te_[cellI] < 0 || !std::isfinite(Te_[cellI]) ||
            Tl_[cellI] < 0 || !std::isfinite(Tl_[cellI]))
        {
            valid = false;
            break;
        }
    }

    return valid;
}

bool twoTemperatureModel::checkEnergyConservation() const
{
    if (!energyInitialized_)
    {
        return true;
    }

    dimensionedScalar currentEnergy = 
        fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);

    scalar energyError = mag
    (
        (currentEnergy.value() - lastTotalEnergy_.value())/
        (mag(lastTotalEnergy_.value()) + SMALL)
    );

    return energyError < dict_.lookupOrDefault<scalar>("energyTolerance", 1e-6);
}

void twoTemperatureModel::updateEnergyTracking() const
{
    lastTotalEnergy_ = fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);
    energyInitialized_ = true;
}


void twoTemperatureModel::solve(const volScalarField& laserSource)
{
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid field values before solve"
            << abort(FatalError);
    }

    // Store initial energy
    updateEnergyTracking();

    // Calculate temperature-dependent properties
    volScalarField ke = electronThermalConductivity();
    volScalarField Ce = electronHeatCapacity();
    volScalarField G = electronPhononCoupling();

    // Solve electron temperature equation
    fvScalarMatrix TeEqn
    (
        fvm::ddt(Ce, Te_)
      - fvm::laplacian(ke, Te_)
     ==
        laserSource
      - G*(Te_ - Tl_)
    );

    TeEqn.relax();
    TeEqn.solve();

    // Solve lattice temperature equation
    fvScalarMatrix TlEqn
    (
        fvm::ddt(Cl_, Tl_)
     ==
        G*(Te_ - Tl_)
    );

    TlEqn.relax();
    TlEqn.solve();

    // Apply temperature bounds
    dimensionedScalar minTemp("minTemp", dimTemperature, SMALL);
    dimensionedScalar maxTemp("maxTemp", dimTemperature, 1e5);

    Te_ = max(min(Te_, maxTemp), minTemp);
    Tl_ = max(min(Tl_, maxTemp), minTemp);

    // Check energy conservation
    if (!checkEnergyConservation())
    {
        WarningInFunction
            << "Energy conservation violation detected" << nl
            << "Error = " << mag((lastTotalEnergy_.value()
                               - fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_).value())
                               /(mag(lastTotalEnergy_.value()) + SMALL))
            << endl;
    }

    // Update energy tracking
    updateEnergyTracking();

    // Validate final state
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid field values after solve"
            << abort(FatalError);
    }

    // Report solution statistics
    Info<< "Two-temperature solve:" << nl
        << "  Te range: " << min(Te_).value() << " - " << max(Te_).value() << " K" << nl
        << "  Tl range: " << min(Tl_).value() << " - " << max(Tl_).value() << " K" << nl
        << "  Max temperature difference: " << max(mag(Te_ - Tl_)).value() << " K" << endl;
}

tmp<volScalarField> twoTemperatureModel::electronThermalConductivity() const
{
    // Create temporary field for electron thermal conductivity
    tmp<volScalarField> tke
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
            mesh_,
            dimensionedScalar("ke", dimPower/dimLength/dimTemperature, 1.0)
        )
    );

    // Get reference to field for modification
    volScalarField& ke = tke.ref();

    // Calculate temperature-dependent conductivity
    // This is a simplified model - could be extended for more complex behavior
    forAll(ke, cellI)
    {
        ke[cellI] = Ce_.value() * Te_[cellI] / (3.0 * G_.value());
    }

    return tke;
}

tmp<volScalarField> twoTemperatureModel::electronHeatCapacity() const
{
    // Create temporary field for electron heat capacity
    tmp<volScalarField> tCe
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
            mesh_,
            Ce_
        )
    );

    return tCe;
}

tmp<volScalarField> twoTemperatureModel::electronPhononCoupling() const
{
    // Create temporary field for electron-phonon coupling
    tmp<volScalarField> tG
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
            mesh_,
            G_
        )
    );

    return tG;
}

bool twoTemperatureModel::valid() const
{
    if (!validateParameters())
    {
        return false;
    }

    if (!validateFields())
    {
        return false;
    }

    // Check field dimensions
    if (Te_.dimensions() != dimTemperature ||
        Tl_.dimensions() != dimTemperature ||
        Ce_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        Cl_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        G_.dimensions() != dimEnergy/dimVolume/dimTime/dimTemperature)
    {
        return false;
    }

    return true;
}

void twoTemperatureModel::write() const
{
    Info<< "Two-temperature model:" << nl
        << "Parameters:" << nl
        << "  Ce = " << Ce_.value() << " J/m³/K" << nl
        << "  Cl = " << Cl_.value() << " J/m³/K" << nl
        << "  G = " << G_.value() << " W/m³/K" << nl
        << "Field statistics:" << nl
        << "  Te range: " << min(Te_).value() << " - " << max(Te_).value() << " K" << nl
        << "  Tl range: " << min(Tl_).value() << " - " << max(Tl_).value() << " K" << nl
        << "  Mean Te: " << average(Te_).value() << " K" << nl
        << "  Mean Tl: " << average(Tl_).value() << " K" << nl;

    if (energyInitialized_)
    {
        dimensionedScalar currentEnergy = fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);
        scalar energyError = mag
        (
            (currentEnergy.value() - lastTotalEnergy_.value())/
            (mag(lastTotalEnergy_.value()) + SMALL)
        );
        
        Info<< "Energy conservation:" << nl
            << "  Current total energy: " << currentEnergy.value() << " J" << nl
            << "  Energy error: " << energyError * 100 << " %" << endl;
    }
}
// Add these implementations after the existing functions in twoTemperatureModel.C

tmp<volScalarField> Foam::twoTemperatureModel::ke() const
{
    // Return electronic thermal conductivity
    tmp<volScalarField> tke
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
            mesh_,
            dimensionedScalar("ke", dimPower/dimLength/dimTemperature, 0.0)
        )
    );

    volScalarField& ke = tke.ref();

    // Calculate temperature-dependent electronic thermal conductivity
    // Using Wiedemann-Franz law and electron temperature dependence
    forAll(mesh_.C(), cellI)
    {
        scalar Te = Te_[cellI];
        ke[cellI] = Ce_.value() * Te / (3.0 * G_.value());
        
        // Apply temperature dependent correction
        if (Te > 1000.0)
        {
            ke[cellI] *= pow(Te/1000.0, 1.5);  // High temperature correction
        }
    }

    return tke;
}

tmp<volScalarField> Foam::twoTemperatureModel::kl() const
{
    // Return lattice thermal conductivity
    tmp<volScalarField> tkl
    (
        new volScalarField
        (
            IOobject
            (
                "kl",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("kl", dimPower/dimLength/dimTemperature, 0.0)
        )
    );

    volScalarField& kl = tkl.ref();

    // Calculate temperature-dependent lattice thermal conductivity
    // Using Drude model with phonon contribution
    forAll(mesh_.C(), cellI)
    {
        scalar Tl = Tl_[cellI];
        kl[cellI] = Cl_.value() * sqrt(Tl) / (3.0 * G_.value());
        
        // Apply temperature dependent correction
        if (Tl > 1000.0)
        {
            kl[cellI] *= pow(1000.0/Tl, 0.5);  // High temperature correction
        }
    }

    return tkl;
}
} // End namespace Foam
