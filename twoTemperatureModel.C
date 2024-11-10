#include "twoTemperatureModel.H"

namespace Foam
{

twoTemperatureModel::twoTemperatureModel(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    Te_
    (
        IOobject("Te", mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    ),
    Tl_
    (
        IOobject("Tl", mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    ),
    Ce_(dimensionedScalar("Ce", dimEnergy/dimVolume/dimTemperature, dict.get<scalar>("Ce"))),
    Cl_(dimensionedScalar("Cl", dimEnergy/dimVolume/dimTemperature, dict.get<scalar>("Cl"))),
    G_(dimensionedScalar("G", dimEnergy/dimVolume/dimTime/dimTemperature, dict.get<scalar>("G")))
{
    // Verify field dimensions
    if (Te_.dimensions() != dimTemperature || Tl_.dimensions() != dimTemperature)
    {
        FatalErrorIn("twoTemperatureModel::twoTemperatureModel")
            << "Incorrect temperature dimensions" << nl
            << "Te dimensions: " << Te_.dimensions() << nl
            << "Tl dimensions: " << Tl_.dimensions()
            << abort(FatalError);
    }
}

void twoTemperatureModel::solve(const volScalarField& laserSource)
{
    dimensionedScalar smallTemp("smallTemp", dimTemperature, SMALL);
    dimensionedScalar maxTemp("maxTemp", dimTemperature, 1e5);

    // Calculate temperature-dependent properties
    volScalarField ke = electronThermalConductivity();
    volScalarField Ce = electronHeatCapacity();
    volScalarField G = electronPhononCoupling();

    // Calculate energy before solving
    dimensionedScalar initElectronEnergy = 
        fvc::domainIntegrate(Ce * Te_);
    dimensionedScalar initLatticeEnergy = 
        fvc::domainIntegrate(Cl_ * Tl_);

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
    Te_ = max(min(Te_, maxTemp), smallTemp);
    Tl_ = max(min(Tl_, maxTemp), smallTemp);

    // Check energy conservation
    dimensionedScalar finalElectronEnergy = 
        fvc::domainIntegrate(Ce * Te_);
    dimensionedScalar finalLatticeEnergy = 
        fvc::domainIntegrate(Cl_ * Tl_);
    dimensionedScalar energyDifference = 
        (finalElectronEnergy + finalLatticeEnergy) -
        (initElectronEnergy + initLatticeEnergy);

    if (mag(energyDifference).value() > SMALL)
    {
        Info<< "Warning: Energy conservation error = " 
            << energyDifference.value() << endl;
    }
}

tmp<volScalarField> twoTemperatureModel::electronThermalConductivity() const
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
            mesh_,
            dimensionedScalar("ke", dimPower/dimLength/dimTemperature, 1.0)
        )
    );
}

tmp<volScalarField> twoTemperatureModel::electronHeatCapacity() const
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
            mesh_,
            Ce_
        )
    );
}

tmp<volScalarField> twoTemperatureModel::electronPhononCoupling() const
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
            mesh_,
            G_
        )
    );
}

} // End namespace Foam
