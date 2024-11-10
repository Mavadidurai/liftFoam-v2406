// radiationModel.C
#include "radiationModel.H"
#include <cmath> // Optional if pow4 not available

namespace Foam
{

void radiationModel::updateProperties() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    forAll(mesh_.C(), cellI)
    {
        absorptivity_[cellI] *= (1 + 0.001 * (T[cellI] - 300));
        emissivity_[cellI] *= (1 + 0.001 * (T[cellI] - 300));
    }
}

void radiationModel::write() const
{
    Info<< "Radiation Model:" << nl
        << "  Average absorptivity: " << gAverage(absorptivity_) << nl
        << "  Average emissivity: " << gAverage(emissivity_) << endl;
}

radiationModel::radiationModel(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    absorptivity_
    (
        IOobject("absorptivity", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("absorptivity", dimless, dict.lookupOrDefault<scalar>("absorptivity", 0.5))
    ),
    emissivity_
    (
        IOobject("emissivity", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("emissivity", dimless, dict.lookupOrDefault<scalar>("emissivity", 0.5))
    )
{}

void radiationModel::correct(const volScalarField& T)
{
    updateProperties();
}

tmp<volScalarField> radiationModel::Ru() const
{
    if (!mesh_.foundObject<volScalarField>("T"))
    {
        FatalErrorIn("radiationModel::Ru")
            << "Temperature field 'T' not found in the mesh."
            << exit(FatalError);
    }
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    tmp<volScalarField> tRu
    (
        new volScalarField
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPower/dimArea, 0)
        )
    );

    volScalarField& Ru = tRu.ref();
    
    forAll(mesh_.C(), cellI)
    {
        Ru[cellI] = emissivity_[cellI] * 5.67e-8 * std::pow(T[cellI], 4)
                  - absorptivity_[cellI] * 5.67e-8 * std::pow(300.0, 4);
    }

    return tRu;
}

}



