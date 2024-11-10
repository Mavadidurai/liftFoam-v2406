#include "multiScaleVisualization.H"

namespace Foam
{

multiScaleVisualization::multiScaleVisualization(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    macroScale_
    (
        IOobject("macroScale", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("macroScale", dimless, 0)
    ),
    mesoScale_
    (
        IOobject("mesoScale", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("mesoScale", dimless, 0)
    ),
    microScale_
    (
        IOobject("microScale", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("microScale", dimless, 0)
    )
{}

void multiScaleVisualization::update(const volScalarField& T, const volScalarField& alpha1)
{
    // Update multi-scale fields based on temperature and phase fraction
    forAll(mesh_.C(), cellI)
    {
        macroScale_[cellI] = T[cellI];
        mesoScale_[cellI] = mag(fvc::grad(T)()[cellI]);
        microScale_[cellI] = mag(fvc::grad(alpha1)()[cellI]);
    }
}

void multiScaleVisualization::write()
{
    macroScale_.write();
    mesoScale_.write();
    microScale_.write();
}

}
