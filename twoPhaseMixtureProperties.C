// twoPhaseMixtureProperties.C
#include "twoPhaseMixtureProperties.H"

namespace Foam
{

twoPhaseMixtureProperties::twoPhaseMixtureProperties(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    ),
    mu_
    (
        IOobject("mu", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("mu", dimDynamicViscosity, 0)
    )
{}

void twoPhaseMixtureProperties::correct(const volScalarField& alpha1)
{
    const scalar rho1 = dict_.get<scalar>("rho1");
    const scalar rho2 = dict_.get<scalar>("rho2");
    const scalar mu1 = dict_.get<scalar>("mu1");
    const scalar mu2 = dict_.get<scalar>("mu2");

    forAll(mesh_.C(), cellI)
    {
        rho_[cellI] = alpha1[cellI] * rho1 + (1 - alpha1[cellI]) * rho2;
        mu_[cellI] = alpha1[cellI] * mu1 + (1 - alpha1[cellI]) * mu2;
    }
}

}
