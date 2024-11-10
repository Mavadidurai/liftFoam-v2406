// curvatureModel.C
#include "curvatureModel.H"

namespace Foam
{

curvatureModel::curvatureModel(const fvMesh& mesh, const volScalarField& alpha1)
:
    mesh_(mesh),
    alpha1_(alpha1)
{}

tmp<volScalarField> curvatureModel::kappa() const
{
    tmp<volScalarField> tKappa
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0)
        )
    );

    volScalarField& kappa = tKappa.ref();

    volVectorField gradAlpha = fvc::grad(alpha1_);
    volScalarField magGradAlpha = mag(gradAlpha);

    kappa = -fvc::div(gradAlpha/(magGradAlpha + SMALL));

    return tKappa;
}

}
