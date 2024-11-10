#include "femtosecondLaserModel.H"

namespace Foam
{

femtosecondLaserModel::femtosecondLaserModel(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    peakIntensity_(dict.get<scalar>("peakIntensity")),
    pulseWidth_(dict.get<scalar>("pulseWidth")),
    wavelength_(dict.get<scalar>("wavelength")),
    direction_(dict.get<vector>("direction")),
    focus_(dict.get<point>("focus"))
{}

void femtosecondLaserModel::update()
{
    // Update laser parameters if needed
}

tmp<volScalarField> femtosecondLaserModel::source() const
{
    tmp<volScalarField> tSource
    (
        new volScalarField
        (
            IOobject
            (
                "laserSource",
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
    const scalar t = mesh_.time().value();

    forAll(mesh_.C(), cellI)
    {
        const point& cellCenter = mesh_.C()[cellI];
        vector r = cellCenter - focus_;
        scalar z = (r & direction_);
        r -= z*direction_;

        scalar R = mag(r);
        scalar temporalTerm = Foam::exp(-4.0*Foam::log(2.0)*sqr(t/pulseWidth_));
        scalar spatialTerm = Foam::exp(-2.0*sqr(R)/(wavelength_*z/Foam::constant::mathematical::pi));

        source[cellI] = peakIntensity_*temporalTerm*spatialTerm;
    }

    return tSource;
}
}
