#include "nonEquilibriumPhaseChangeModel.H"

namespace Foam
{

nonEquilibriumPhaseChangeModel::nonEquilibriumPhaseChangeModel
(
    const fvMesh& mesh,
    const volScalarField& T,
    volScalarField& alpha1,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    Tm_("Tm", dimTemperature, dict.get<scalar>("Tm")),
    Tv_("Tv", dimTemperature, dict.get<scalar>("Tv")),
    Lm_("Lm", dimEnergy/dimMass, dict.get<scalar>("Lm")),
    Lv_("Lv", dimEnergy/dimMass, dict.get<scalar>("Lv")),
    phaseIndicator_
    (
        IOobject
        (
            "phaseIndicator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    )
{
    // Verify dimensions
    if (T.dimensions() != dimTemperature)
    {
        FatalErrorIn("nonEquilibriumPhaseChangeModel::nonEquilibriumPhaseChangeModel")
            << "Incorrect temperature dimensions: " << T.dimensions()
            << abort(FatalError);
    }

    // Verify temperature ordering
    if (Tm_.value() >= Tv_.value())
    {
        FatalErrorIn("nonEquilibriumPhaseChangeModel::nonEquilibriumPhaseChangeModel")
            << "Invalid temperature values: Melting temperature must be less than vaporization temperature\n"
            << "Tm = " << Tm_.value() << " K, Tv = " << Tv_.value() << " K"
            << abort(FatalError);
    }
}

void nonEquilibriumPhaseChangeModel::correct(const volScalarField& T)
{
    if (T.dimensions() != dimTemperature)
    {
        FatalErrorIn("nonEquilibriumPhaseChangeModel::correct")
            << "Incorrect temperature dimensions"
            << abort(FatalError);
    }

    forAll(mesh_.C(), cellI)
    {
        if (T[cellI] < Tm_.value())
        {
            phaseIndicator_[cellI] = 0;  // Solid
        }
        else if (T[cellI] >= Tm_.value() && T[cellI] < Tv_.value())
        {
            phaseIndicator_[cellI] = 1;  // Liquid
        }
        else
        {
            phaseIndicator_[cellI] = 2;  // Vapor
        }
    }
}

tmp<volScalarField> nonEquilibriumPhaseChangeModel::phaseChangeRate() const
{
    dimensionedScalar deltaT("deltaT", dimTime, mesh_.time().deltaTValue());

    tmp<volScalarField> tRate
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
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    volScalarField& rate = tRate.ref();

    forAll(mesh_.C(), cellI)
    {
        if (phaseIndicator_[cellI] == 1)  // Melting
        {
            rate[cellI] = (Lm_/deltaT).value();
        }
        else if (phaseIndicator_[cellI] == 2)  // Vaporization
        {
            rate[cellI] = (Lv_/deltaT).value();
        }
    }

    return tRate;
}

tmp<fvScalarMatrix> nonEquilibriumPhaseChangeModel::Sh(const volScalarField& h) const
{
    if (h.dimensions() != dimEnergy/dimMass)
    {
        FatalErrorIn("nonEquilibriumPhaseChangeModel::Sh")
            << "Incorrect enthalpy dimensions: " << h.dimensions()
            << abort(FatalError);
    }

    return tmp<fvScalarMatrix>();
}

void nonEquilibriumPhaseChangeModel::write() const
{
    Info<< "Non-Equilibrium Phase Change Model parameters:" << nl
        << "  Melting Temperature: " << Tm_.value() << " K" << nl
        << "  Vaporization Temperature: " << Tv_.value() << " K" << nl
        << "  Latent Heat of Melting: " << Lm_.value() << " J/kg" << nl
        << "  Latent Heat of Vaporization: " << Lv_.value() << " J/kg" << endl;
}

}
