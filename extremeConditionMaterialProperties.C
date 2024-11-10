#include "extremeConditionMaterialProperties.H"

namespace Foam
{

extremeConditionMaterialProperties::extremeConditionMaterialProperties(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    rho0_(dict.get<scalar>("rho0")),
    Cp0_(dict.get<scalar>("Cp0")),
    k0_(dict.get<scalar>("k0")),
    Te0_(dict.get<scalar>("Te0")),
    A_(dict.get<scalar>("A")),
    B_(dict.get<scalar>("B"))
{}

void extremeConditionMaterialProperties::update(const volScalarField& T)
{
    // Update material properties if needed
}

tmp<volScalarField> extremeConditionMaterialProperties::rho(const volScalarField& T, const volScalarField& p) const
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
            rho0_ * (1 - A_ * (T - 300))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::Cp(const volScalarField& T, const volScalarField& p) const
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
            Cp0_ * (1 + B_ * (T - 300))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::k(const volScalarField& T, const volScalarField& p) const
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
            k0_ * (1 + 0.01 * (T - 300))
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::electronThermalConductivity(const volScalarField& Te) const
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
            k0_ * pow(Te / Te0_, 1.5)
        )
    );
}

tmp<volScalarField> extremeConditionMaterialProperties::electronHeatCapacity(const volScalarField& Te) const
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
            A_ * Te
        )
    );
}

}
