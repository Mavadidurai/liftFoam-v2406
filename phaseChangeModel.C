#include "phaseChangeModel.H"
#include "fvc.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModelBase, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const fvMesh& mesh,
    const volScalarField& T,
    volScalarField& alpha1,
    const dictionary& dict
)
:
    phaseChangeModelBase(),
    mesh_(mesh),
    T_(T),
    alpha1_(alpha1),
    dict_(dict),
    meltingTemperature_
    (
        "meltingTemperature",
        dimTemperature,
        dict.get<scalar>("meltingTemperature")
    ),
    latentHeat_
    (
        "latentHeat",
        dimEnergy/dimMass,
        dict.get<scalar>("latentHeat")
    ),
    undercoolingCoeff_
    (
        "undercoolingCoeff",
        dimTemperature/dimLength,
        dict.get<scalar>("undercoolingCoeff")
    ),
    gradualMeltingRate_
    (
        "gradualMeltingRate",
        dimless/dimTime,
        dict.get<scalar>("gradualMeltingRate")
    ),
    gradualSolidificationRate_
    (
        "gradualSolidificationRate",
        dimless/dimTime,
        dict.get<scalar>("gradualSolidificationRate")
    ),
    interfaceWidth_
    (
        "interfaceWidth",
        dimLength,
        dict.lookupOrDefault<scalar>("interfaceWidth", 5e-6)
    ),
    dTmdp_
    (
        "dTmdp",
        dimTemperature/dimPressure,
        dict.lookupOrDefault<scalar>("dTmdp", 0.0)
    ),
    pRef_
    (
        "pRef",
        dimPressure,
        dict.lookupOrDefault<scalar>("pRef", 1e5)
    ),
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
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1000.0)
    ),
    phaseChangeRate_
    (
        new volScalarField
        (
            IOobject
            (
                "phaseChangeRate",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeModel::correct(const volScalarField& T)
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    volScalarField Tm_eff = meltingTemperature_ + dTmdp_ * (p - pRef_);

    forAll(mesh_.C(), cellI)
    {
        scalar phaseFunction = 0.5 * (1.0 + tanh((T[cellI] - Tm_eff[cellI])/interfaceWidth_.value()));
        phaseIndicator_[cellI] = phaseFunction;

        scalar localDensity = rho_[cellI];
        
        if(phaseChangeRate_.valid())
        {
            phaseChangeRate_.ref()[cellI] = localDensity * latentHeat_.value() * 
                mag(fvc::ddt(phaseIndicator_)().primitiveField()[cellI]);
        }
    }

    dimensionedScalar totalLatentHeat = fvc::domainIntegrate(phaseChangeRate_() * mesh_.time().deltaT());
    Info<< typeName << "::correct(): " 
        << "Total latent heat = " << totalLatentHeat.value() << nl
        << "Phase change: average phase indicator = " 
        << gAverage(phaseIndicator_.primitiveField()) << endl;
}

void Foam::phaseChangeModel::calculatePhaseChangeRate() const
{
    if (!phaseChangeRate_.valid())
    {
        phaseChangeRate_ = tmp<volScalarField>
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
    }

    volScalarField& pcr = phaseChangeRate_.ref();
    tmp<volVectorField> gradT = fvc::grad(T_);
    const volVectorField& gradTRef = gradT();

    forAll(mesh_.C(), cellI)
    {
        scalar Tm_effective = meltingTemperature_.value() - 
                            undercoolingCoeff_.value() * mag(gradTRef[cellI]);
        
        if (T_[cellI] > Tm_effective)
        {
            pcr[cellI] = gradualMeltingRate_.value() * latentHeat_.value();
        }
        else if (T_[cellI] < Tm_effective)
        {
            pcr[cellI] = -gradualSolidificationRate_.value() * latentHeat_.value();
        }
        else
        {
            pcr[cellI] = 0;
        }
    }
}

Foam::tmp<Foam::volScalarField> 
Foam::phaseChangeModel::phaseChangeRate() const
{
    calculatePhaseChangeRate();
    return phaseChangeRate_;
}

void Foam::phaseChangeModel::write() const
{
    Info<< "Phase Change Model parameters:" << nl
        << "  Melting Temperature: " << meltingTemperature_.value() << " K" << nl
        << "  Latent Heat: " << latentHeat_.value() << " J/kg" << nl
        << "  Interface Width: " << interfaceWidth_.value() << " m" << nl
        << "  dT/dp: " << dTmdp_.value() << " K/Pa" << nl
        << "  Reference pressure: " << pRef_.value() << " Pa" << endl;

    if(phaseChangeRate_.valid())
    {
        phaseChangeRate_().write();
    }
    phaseIndicator_.write();
}

// ************************************************************************* //

namespace Foam
{
    addToRunTimeSelectionTable
    (
        phaseChangeModelBase,        // baseType
        phaseChangeModel,            // modelType
        dictionary
    );
}

// ************************************************************************* //
