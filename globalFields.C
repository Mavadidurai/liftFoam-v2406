#include "globalFields.H"

namespace Foam
{
    GlobalFields* globalFieldsPtr = nullptr;

    bool GlobalFields::initialize()
    {
        if (initialized_)
        {
            return true;
        }

        try
        {
            Info<< "Initializing fields" << endl;

            // Initialize velocity field
            UPtr.reset(new volVectorField
            (
                IOobject
                (
                    "U",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            ));

            // Initialize flux field
            phiPtr.reset(new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(U()) & mesh_.Sf()
            ));

            // Initialize pressure fields
            p_rghPtr.reset(new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            ));

            // Initialize mixture
            mixturePtr.reset
            (
                new immiscibleIncompressibleTwoPhaseMixture(U(), phi())
            );

            // Initialize phase fractions
            alpha1Ptr.reset(&mixturePtr->alpha1());
            alpha2Ptr.reset(&mixturePtr->alpha2());

            // Initialize density field
            rhoPtr.reset(new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                alpha1()*mixturePtr->rho1() + alpha2()*mixturePtr->rho2()
            ));

            // Initialize temperature fields
            TePtr.reset(new volScalarField
            (
                IOobject
                (
                    "Te",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            ));

            TlPtr.reset(new volScalarField
            (
                IOobject
                (
                    "Tl",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            ));

            // Initialize turbulence model
            turbulencePtr.reset
            (
                incompressible::turbulenceModel::New
                (
                    U(),
                    phi(),
                    mixturePtr()
                ).ptr()
            );

            initialized_ = true;
            Info<< "Fields initialized successfully" << endl;
            return true;
        }
        catch(Foam::error& e)
        {
            Info<< "Error initializing fields: " << e.message() << endl;
            return false;
        }
    }

    bool GlobalFields::validate()
    {
        if (!initialized_)
        {
            Info<< "Error: Fields not initialized" << endl;
            return false;
        }

        if (!mixturePtr.valid() || !turbulencePtr.valid())
        {
            Info<< "Error: Invalid mixture or turbulence model" << endl;
            return false;
        }

        if (U().size() != mesh_.nCells() || p().size() != mesh_.nCells())
        {
            Info<< "Error: Field sizes do not match mesh" << endl;
            return false;
        }

        return true;
    }

    void GlobalFields::clear()
    {
        UPtr.clear();
        phiPtr.clear();
        pPtr.clear();
        p_rghPtr.clear();
        rhoPtr.clear();
        alpha1Ptr.clear();
        alpha2Ptr.clear();
        TePtr.clear();
        TlPtr.clear();
        mixturePtr.clear();
        turbulencePtr.clear();
        initialized_ = false;
    }
}
