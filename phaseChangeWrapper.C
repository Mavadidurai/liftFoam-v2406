#include "phaseChangeWrapper.H"
#include "phaseChangeModel.H"
#include "nonEquilibriumPhaseChangeModel.H"

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeWrapper, 0);

    phaseChangeWrapper::phaseChangeWrapper
    (
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    )
    :
        mesh_(mesh),
        dict_(dict)
    {
        const word modelType(dict.get<word>("phaseChangeModel"));

        Info<< "Selecting phase change model: " << modelType << endl;

        if (modelType == "nonEquilibrium")
        {
            model_.reset(new nonEquilibriumPhaseChangeModel(mesh, T, alpha1, dict));
        }
        else if (modelType == "equilibrium")
        {
            model_.reset(new phaseChangeModel(mesh, T, alpha1, dict));
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Unknown phaseChangeModel type "
                << modelType << nl << nl
                << "Valid phaseChangeModel types are:" << nl
                << "    equilibrium" << nl
                << "    nonEquilibrium" << exit(FatalIOError);
        }
    }

    void phaseChangeWrapper::correct(const volScalarField& T)
    {
        model_->correct(T);
    }

    tmp<volScalarField> phaseChangeWrapper::phaseChangeRate() const
    {
        return model_->phaseChangeRate();
    }

    tmp<fvScalarMatrix> phaseChangeWrapper::Sh(const volScalarField& h) const
    {
        return model_->Sh(h);
    }

    void phaseChangeWrapper::write() const
    {
        model_->write();
    }
}
