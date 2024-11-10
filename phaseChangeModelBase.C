#include "phaseChangeModelBase.H"

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModelBase, 0);
    defineRunTimeSelectionTable(phaseChangeModelBase, dictionary);

    autoPtr<phaseChangeModelBase> phaseChangeModelBase::New
    (
        const word& modelType,
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    )
    {
        Info<< "Selecting phase change model: " << modelType << endl;

        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "phaseChangeModel",
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<phaseChangeModelBase>
        (
            ctorPtr(mesh, T, alpha1, dict)
        );
    }
}
