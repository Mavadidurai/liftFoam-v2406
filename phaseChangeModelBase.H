#ifndef phaseChangeModelBase_H
#define phaseChangeModelBase_H

#include "fvCFD.H"
#include "fvMatrix.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phaseChangeModelBase
{
protected:
    //- Constructor
    phaseChangeModelBase() = default;

public:
    //- Runtime type information
    TypeName("phaseChangeModelBase");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModelBase,
        dictionary,
        (
            const fvMesh& mesh,
            const volScalarField& T,
            volScalarField& alpha1,
            const dictionary& dict
        ),
        (mesh, T, alpha1, dict)
    );

    //- Selector
    static autoPtr<phaseChangeModelBase> New
    (
        const word& modelType,
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    );

    //- Destructor
    virtual ~phaseChangeModelBase() = default;

    // Member Functions
    virtual void correct(const volScalarField& T) = 0;
    virtual tmp<fvScalarMatrix> Sh(const volScalarField& h) const = 0;
    virtual void write() const = 0;
    virtual tmp<volScalarField> phaseChangeRate() const = 0;
};

} // End namespace Foam

#endif
