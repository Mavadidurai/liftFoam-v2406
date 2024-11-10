#ifndef phaseChangeWrapper_H
#define phaseChangeWrapper_H

#include "fvCFD.H"
#include "phaseChangeModelBase.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phaseChangeWrapper
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;
    autoPtr<phaseChangeModelBase> model_;

    //- Disallow copy construct
    phaseChangeWrapper(const phaseChangeWrapper&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const phaseChangeWrapper&) = delete;

public:
    //- Runtime type information
    TypeName("phaseChangeWrapper");

    // Constructors

        //- Construct from components
        phaseChangeWrapper
        (
            const fvMesh& mesh,
            const volScalarField& T,
            volScalarField& alpha1,
            const dictionary& dict
        );

    //- Destructor
    virtual ~phaseChangeWrapper() = default;

    // Member Functions

    virtual void correct(const volScalarField& T);
    virtual tmp<volScalarField> phaseChangeRate() const;
    virtual tmp<fvScalarMatrix> Sh(const volScalarField& h) const;
    virtual void write() const;
};

}

#endif

