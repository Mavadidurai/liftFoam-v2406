// radiationModel.H
#ifndef radiationModel_H
#define radiationModel_H

#include "fvCFD.H"
#include "fvMesh.H"

namespace Foam
{

class radiationModel
{
protected:
    const fvMesh& mesh_;
    dictionary dict_;
    mutable volScalarField absorptivity_;  // Make mutable since we modify in const method
    mutable volScalarField emissivity_;    // Make mutable since we modify in const method

public:
    radiationModel(const fvMesh& mesh, const dictionary& dict);
    virtual ~radiationModel() {}

    virtual void correct(const volScalarField& T);
    virtual tmp<volScalarField> Ru() const;
    virtual void updateProperties() const;  // Keep const
    virtual void write() const;
};

}

#endif
