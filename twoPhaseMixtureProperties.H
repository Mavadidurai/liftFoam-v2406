// twoPhaseMixtureProperties.H
#ifndef twoPhaseMixtureProperties_H
#define twoPhaseMixtureProperties_H

#include "fvCFD.H"
#include "fvMesh.H"

namespace Foam
{

class twoPhaseMixtureProperties
{
protected:
    const fvMesh& mesh_;
    dictionary dict_;
    volScalarField rho_;
    volScalarField mu_;

public:
    twoPhaseMixtureProperties(const fvMesh& mesh, const dictionary& dict);
 bool valid() const
    {
        return rho_.size() > 0 && mu_.size() > 0;
    }
    void correct(const volScalarField& alpha1);
    const volScalarField& rho() const { return rho_; }
    const volScalarField& mu() const { return mu_; }
};
}

#endif

