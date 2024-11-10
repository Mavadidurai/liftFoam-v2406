// curvatureModel.H
#ifndef curvatureModel_H
#define curvatureModel_H

#include "fvCFD.H"

namespace Foam
{

class curvatureModel
{
private:
    const fvMesh& mesh_;
    const volScalarField& alpha1_;

public:
    curvatureModel(const fvMesh& mesh, const volScalarField& alpha1);

    tmp<volScalarField> kappa() const;
};

}

#endif
