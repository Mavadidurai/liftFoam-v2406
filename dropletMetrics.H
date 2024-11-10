// dropletMetrics.H
#ifndef dropletMetrics_H
#define dropletMetrics_H

#include "fvCFD.H"

namespace Foam
{

class dropletMetrics
{
private:
    const fvMesh& mesh_;
    const volScalarField& alpha1_;

public:
    dropletMetrics(const fvMesh& mesh, const volScalarField& alpha1);

    scalar aspectRatio() const;
    scalar circularity() const;
};

}

#endif
