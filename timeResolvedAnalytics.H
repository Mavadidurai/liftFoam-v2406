#ifndef timeResolvedAnalytics_H
#define timeResolvedAnalytics_H

#include "fvCFD.H"

namespace Foam
{

class timeResolvedAnalytics
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;
    
    DynamicList<scalar> timeHistory_;
    DynamicList<scalar> temperatureHistory_;
    DynamicList<scalar> pressureHistory_;

public:
    timeResolvedAnalytics(const fvMesh& mesh, const dictionary& dict);

    void calculate(const volScalarField& T, const volScalarField& p);
    void write();
};

}

#endif
