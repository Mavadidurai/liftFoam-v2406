#ifndef ultraFastShockWaveModel_H
#define ultraFastShockWaveModel_H

#include "fvCFD.H"

namespace Foam
{

class ultraFastShockWaveModel
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    volScalarField shockStrength_;
    scalar shockThreshold_;

public:
    ultraFastShockWaveModel(const fvMesh& mesh, const dictionary& dict);
bool valid() const
    {
        return shockStrength_.size() > 0;
    }
    void update(const volScalarField& p, const volVectorField& U);
    tmp<volScalarField> shockStrength() const;
    void propagate();
    tmp<volVectorField> velocityContribution() const;
};

}

#endif
