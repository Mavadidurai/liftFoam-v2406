#ifndef femtosecondDropletModel_H
#define femtosecondDropletModel_H

#include "dropletModel.H"
#include "fvCFD.H"

namespace Foam
{
class femtosecondDropletModel : public dropletModel
{
private:
    volVectorField ejectionVelocity_;
    scalar criticalTemperature_;
    scalar criticalPressure_;

public:
    femtosecondDropletModel(const fvMesh& mesh, const dictionary& dict, const volScalarField& rho, const volVectorField& U);

    virtual void update(const volScalarField& T, const volScalarField& p, const volVectorField& U) override;

    tmp<volVectorField> ejectionVelocity() const;
};

}

#endif
