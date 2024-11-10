#ifndef dropletModel_H
#define dropletModel_H

#include "fvCFD.H"
#include "volFields.H"

namespace Foam
{

class dropletModel
{
protected:
    const fvMesh& mesh_;
    const volScalarField& rho_;
    const volVectorField& U_;
    scalar L_;
    scalar criticalWe_;
    scalar viscosity_;
    scalar surfaceTension_;
    volScalarField dropletIndicator_;

public:
    dropletModel(const fvMesh& mesh, const dictionary& dict, const volScalarField& rho, const volVectorField& U);

    virtual ~dropletModel() = default;  // Add virtual destructor

    virtual void update(const volScalarField& T, const volScalarField& p, const volVectorField& U);
    bool valid() const
    {
        return dropletIndicator_.size() > 0 && surfaceTension_ > 0 && viscosity_ > 0;
    }
    void checkBreakup();
    void checkCoalescence();
    bool isDropletFormed() const;
    void calculateTrajectory();
};

}

#endif
