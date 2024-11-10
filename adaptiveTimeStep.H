// adaptiveTimeStep.H
#ifndef adaptiveTimeStep_H
#define adaptiveTimeStep_H

#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedScalar.H"

namespace Foam
{

class adaptiveTimeStep
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;
    
    scalar maxCo_;
    scalar maxDeltaT_;
    scalar minDeltaT_;
    scalar cfl_;
    scalar maxTChange_;

public:
    // Constructor
    adaptiveTimeStep(const fvMesh& mesh, const dictionary& dict);

    // Check if configuration is valid
    bool valid() const
    {
        return maxCo_ > 0 && maxDeltaT_ > minDeltaT_ && minDeltaT_ > 0;
    }

    // Calculate adaptive time step - with dimensionedScalar
    scalar calculateDeltaT
    (
        const volScalarField& T,
        const volVectorField& U,
        const dimensionedScalar& maxDeltaT
    ) const;

    // Calculate adaptive time step - without dimensionedScalar
    scalar calculateDeltaT
    (
        const volScalarField& T,
        const volVectorField& U
    ) const;

    // Access functions
    scalar maxCo() const { return maxCo_; }
};

}

#endif
