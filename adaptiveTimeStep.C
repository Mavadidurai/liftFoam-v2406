#include "adaptiveTimeStep.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "findRefCell.H"
#include "adjustPhi.H"

namespace Foam
{

adaptiveTimeStep::adaptiveTimeStep
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    maxCo_(dict.lookupOrDefault<scalar>("maxCo", 0.5)),
    maxDeltaT_(dict.lookupOrDefault<scalar>("maxDeltaT", GREAT)),
    minDeltaT_(dict.lookupOrDefault<scalar>("minDeltaT", SMALL)),
    cfl_(dict.lookupOrDefault<scalar>("cfl", 1.0)),
    maxTChange_(dict.lookupOrDefault<scalar>("maxTChange", 10.0))
{}

scalar adaptiveTimeStep::calculateDeltaT
(
    const volScalarField& T,
    const volVectorField& U
) const
{
    scalar deltaT = maxDeltaT_;

    if (mesh_.nInternalFaces())
    {
        // Get max velocity magnitude
        tmp<volScalarField> magU = mag(U);
        scalar maxSpeed = gMax(magU().internalField()) + SMALL;

        // Get minimum cell dimension
        scalar minDx = gMin(pow(mesh_.V().field(), 1.0/3.0));

        // Calculate Courant number
        scalar CoNum = 0.5*maxSpeed*mesh_.time().deltaTValue()/minDx;

        if (CoNum > SMALL)
        {
            // Limit based on Courant number
            deltaT = min(maxCo_/CoNum*mesh_.time().deltaTValue(), maxDeltaT_);
        }

        // Temperature change limitation
        if (T.time().timeIndex() > 1)
        {
            volScalarField dTdt = fvc::ddt(T);
            scalar maxTempChange = gMax(mag(dTdt.primitiveField()));
            scalar maxDeltaTTemp = maxTChange_/(maxTempChange + SMALL);
            deltaT = min(deltaT, maxDeltaTTemp);
        }

        // Don't allow the time step to change more than a factor of 2
        deltaT = min(deltaT, 1.2*mesh_.time().deltaTValue());

        // Bound time step between limits
        deltaT = max(min(deltaT, maxDeltaT_), minDeltaT_);
    }

    return deltaT;
}

scalar adaptiveTimeStep::calculateDeltaT
(
    const volScalarField& T,
    const volVectorField& U,
    const dimensionedScalar& maxDeltaT
) const
{
    // Get laser parameters
    const scalar laserPulseDuration = 
        dict_.lookupOrDefault<scalar>("laserPulseDuration", 214e-15);
    const scalar electronRelaxationTime = 
        dict_.lookupOrDefault<scalar>("electronRelaxationTime", 1e-14);
    
    scalar deltaT = maxDeltaT.value();

    if (mesh_.nInternalFaces())
    {
        // Calculate flow-based time step (similar to base function)
        tmp<volScalarField> magU = mag(U);
        scalar maxSpeed = gMax(magU().internalField()) + SMALL;
        scalar minDx = gMin(pow(mesh_.V().field(), 1.0/3.0));
        scalar CoNum = 0.5*maxSpeed*mesh_.time().deltaTValue()/minDx;
        scalar dtFlow = (CoNum > SMALL) ? maxCo_/CoNum*mesh_.time().deltaTValue() : maxDeltaT.value();

        // Calculate physics-based time steps
        scalar dtLaser = laserPulseDuration/10;  // Resolve laser pulse
        scalar dtElectron = electronRelaxationTime/10;  // Resolve electron dynamics
        
        // Take minimum of all time scales
        deltaT = min(dtLaser, dtElectron);
        deltaT = min(deltaT, dtFlow);

        // Temperature change limitation
        if (T.time().timeIndex() > 1)
        {
            volScalarField dTdt = fvc::ddt(T);
            scalar maxTempChange = gMax(mag(dTdt.primitiveField()));
            scalar maxDeltaTTemp = maxTChange_/(maxTempChange + SMALL);
            deltaT = min(deltaT, maxDeltaTTemp);
        }

        // Don't allow the time step to change more than a factor of 2
        deltaT = min(deltaT, 1.2*mesh_.time().deltaTValue());

        // Additional limits
        deltaT = min(deltaT, maxDeltaT.value());
        deltaT = max(deltaT, minDeltaT_);

        Info<< "Time scales:" << nl
            << "  Laser dt: " << dtLaser << nl
            << "  Electron dt: " << dtElectron << nl
            << "  Flow dt: " << dtFlow << nl
            << "  Final dt: " << deltaT << endl;
    }

    return deltaT;
}

} // End namespace Foam
