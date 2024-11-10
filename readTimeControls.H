#ifndef readTimeControls_H
#define readTimeControls_H

{
    // Read from controlDict
    adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    // General time controls from controlDict
    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    maxAlphaCo = runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 1.0);
    maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    minDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

    // Read PIMPLE controls
    const dictionary& pimpleDict = pimple.dict();
    maxCo = pimpleDict.lookupOrDefault<scalar>("maxCo", maxCo);
    maxAlphaCo = pimpleDict.lookupOrDefault<scalar>("maxAlphaCo", maxAlphaCo);
    maxDeltaT = pimpleDict.lookupOrDefault<scalar>("maxDeltaT", maxDeltaT);

    // LIFT-specific time controls
    const scalar laserPulseDuration = 
        pimpleDict.lookupOrDefault<scalar>("laserPulseDuration", 214e-15);
    const scalar electronRelaxationTime = 
        pimpleDict.lookupOrDefault<scalar>("electronRelaxationTime", 1e-14);

    // Calculate physics-based time step
    const scalar physicsTimeScale = min
    (
        laserPulseDuration/10.0,
        electronRelaxationTime/10.0
    );

    // Update minimum time step considering all constraints
    minDeltaT = max
    (
        minDeltaT,
        physicsTimeScale
    );

    // Verify time step limits
    if (maxDeltaT < minDeltaT)
    {
        FatalErrorInFunction
            << "Invalid time step limits:" << nl
            << "maxDeltaT = " << maxDeltaT << nl
            << "minDeltaT = " << minDeltaT << nl
            << "physicsTimeScale = " << physicsTimeScale << nl
            << "maxDeltaT must be greater than minDeltaT"
            << abort(FatalError);
    }

    // Output current settings
    Info<< "Time control settings:" << nl
        << "  adjustTimeStep        : " << adjustTimeStep << nl
        << "  maxCo                : " << maxCo << nl
        << "  maxAlphaCo           : " << maxAlphaCo << nl
        << "  maxDeltaT            : " << maxDeltaT << nl
        << "  minDeltaT            : " << minDeltaT << nl
        << "  laserPulseDuration   : " << laserPulseDuration << nl
        << "  electronRelaxTime    : " << electronRelaxationTime << nl
        << "  physicsTimeScale     : " << physicsTimeScale << nl
        << endl;
}

#endif
