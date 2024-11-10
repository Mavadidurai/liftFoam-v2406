#ifndef setDeltaT_H
#define setDeltaT_H

{
    // Calculate Courant number factors
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar alphaDeltaTFact = maxAlphaCo/(alphaCoNum + SMALL);

    // Calculate initial time step factor
    scalar deltaTFact = min(min(maxDeltaTFact, alphaDeltaTFact), 1.0);

    // Get time scales from dictionary
    const dictionary& pimpleDict = pimple.dict();
    const scalar laserPulseDuration = 
        pimpleDict.lookupOrDefault<scalar>("laserPulseDuration", 214e-15);
    const scalar electronRelaxationTime = 
        pimpleDict.lookupOrDefault<scalar>("electronRelaxationTime", 1e-14);

    // Calculate physics-based time step
    const scalar physicsTimeStep = min
    (
        laserPulseDuration/10.0,
        electronRelaxationTime/10.0
    );

    // Limit rate of change of time step
    const scalar maxChange = 1.2;
    deltaTFact = min(deltaTFact, maxChange);

    // Calculate new time step considering all constraints
    scalar deltaT = min
    (
        deltaTFact*runTime.deltaTValue(),
        maxDeltaT
    );

    // Apply physics-based limit
    deltaT = min(deltaT, physicsTimeStep);

    // Apply absolute minimum limit
    deltaT = max(deltaT, minDeltaT);

    // Set the new time step
    runTime.setDeltaT(deltaT);

    Info<< "Time step adjustment:" << nl
        << "    Courant Number max: " << CoNum << nl
        << "    Interface Courant Number max: " << alphaCoNum << nl
        << "    Physics time scale: " << physicsTimeStep << nl
        << "    New deltaT = " << runTime.deltaTValue() << nl
        << endl;
}

#endif
