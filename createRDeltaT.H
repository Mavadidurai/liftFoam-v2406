{
    volScalarField& rDeltaT = rDeltaT.ref();
// Define maxDeltaT if it's not already defined
scalar maxDeltaT = pimpleDict.get<scalar>("maxDeltaT");

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(pimpleDict.get<scalar>("maxCo"));

    // Maximum time scale ratio
    scalar maxAlphaCo(pimpleDict.get<scalar>("maxAlphaCo"));

    // Minimum time step
    scalar minDeltaT(pimpleDict.get<scalar>("minDeltaT"));

    // Read the maximum damping factor
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Read the maximum reciprocal time step
    dimensionedScalar rDeltaTMax
    (
        ((pimpleDict.found("rDeltaTMax"))
      ? dimensionedScalar("rDeltaTMax", dimless/dimTime, pimpleDict)
      : dimensionedScalar("rDeltaTMax", dimless/dimTime, GREAT))
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi))()()
       /((2*maxCo)*mesh.V()*rho())
    );

    // Further limit the reciprocal time-step
    // in the vicinity of the interface
    rDeltaT.ref() = max
    (
        rDeltaT(),
        pos0(alpha1 - 0.01)*pos0(0.99 - alpha1)
       *fvc::surfaceSum(mag(phi*mixture.nHatf()))()()
       /((2*maxAlphaCo)*mesh.V()*rho())
    );

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTSmoothingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (1.0 - rDeltaTSmoothingCoeff)*rDeltaT.oldTime()
        );
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Bound the reciprocal time-step
    rDeltaT.ref() = max(rDeltaT(), 1/dimensionedScalar("minDeltaT", dimTime, minDeltaT));
    rDeltaT.ref() = min(rDeltaT(), rDeltaTMax);
}
