// createTimeControls.H
{
    // Read time control parameters from PIMPLE dictionary
    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    scalar maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.5);
    scalar maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", 0.2);
    scalar maxDeltaT = pimpleDict.getOrDefault<scalar>("maxDeltaT", runTime.deltaTValue());

    // Initialize Courant numbers
    scalar CoNum = 0.0;
    scalar alphaCoNum = 0.0;

    if (mesh.nInternalFaces())
    {
        const surfaceScalarField& phi = 
            mesh.lookupObject<surfaceScalarField>("phi");

        // Calculate face flux magnitudes
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        // Calculate Courant number for velocity field
        CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        // Calculate Courant number at the interface
        alphaCoNum = 0.5*gMax
        (
            sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
        )*runTime.deltaTValue();
    }

    // Report Courant numbers
    Info<< "Courant Number mean: " << CoNum
        << " max: " << alphaCoNum << endl;

    // Adjust time step if required
    if (adjustTimeStep)
    {
        // Calculate time step factor
        scalar deltaTFact = min
        (
            min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL)),
            1.2  // Limit expansion rate
        );

        // Set new time step
        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );

        Info<< "New timestep: " << runTime.deltaTValue() << " s" << endl;
    }
}
