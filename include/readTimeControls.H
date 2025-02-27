#ifndef readTimeControls_H
#define readTimeControls_H

{
    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    scalar maxCo = 0.0;
    scalar maxAlphaCo = 0.0;
    scalar maxDeltaT = 0.0;

    if (adjustTimeStep)
    {
        maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.5);
        maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", 0.2);
        maxDeltaT = pimpleDict.getOrDefault<scalar>("maxDeltaT", runTime.deltaTValue());
    }

    scalar CoNum = 0.0;
    scalar alphaCoNum = 0.0;

    if (mesh.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        alphaCoNum = 0.5*gMax
        (
            sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
        )*runTime.deltaTValue();
    }

    Info<< "Courant Number mean: " << CoNum
        << " max: " << alphaCoNum << endl;

    // Set deltaTFact
    scalar deltaTFact = 1.0;
    if (adjustTimeStep)
    {
        deltaTFact = min
        (
            min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL)),
            1.2
        );
    }

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );
}

#endif
