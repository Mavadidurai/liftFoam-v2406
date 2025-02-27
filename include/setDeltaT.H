#ifndef setDeltaT_H
#define setDeltaT_H

{
    // Ensure phi exists before using it
    if (!mesh.foundObject<surfaceScalarField>("phi"))
    {
        FatalErrorIn("setDeltaT.H")
            << "phi field not found in registry before Courant number calculation!"
            << abort(FatalError);
    }

    // Calculate Courant number properly using phi
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar alphaCoNum = 0.0;  // Added declaration
    
    if (mesh.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        CoNum = 0.5 * gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();
        meanCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field())) * runTime.deltaTValue();

        // Calculate interface Courant number
        alphaCoNum = 0.5 * gMax
        (
            sumPhi /
            (
                mesh.V().field()
              * mesh.surfaceInterpolation::deltaCoeffs().primitiveField()
            )
        ) * runTime.deltaTValue();
    }

    // Get PIMPLE controls
    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    
    // Read or use defaults
    const scalar maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.5);
    const scalar maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", 0.2);
    const scalar maxDeltaT = pimpleDict.getOrDefault<scalar>("maxDeltaT", 1e-6); // Increased from 1e-12 to 1e-6

    // Calculate time factors with better stability
    scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
    scalar alphaCoDeltaTFact = maxAlphaCo / (alphaCoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, alphaCoDeltaTFact), 1.2);

    // Set new deltaT with bounds checking
    runTime.setDeltaT
    (
        min
        (
            deltaTFact * runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum 
        << " interface Courant Number: " << alphaCoNum << endl;
}

#endif
