#ifndef setDeltaT_H
#define setDeltaT_H

// setDeltaT.H

// Read from the PIMPLE dictionary
const dictionary& pimpleDict = pimple.dict();

// Get time step adjuster
bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

// Get maximum Courant number and alpha Courant number
scalar maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.5);
scalar maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", 0.25);
scalar maxDeltaT = pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT);

// Calculate Courant numbers
scalar CoNum = 0.0;
scalar alphaCoNum = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    // Velocity Courant number
    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    // Interface Courant number
    alphaCoNum = 0.5*gMax
    (
        sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
    )*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << CoNum
    << " max: " << alphaCoNum << endl;

// Set the new deltaT based on Courant numbers
if (adjustTimeStep)
{
    scalar maxDeltaTFact = 
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

    // Limit the time step change
    maxDeltaTFact = min(maxDeltaTFact, 1.2);

    runTime.setDeltaT
    (
        min
        (
            maxDeltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

#endif
