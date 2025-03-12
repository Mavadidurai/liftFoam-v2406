#ifndef alphaEqnSubCycle_H
#define alphaEqnSubCycle_H
// alphaEqnSubCycle.H

{
    // Apply subcycling to improve stability for alpha equation
    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", rhoPhi.dimensions(), 0)
        );

        // Loop through sub-cycles
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            // Solve alpha equation for this subcycle
            #include "alphaEqn.H"
            
            // Accumulate weighted rhoPhi for the timestep
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        // Correct the flux using the accumulated values
        rhoPhi = rhoPhiSum;
    }
    else
    {
        // No subcycling, just solve the alpha equation once
        #include "alphaEqn.H"
    }
}

#endif
