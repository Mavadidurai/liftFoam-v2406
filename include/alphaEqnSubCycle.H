#ifndef alphaEqnSubCycle_H
#define alphaEqnSubCycle_H

{
   const label nAlphaSubCycles = pimple.dict().lookupOrDefault("nAlphaSubCycles", 1);

    if (nAlphaSubCycles > 1)
    {
        // Store original phi
        surfaceScalarField phiCopy = phi;
        
        dimensionedScalar totalDeltaT = runTime.deltaT();

        // Initialize sum field with proper dimensions
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE  // Changed from default
            ),
            mesh,
            dimensionedScalar("zero", rhoPhi.dimensions(), 0.0)
        );

        // Correct subcycling loop syntax
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        while (!(++alphaSubCycle).end())  // Corrected loop structure
        {
            // Update phi for subcycle
            phi = phiCopy * (runTime.deltaT()/totalDeltaT);
            
            #include "alphaEqn.H"

            rhoPhiSum += (runTime.deltaT()/totalDeltaT) * rhoPhi;
        }

        // Restore original phi
phi = phiCopy;
mesh.objectRegistry::store(new surfaceScalarField(phi));

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }
}

#endif

