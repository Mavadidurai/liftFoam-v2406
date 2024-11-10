#ifndef ALPHAEQNSUBCYCLE_H
#define ALPHAEQNSUBCYCLE_H

{
    #include "alphaControls.H"
    
    // Store initial mass
    const dimensionedScalar totalDeltaT = runTime.deltaT();
    const volScalarField alpha10("alpha10", alpha1);
    const dimensionedScalar mass0 = fvc::domainIntegrate(alpha1);

    // Initialize rhoPhiSum for sub-cycling
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
        dimensionedScalar(rhoPhi.dimensions(), 0)
    );

    // Calculate interface compression magnitudes
    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));
    
    // Get non-const reference to mixture for interface calculations
    immiscibleIncompressibleTwoPhaseMixture& mixtureRef =
        const_cast<immiscibleIncompressibleTwoPhaseMixture&>(mixture);
    
    surfaceScalarField phir(phic*mixtureRef.nHatf());

    // Create limiter fields for MULES
    volScalarField limiter1
    (
        IOobject
        (
            "limiter1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    volScalarField limiter2
    (
        IOobject
        (
            "limiter2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Alpha correction loop
    for (label aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        // Create source terms for phase fraction equation
        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        // Calculate alpha flux
        surfaceScalarField phiAlpha
        (
            IOobject
            (
                "phiAlpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi*fvc::interpolate(alpha1)
        );

        // Add compression flux
        surfaceScalarField phiAlphaCompr
        (
            IOobject
            (
                "phiAlphaCompr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiAlpha + fvc::flux(-phir, scalar(1) - alpha1, alpharScheme)
        );

        // Solve using MULES
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            phiAlphaCompr,
            Su,
            Sp,
            limiter1,
            limiter2
        );

        // Apply explicit bounds after MULES
        alpha1.max(0);
        alpha1.min(1);
    }

    // Sub-cycle the phase fraction equation
    if (nAlphaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }

    // Check mass conservation
    const dimensionedScalar mass = fvc::domainIntegrate(alpha1);
    const dimensionedScalar massError = 
        (mass - mass0)/max(mass0, dimensionedScalar("small", dimless, SMALL));

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << "  Mass error = " << massError.value()
        << endl;

    // Update mixture properties using non-const reference
    mixtureRef.correct();
}

#endif
