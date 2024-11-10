#ifndef UEqn_H
#define UEqn_H

{
    Info<< "Solving momentum equation" << endl;

    // Calculate effective viscosity
    volScalarField muEff(turbulencePtr->muEff());
    
    // Create momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
    );

    UEqn.relax();

    // MRF contribution (use correct API)
    {
        volScalarField rAU(1.0/UEqn.A());
        
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
        
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
        );

        MRF.makeRelative(phiHbyA);

        adjustPhi(phiHbyA, U, p_rgh);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, U, phiHbyA, rAU);
    }

    // Solve momentum predictor
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
}

#endif
