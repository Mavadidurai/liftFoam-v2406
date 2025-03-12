#ifndef UEQN_H
#define UEQN_H

// UEqn.H

// Solve momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
);

// Add turbulence contribution
turbulence->correct();

// Add effective viscosity contribution from turbulence model
volScalarField muEff(mixture.mu());
if (turbulence->nut().valid())
{
    muEff = muEff + rho*turbulence->nut();
}

UEqn.relax();

// Apply MRF rotation correction if present
MRF.correctBoundaryVelocity(U);

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
      + fvOptions(rho, U)
    );

    fvOptions.correct(U);
}
#endif

