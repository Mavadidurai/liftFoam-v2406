#ifndef UEQN_H
#define UEQN_H

{
    // Add validity check for phi
    if (!mesh.foundObject<surfaceScalarField>("phi"))
    {
        FatalErrorIn("UEqn.H")
            << "phi field not found - must be initialized before UEqn"
            << exit(FatalError);
    }

    // Ensure rho is positive before usage
    if (min(rho).value() <= 0)
    {
        FatalErrorIn("UEqn.H")
            << "Non-positive density detected in rho!"
            << abort(FatalError);
    }

    // Calculate effective viscosity with null check for turbulence
    volScalarField muEff("muEff", mixture.mu());
    if (turbulence.valid())
    {
        muEff += rho * turbulence->mut();
    }

    // Handle mesh motion
    if (mesh.moving())
    {
        MRF.correctBoundaryVelocity(U);
    }

    // Store momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf * fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        ));

        // Check for NaN or divergence in U (component-wise check)
        if (gMin(U.component(vector::X)().primitiveField()) == -GREAT || gMax(U.component(vector::X)().primitiveField()) == GREAT ||
            gMin(U.component(vector::Y)().primitiveField()) == -GREAT || gMax(U.component(vector::Y)().primitiveField()) == GREAT ||
            gMin(U.component(vector::Z)().primitiveField()) == -GREAT || gMax(U.component(vector::Z)().primitiveField()) == GREAT)
        {
            FatalErrorIn("UEqn.H")
                << "Velocity field has diverged (NaN or Inf values detected in one of the components)!"
                << abort(FatalError);
        }

        fvOptions.correct(U);
    }
}

#endif

