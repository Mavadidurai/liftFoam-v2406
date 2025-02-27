#ifndef PEQN_H
#define PEQN_H

// pEqn.H
{
    // Face velocity for momentum predictor
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(IOobject::groupName("HbyA", alpha1.group()), U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
    );

    // Update MRF fluxes based on absolute velocity
    MRF.makeRelative(phiHbyA);

    // Add phase change contribution
    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& pcr = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
        phiHbyA += fvc::interpolate(pcr/rho1)*mesh.magSf();
    }

    // Add laser source momentum contribution if significant
    if (mesh.foundObject<volScalarField>("laserSource"))
    {
        const volScalarField& laser = 
            mesh.lookupObject<volScalarField>("laserSource");
        if (max(laser).value() > SMALL)
        {
            phiHbyA += fvc::interpolate(laser/rho1)*mesh.magSf();
        }
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Cache old time-step phi
    surfaceScalarField phiHbyA0("phiHbyA0", phiHbyA);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
          - fvc::ddt(rho) + fvc::ddt(rho.oldTime())
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            // Update flux
            phi = phiHbyA - p_rghEqn.flux();

            // Update velocity
            U = HbyA - rAU*fvc::reconstruct((p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    // Update pressure including hydrostatic contribution
    p = p_rgh + rho*gh;

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    mixture.correct();

    // Update density
    rho = alpha1*rho1 + alpha2*rho2;

    // Correct mass flux
    rhoPhi = fvc::interpolate(rho)*phi;
}

#endif
