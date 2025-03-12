#ifndef PEQN_H
#define PEQN_H
// pEqn.H

{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    
    volVectorField HbyA(U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
    );

    // Apply MRF correction
    MRF.makeRelative(phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Add phase change contribution
    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& pcr = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
            
        if (max(mag(pcr)).value() > SMALL)
        {
            dimensionedScalar rho1Inv("rho1Inv", dimVolume/dimMass, 1.0/mixture.rho1().value());
            phiHbyA += fvc::interpolate(pcr*rho1Inv)*mesh.magSf();
        }
    }

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
            
            U = HbyA - rAU*fvc::reconstruct(p_rghEqn.flux()/rAUf);
            U.correctBoundaryConditions();
        }
    }

    // Update pressure
    p = p_rgh + rho*gh;

    // Update mass flux
    rhoPhi = phi*fvc::interpolate(rho);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
#endif
