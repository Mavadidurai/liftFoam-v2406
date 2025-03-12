// correctPhi.H
#ifndef CORRECT_PHI_H
#define CORRECT_PHI_H

{
    #include "continuityErrs.H"

    // Define rAUf if not already defined
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    // Flux predictor using absolute flux from cell center velocities
    phi = mesh.Sf() & fvc::interpolate(U);

    // MRF fluxes
    if (mesh.changing())
    {
        // Relative MRF flux
        phi += fvc::meshPhi(rho, U);
    }

    // Phase change contribution to flux
    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& pcr = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
        if (gMax(mag(pcr)) > SMALL)
        {
            phi += fvc::interpolate(pcr/rho1)*mesh.magSf();
        }
    }

    // Non-conservative fluxes by optional user-defined sources 
    {
        volScalarField pcorr
        (
            IOobject
            (
                "pcorr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0)
        );

        // Solve for correction flux
        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pcorrEqn
            (
                fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
            );

            pcorrEqn.setReference(pRefCell, pRefValue);
            pcorrEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi -= pcorrEqn.flux();
            }
        }
    }

    // Make flux relative to mesh motion
    fvc::makeRelative(phi, U);
}
#endif
