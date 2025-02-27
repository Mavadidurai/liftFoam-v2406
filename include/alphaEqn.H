#ifndef ALPHAEQN_H
#define ALPHAEQN_H

// alphaEqn.H
{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("interfaceCompression");

    // Store old time alpha1
    volScalarField alpha10("alpha10", alpha1);

    // Calculate interface compression
    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(mixture.cAlpha()*phic, max(phic));
    
    surfaceScalarField phir(phic*mixture.nHatf());

    // Calculate alpha flux
    surfaceScalarField phiAlpha
    (
        "phiAlpha",
        fvc::flux(phi, alpha1)
    );

    // Add compression terms
    phiAlpha += fvc::flux(phir, alpha1*(1.0 - alpha1));

    // Add phase change contribution if enabled
    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& phaseChangeRate = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
        phiAlpha += fvc::interpolate(phaseChangeRate/mixture.rho1());
    }

    // Solve explicit alpha equation using MULES
    MULES::explicitSolve
    (
        geometricOneField(),  // RhoType
        alpha1,              // field
        phi,                 // flux
        phiAlpha,           // alphaFlux  
        zeroField(),        // Sp
        zeroField()         // Su
    );

    // Apply boundedness
    alpha1.correctBoundaryConditions();
    alpha1.max(0.0);
    alpha1.min(1.0);

    // Update mixture properties
    mixture.correct();

    // Update density
    rho = alpha1*mixture.rho1() + (1.0 - alpha1)*mixture.rho2();
    rhoPhi = fvc::interpolate(rho)*phi;
}
#endif
