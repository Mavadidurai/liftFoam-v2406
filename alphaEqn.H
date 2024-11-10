#ifndef ALPHAEQN_H
#define ALPHAEQN_H

{
    word alphaScheme("div(phi,alpha)");
    
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

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));
    
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh
        ),
        phi
    );

    // Add compression term
    surfaceScalarField phir(phic*mixture.nHatf());

    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phiAlpha, alpha1)
      + fvm::div(phir, alpha1)
    );

    alpha1Eqn.solve();

    // Ensure boundedness
    alpha1.max(0.0);
    alpha1.min(1.0);
}

#endif
