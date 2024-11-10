#ifndef alphaSolve_H
#define alphaSolve_H

inline void solveAlpha
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiAlpha,
    const dimensionedScalar& deltaT
)
{
    // Create bound fields
    scalarField alphaMin(alpha1.size(), 0.0);
    scalarField alphaMax(alpha1.size(), 1.0);

    // Solve using MULES with proper bound fields
    MULES::explicitSolve
    (
        geometricOneField(),  // One field for constant density
        alpha1,              // Phase fraction field to solve
        phi,                 // Flux field
        phiAlpha,           // Limited phase flux
        alphaMin,           // Lower bound field (0)
        alphaMax,           // Upper bound field (1)
        1,                  // Number of corrections
        0                   // Debug level
    );

    // Apply explicit bounds after solution
    alpha1.boundaryFieldRef().evaluate();
    
    // Ensure strict bounds
    alpha1.max(0);
    alpha1.min(1);

    Info<< "Phase-1 volume fraction = " 
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}

#endif
