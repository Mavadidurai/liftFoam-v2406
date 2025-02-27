#ifndef initializeFields_H 
#define initializeFields_H

{
    // Initialize mixture 
    mixture.correct();

    // Initialize density and flux fields
    rho = alpha1*mixture.rho1() + (1.0 - alpha1)*mixture.rho2();
    rhoPhi = fvc::interpolate(rho)*phi;

    // Initialize temperature fields
    Te.correctBoundaryConditions();
    Tl.correctBoundaryConditions();
    
    // Initialize phase change and laser fields
    phaseChangeRate.correctBoundaryConditions();
    laserSource.correctBoundaryConditions();

    // Time controls 
    scalar maxCo = pimple.dict().lookupOrDefault("maxCo", 0.5);
    scalar maxAlphaCo = pimple.dict().lookupOrDefault("maxAlphaCo", 0.2);

    // Apply the time controls
    if (mesh.nInternalFaces())
    {
        scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
        scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
        scalar alphaCoNum = 0.5*gMax
        (
            sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
        )*runTime.deltaTValue();

        if (CoNum > maxCo || alphaCoNum > maxAlphaCo)
        {
            scalar deltaTFact = min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));
            runTime.setDeltaT(deltaTFact*runTime.deltaTValue());
        }
    }
}

#endif
