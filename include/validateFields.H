#ifndef validateFields_H
#define validateFields_H

{
    Info<< "Validating fields..." << endl;

    // Check phase fraction bounds
    if (min(alpha1).value() < -SMALL || max(alpha1).value() > 1.0 + SMALL)
    {
        Info<< "Warning: Phase fraction out of bounds: " << min(alpha1).value()
            << " to " << max(alpha1).value() << endl;
        
        // Apply limits
        alpha1.max(0.0);
        alpha1.min(1.0);
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
    }

    // Check temperature positivity
    if (min(T).value() < 0)
    {
        Info<< "Warning: Negative temperature detected: " << min(T).value() << endl;
        
        // Apply limits
        T.max(dimensionedScalar("minT", dimTemperature, 1.0));
        T.correctBoundaryConditions();
    }

    // Check femtosecond temperature fields if in femtosecond regime
    if (regime == "femtosecond")
    {
        if (min(Te).value() < 0 || min(Tl).value() < 0)
        {
            Info<< "Warning: Negative femtosecond temperatures detected: "
                << "Te_min = " << min(Te).value() << ", "
                << "Tl_min = " << min(Tl).value() << endl;
            
            // Apply limits
            Te.max(dimensionedScalar("minTe", dimTemperature, 1.0));
            Tl.max(dimensionedScalar("minTl", dimTemperature, 1.0));
            Te.correctBoundaryConditions();
            Tl.correctBoundaryConditions();
        }
    }

    // Check density positivity
    if (min(rho).value() <= 0)
    {
        Info<< "Warning: Non-positive density detected: " << min(rho).value() << endl;
        
        // Update density
        rho = max(alpha1*rho1 + alpha2*rho2, 
                  dimensionedScalar("minRho", dimDensity, SMALL));
        rho.correctBoundaryConditions();
    }

    Info<< "Fields validated" << endl;
}
#endif
