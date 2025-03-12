#ifndef initializeFields_H 
#define initializeFields_H

{
    Info<< "Initializing fields..." << endl;
    
    // Initialize mixture properties
    mixture.correct();

    // Initialize density and flux fields
    rho = alpha1*mixture.rho1() + (1.0 - alpha1)*mixture.rho2();
    rhoPhi = fvc::interpolate(rho)*phi;

    // Initialize temperature fields
    T.correctBoundaryConditions();
    
    // Initialize femtosecond-specific fields if in femtosecond regime
    if (regime == "femtosecond")
    {
        // Ensure Te and Tl exist
        if (mesh.foundObject<volScalarField>("Te"))
        {
            volScalarField& Te = const_cast<volScalarField&>(
                mesh.lookupObject<volScalarField>("Te"));
            Te.correctBoundaryConditions();
        }
        
        if (mesh.foundObject<volScalarField>("Tl"))
        {
            volScalarField& Tl = const_cast<volScalarField&>(
                mesh.lookupObject<volScalarField>("Tl"));
            Tl.correctBoundaryConditions();
        }
        
        if (mesh.foundObject<volScalarField>("laserSource"))
        {
            volScalarField& laserSource = const_cast<volScalarField&>(
                mesh.lookupObject<volScalarField>("laserSource"));
            laserSource.correctBoundaryConditions();
        }
        
        if (mesh.foundObject<volScalarField>("phaseChangeRate"))
        {
            volScalarField& phaseChangeRate = const_cast<volScalarField&>(
                mesh.lookupObject<volScalarField>("phaseChangeRate"));
            phaseChangeRate.correctBoundaryConditions();
        }
    }

    // Initialize pressure fields
    p = p_rgh + rho*gh;
    p.correctBoundaryConditions();
    p_rgh.correctBoundaryConditions();

    // Ensure all fields are properly registered in mesh
    if (!mesh.foundObject<surfaceScalarField>("phi"))
    {
        mesh.objectRegistry::store(new surfaceScalarField(phi));
    }

    // Initialize flux field
    phi = linearInterpolate(U) & mesh.Sf();
    phi.correctBoundaryConditions();

    // Initialize mass flux
    rhoPhi = fvc::interpolate(rho)*phi;

    Info<< "Fields initialized" << endl;
}

#endif
