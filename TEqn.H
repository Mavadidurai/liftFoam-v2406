#ifndef TEqn_H
#define TEqn_H

{
    // Store initial temperature fields for energy conservation check
    volScalarField Te_old = Te;
    volScalarField Tl_old = Tl;

    // Get material properties from dictionary
    const dictionary& materialDict = mesh.lookupObject<IOdictionary>("materialProperties");
    
    // Electronic specific heat
    dimensionedScalar Ce
    (
        "Ce",
        dimEnergy/dimVolume/dimTemperature,
        materialDict.get<scalar>("electronSpecificHeat")
    );
    
    // Lattice specific heat
    dimensionedScalar Cl
    (
        "Cl",
        dimEnergy/dimVolume/dimTemperature,
        materialDict.get<scalar>("latticeSpecificHeat")
    );
    
    // Electron-phonon coupling coefficient
    dimensionedScalar G
    (
        "G",
        dimEnergy/dimVolume/dimTime/dimTemperature,
        materialDict.get<scalar>("electronPhononCoupling")
    );

    // Calculate initial energy
    dimensionedScalar initialEnergy = fvc::domainIntegrate
    (
        Ce*Te + Cl*Tl + 0.5*rho*magSqr(U)
    );

    // Electron temperature equation
    fvScalarMatrix TeEqn
    (
        fvm::ddt(Ce, Te)
      - fvm::laplacian(thermo.alphahe(), Te)
     ==
        laser.source()
      - G*(Te - Tl)
      + fvOptions(rho, Te)
    );

    TeEqn.relax();

    // Solve electron temperature equation and check convergence
    scalar TeResidual = TeEqn.solve().initialResidual();
    
    if (TeResidual > mesh.solution().dict().get<scalar>("temperatureTolerance"))
    {
        Info<< "Warning: High electron temperature equation residual: " 
            << TeResidual << endl;
    }

    // Lattice temperature equation
    fvScalarMatrix TlEqn
    (
        fvm::ddt(Cl, Tl)
      + fvm::div(phi, Tl)
      - fvm::laplacian(thermo.alphac(), Tl)
     ==
        G*(Te - Tl)
      + phaseChange.source()
      + radiation.source()
      + fvOptions(rho, Tl)
    );

    TlEqn.relax();

    // Solve lattice temperature equation and check convergence
    scalar TlResidual = TlEqn.solve().initialResidual();
    
    if (TlResidual > mesh.solution().dict().get<scalar>("temperatureTolerance"))
    {
        Info<< "Warning: High lattice temperature equation residual: " 
            << TlResidual << endl;
    }

    // Apply temperature bounds
    dimensionedScalar minTemp("minTemp", dimTemperature, 300.0);
    dimensionedScalar maxTemp("maxTemp", dimTemperature, 1e5);
    
    Te = max(min(Te, maxTemp), minTemp);
    Tl = max(min(Tl, maxTemp), minTemp);

    // Calculate final energy and check conservation
    dimensionedScalar finalEnergy = fvc::domainIntegrate
    (
        Ce*Te + Cl*Tl + 0.5*rho*magSqr(U)
    );
    
    scalar energyError = mag
    (
        (finalEnergy.value() - initialEnergy.value()) /
        (initialEnergy.value() + SMALL)
    );
    
    if (energyError > mesh.solution().dict().get<scalar>("energyTolerance"))
    {
        Info<< "Warning: Energy conservation error in temperature solution: " 
            << energyError << endl;
    }

    // Calculate electron-lattice coupling strength
    volScalarField couplingStrength = mag(G*(Te - Tl));
    scalar maxCoupling = max(couplingStrength).value();
    
    if (maxCoupling > mesh.solution().dict().get<scalar>("maxCouplingStrength"))
    {
        Info<< "Warning: Strong electron-lattice coupling detected: " 
            << maxCoupling << endl;
    }

    // Report temperature statistics
    Info<< "Temperature ranges:" << nl
        << "  Electron temperature: " 
        << min(Te).value() << " - " << max(Te).value() << " K" << nl
        << "  Lattice temperature: " 
        << min(Tl).value() << " - " << max(Tl).value() << " K" << nl
        << "  Temperature difference: " 
        << max(mag(Te - Tl)).value() << " K" << nl
        << "  Energy conservation error: " << energyError << endl;

    // Update thermo properties
    thermo.correct();
}

#endif
