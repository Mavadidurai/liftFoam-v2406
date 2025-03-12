/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
        Solver for Laser-Induced Forward Transfer (LIFT) process.
        Handles both regular and femtosecond regimes.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// Include model headers 
#include "twoTemperatureModel.H"
#include "phaseChangeModel.H"
#include "femtosecondLaserModel.H"
#include "dropletModel.H"

// Forward declaration for createMandatoryFields
namespace Foam
{
    void createMandatoryFields(fvMesh& mesh);
    bool fieldExists(const fvMesh& mesh, const word& fieldName);
}

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for Laser-Induced Forward Transfer (LIFT) process."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    
    // Create cellLevel and pointLevel fields for dynamic mesh refinement
    Info << "Creating cellLevel and pointLevel fields for dynamic mesh refinement" << endl;
    
    // Essential fields that must exist before mesh operations
    if (!fieldExists(mesh, "cellLevel"))
    {
        auto cellLevelPtr = new volScalarField
        (
            IOobject
            (
                "cellLevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );
        mesh.objectRegistry::store(cellLevelPtr);
    }

    if (!fieldExists(mesh, "pointLevel")) 
    {
        auto pointLevelPtr = new volScalarField
        (
            IOobject
            (
                "pointLevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );
        mesh.objectRegistry::store(pointLevelPtr);
    }

    #include "createFields.H"
    
    // Read regime from dictionary
    const word regime(mesh.solutionDict().lookupOrDefault<word>("regime", "regular"));
    const bool isFemtosecond = (regime == "femtosecond");

    Info<< "Running in " << regime << " regime\n" << endl;

    // Set pressure reference
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p_rgh.name());

    #include "createTimeControls.H"
    
    // Calculate initial Courant Number
    scalar initialCoNum = 0.0;
    scalar initialAlphaCoNum = 0.0;
    scalar meanCoNum = 0.0;
    
    if (mesh.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        initialCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
        meanCoNum = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
        
        // Calculate initial interface Courant number
        initialAlphaCoNum = 0.5*gMax
        (
            sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
        )*runTime.deltaTValue();
    }
    
    Info<< "Initial Courant Number = " << initialCoNum << endl;
    Info<< "Initial interface Courant Number = " << initialAlphaCoNum << endl;
    Info<< "Initial mean Courant Number = " << meanCoNum << endl;
    
    // Initialization of important variables before using setInitialDeltaT
    bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
    scalar maxCo = 0.5;
    scalar maxAlphaCo = 0.2;
    scalar maxDeltaT = runTime.endTime().value();
    
    if (adjustTimeStep)
    {
        maxCo = pimple.dict().getOrDefault<scalar>("maxCo", 0.5);
        maxAlphaCo = pimple.dict().getOrDefault<scalar>("maxAlphaCo", 0.2);
        maxDeltaT = pimple.dict().getOrDefault<scalar>("maxDeltaT", maxDeltaT);
    }
    
    // Custom implementation of setInitialDeltaT without using the include file
    if (adjustTimeStep)
    {
        if (runTime.timeIndex() == 0 && initialCoNum > SMALL)
        {
            scalar maxDeltaTFact = min(maxCo/(initialCoNum + SMALL), 
                                       maxAlphaCo/(initialAlphaCoNum + SMALL));
            scalar deltaTFact = min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact);
            runTime.setDeltaT
            (
                min
                (
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            );
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }
    }
    
    // Include initialization and validation
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
    
    // Field validation steps
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

    // Set up transport properties dictionary
    const dictionary& dict = 
        mesh.lookupObject<IOdictionary>("transportProperties");

    // Initialize models
    autoPtr<phaseChangeModel> phaseChangePtr(nullptr);
    if (isFemtosecond)
    {
        Info<< "Initializing femtosecond phase change model" << endl;
        const word modelType = "femtosecondPhaseChange";
        phaseChangePtr = phaseChangeModel::New(modelType, mesh, T, alpha1, dict);
    }
    else
    {
        Info<< "Initializing standard phase change model" << endl;
        const word modelType = "standardPhaseChange";
        phaseChangePtr = phaseChangeModel::New(modelType, mesh, T, alpha1, dict);
    }

    Info<< "Initializing laser model" << endl;
    autoPtr<femtosecondLaserModel> laserPtr(new femtosecondLaserModel(mesh, dict));
    
    Info<< "Initializing droplet model" << endl;
    autoPtr<dropletModel> dropletPtr(new dropletModel(mesh, dict, rho, U, isFemtosecond));

    // Model validation
    if (!phaseChangePtr->valid())
    {
        FatalError << "Phase change model is invalid" << exit(FatalError);
    }
    
    if (!laserPtr->valid())
    {
        FatalError << "Laser model is invalid" << exit(FatalError);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        
        // Custom implementation of setDeltaT.H to avoid redeclaration of CoNum
        if (adjustTimeStep)
        {
            // Get the interface Courant number from CourantNo.H
            scalar localAlphaCoNum = 0.0;
            
            if (mesh.nInternalFaces())
            {
                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phi))().primitiveField()
                );

                localAlphaCoNum = 0.5*gMax
                (
                    sumPhi/(mesh.V().field()*mesh.surfaceInterpolation::deltaCoeffs().primitiveField())
                )*runTime.deltaTValue();
            }
            
            // Calculate time step based on the current Courant number
            scalar maxDeltaTFact = min(maxCo/(CoNum + SMALL), maxAlphaCo/(localAlphaCoNum + SMALL));
            
            // Limit time step change
            maxDeltaTFact = min(maxDeltaTFact, 1.2);
            
            // Set the new time step
            runTime.setDeltaT
            (
                min
                (
                    maxDeltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            );
            
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Mesh motion flag
        bool moveMesh = mesh.update();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMesh)
            {
                // Update mesh and fluxes
                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
                }

                // Update mixture properties
                mixture.correct();
                
                // Update density for femtosecond mode
                if (isFemtosecond)
                {
                    rho = alpha1*mixture.rho1() + (1.0 - alpha1)*mixture.rho2();
                    solve(fvm::ddt(rho) + fvc::div(phi*fvc::interpolate(rho)));
                }
            }

            // Phase fraction equation with subcycling
            const dictionary& alphaControls = mesh.solutionDict(alpha1.name());
            const scalar cAlpha = alphaControls.getOrDefault<scalar>("cAlpha", 1.0);
            const label nAlphaSubCycles = alphaControls.getOrDefault<label>("nAlphaSubCycles", 1);
            
            // Apply subcycling to improve stability for alpha equation
            if (nAlphaSubCycles > 1)
            {
                dimensionedScalar totalDeltaT = runTime.deltaT();
                surfaceScalarField rhoPhiSum
                (
                    IOobject
                    (
                        "rhoPhiSum",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", rhoPhi.dimensions(), 0)
                );

                // Loop through sub-cycles
                for
                (
                    subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                    !(++alphaSubCycle).end();
                )
                {
                    // Solve alpha equation for this subcycle
                    #include "alphaEqn.H"
                    
                    // Accumulate weighted rhoPhi for the timestep
                    rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
                }

                // Correct the flux using the accumulated values
                rhoPhi = rhoPhiSum;
            }
            else
            {
                // No subcycling, just solve the alpha equation once
                #include "alphaEqn.H"
            }

            // Update mixture
            mixture.correct();

            // Momentum equation
            #include "UEqn.H"

            // Temperature equation including laser heating
            laserPtr->correct();
            {
                // Get thermal diffusivity - use a default value
                scalar thermalDiffValue = 1e-5;  // Default value in m²/s
                
                // Try to read from the same dictionary we passed to the model initializations
                if (dict.found("thermalDiffusivity"))
                {
                    thermalDiffValue = dict.get<scalar>("thermalDiffusivity");
                }
                
                dimensionedScalar thermalDiff("alpha", dimViscosity, thermalDiffValue);
                
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T) - fvm::laplacian(thermalDiff, T)
                  ==
                    laserPtr->source()
                );

                if (isFemtosecond)
                {
                    // Add two-temperature coupling for femtosecond mode
                    TEqn += phaseChangePtr->electronSource();
                }

                TEqn.relax();
                TEqn.solve();
            }

            // Update phase change and droplets
            phaseChangePtr->correct(T);
            dropletPtr->update(T, p, U);

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


