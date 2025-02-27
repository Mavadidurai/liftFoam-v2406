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

// LIFT-specific includes
#include "twoTemperatureModel.H"
#include "phaseChangeModel.H"
#include "femtosecondLaserModel.H"
#include "dropletModel.H"

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
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    
    // Read regime from dictionary
    const word regime(mesh.solutionDict().lookupOrDefault<word>("regime", "regular"));
    const bool isFemtosecond = (regime == "femtosecond");

    Info<< "Running in " << regime << " regime\n" << endl;

    // Initialize base fields
    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    
    // Set pressure reference
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p_rgh.name());

    // Initialize density if in femtosecond mode
    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    );

    if (isFemtosecond)
    {
        rho = alpha1*mixture.rho1() + (1.0 - alpha1)*mixture.rho2();
    }
            
    #include "createTimeControls.H"
    #include "initializeFields.H"
    #include "validateFields.H"

    // Initialize models
    autoPtr<phaseChangeModel> phaseChange(phaseChangeModel::New(mesh, T, alpha1, dict));
    autoPtr<femtosecondLaserModel> laser(femtosecondLaserModel::New(mesh, dict));
    autoPtr<dropletModel> droplet(dropletModel::New(mesh, dict, rho, U, isFemtosecond));

    // Model validation
    turbulence->validate();
    phaseChange->validate();
    laser->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

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
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            // Update mixture
            mixture.correct();

            // Momentum equation
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(rhoPhi, U)
              - fvm::laplacian(mixture.mu(), U)
              - (fvc::grad(U) & fvc::grad(mixture.mu()))
            );

            UEqn.relax();

            if (pimple.momentumPredictor())
            {
                solve
                (
                    UEqn
                  ==
                    fvc::reconstruct
                    (
                        (
                            mixture.surfaceTensionForce()
                          - ghf*fvc::snGrad(rho)
                          - fvc::snGrad(p_rgh)
                        ) * mesh.magSf()
                    )
                );
            }

            // Temperature equation including laser heating
            laser->correct();
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T) - fvm::laplacian(mixture.alpha(), T)
                  ==
                    laser->source()/mixture.Cp()
                );

                if (isFemtosecond)
                {
                    // Add two-temperature coupling for femtosecond mode
                    TEqn += phaseChange->electronSource()/mixture.Cp();
                }

                TEqn.relax();
                TEqn.solve();
            }

            // Update phase change and droplets
            phaseChange->correct(T);
            droplet->update(T, p, U);

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                // Pressure equation
                volScalarField rAU(1.0/UEqn.A());
                surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
                
                volVectorField HbyA(U);
                HbyA = rAU*UEqn.H();

                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                  + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
                );

                if (p_rgh.needReference())
                {
                    fvc::makeRelative(phiHbyA, U);
                    adjustPhi(phiHbyA, U, p_rgh);
                    fvc::makeAbsolute(phiHbyA, U);
                }

                // Update pressure
                fvScalarMatrix p_rghEqn
                (
                    fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
                );

                p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
                p_rghEqn.solve();

                if (pimple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - p_rghEqn.flux();
                    U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
                    U.correctBoundaryConditions();
                }
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
