/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Application
    liftFoam

Description
    Solver for femtosecond Laser-Induced Forward Transfer (LIFT) process.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "fieldOperationsFix.H"
#include "LIFTModel.H"
#include "DimensionValidator.H"
#include "LiftFoamConfigHandler.H"

using namespace Foam;

int main(int argc, char *argv[])
{

    argList::addNote
    (
        "Solver for femtosecond Laser-Induced Forward Transfer (LIFT) process."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    // Create PIMPLE control
    pimpleControl pimple(mesh);

    // Create MRF zones
    IOMRFZoneList MRF(mesh);
    
    // Initialize control variables
    bool moveMesh = mesh.moving();
    bool correctPhi = pimple.dict().lookupOrDefault<bool>("correctPhi", false);
    bool checkMeshCourantNo = 
        pimple.dict().lookupOrDefault<bool>("checkMeshCourantNo", false);

    // Initialize time control variables
    bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
    scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.5);
    scalar maxAlphaCo = runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 0.2);
    scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    scalar minDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // Initialize LIFT model
    Info<< "\nInitializing LIFT model and configuration\n" << endl;
    
    LiftFoamConfigHandler config(runTime, "liftFoamDict");
    
    if (!config.validateConfiguration())
    {
        FatalError
            << "Invalid configuration in liftFoamDict"
            << abort(FatalError);
    }

    // Create LIFT model
    LIFTModel liftModel(mesh, config.dict());
    
    // Verify model initialization
    if (!liftModel.valid())
    {
        FatalError
            << "LIFT model initialization failed"
            << abort(FatalError);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update MRF velocities
        MRF.correctBoundaryVelocity(U);

        // PIMPLE loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMesh)
            {
                if (correctPhi)
                {
                    #include "correctPhi.H"
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }

            // Alpha sub-cycling
            {
                #include "alphaControls.H"
                #include "alphaEqnSubCycle.H"
            }

            // Momentum predictor
            tmp<fvVectorMatrix> tUEqn;
            {
                #include "UEqn.H"
            }

            // Pressure-velocity correction
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulencePtr->correct();
            }
        }

        // Update LIFT model
        try
        {
            liftModel.solve();
            liftModel.updateFields(alpha1, U, p, p_rgh, rho);
        }
        catch (const Foam::error& err)
        {
            FatalError
                << "Error in LIFT model solve: " << err.message()
                << abort(FatalError);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
