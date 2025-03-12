// createPhi.H
Info<< "Creating face flux field phi directly\n" << endl;

// Create phi field with correct dimensions and register it to the registry
surfaceScalarField* phiPtr = new surfaceScalarField
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume/dimTime, 0.0)
);

// Explicitly store it in the registry
mesh.objectRegistry::store(phiPtr);

// Get a reference to the stored phi field
surfaceScalarField& phi = *phiPtr;

// Initialize with velocity field
phi = fvc::flux(U);

// Report the result
if (mesh.foundObject<surfaceScalarField>("phi"))
{
    Info<< "Successfully created and registered phi field\n" << endl;
}
else
{
    Info<< "Warning: Failed to register phi field\n" << endl;
}
