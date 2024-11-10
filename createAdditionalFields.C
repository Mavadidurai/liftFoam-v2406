#include "fvCFD.H"
#include "createAdditionalFields.H"

void createAdditionalFields(fvMesh& mesh)
{
    Info << "Creating additional fields for LIFT process\n" << endl;

    // Lambda function to create a scalar field if it doesn't already exist
    auto createFieldIfNeeded = [&](const word& fieldName,
                                   const dimensionSet& dims,
                                   const scalar& initialValue)
    {
        if (!mesh.foundObject<volScalarField>(fieldName))
        {
            Info << "Creating field " << fieldName << endl;
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(fieldName, dims, initialValue)
                )
            );
        }
    };

    // Create necessary fields if they don't already exist in the mesh
    createFieldIfNeeded("Te", dimTemperature, 300.0);                       // Electron temperature
    createFieldIfNeeded("Tl", dimTemperature, 300.0);                       // Lattice temperature
    createFieldIfNeeded("laserSource", dimPower/dimVolume, 0.0);            // Laser source term
    createFieldIfNeeded("phaseChangeRate", dimMass/dimVolume/dimTime, 0.0); // Phase change rate
    createFieldIfNeeded("phaseIndicator", dimless, 0.0);                    // Phase indicator
    createFieldIfNeeded("alpha1", dimless, 0.0);                            // Phase fraction for titanium
    createFieldIfNeeded("rho", dimDensity, 1.225);                          // Density field
    createFieldIfNeeded("p_rgh", dimPressure, 0.0);                         // Pressure field p_rgh
    createFieldIfNeeded("p", dimPressure, 1e5);                             // Absolute pressure field p
}

