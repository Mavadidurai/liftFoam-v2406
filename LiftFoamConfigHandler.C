#include "LiftFoamConfigHandler.H"

namespace Foam
{

LiftFoamConfigHandler::LiftFoamConfigHandler
(
    const Time& runTime,
    const fileName& dictPath
)
:
    runTime_(runTime),
    dict_
    (
        IOobject
        (
            dictPath,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    Info<< "Reading configuration from " << dictPath << endl;
}

bool LiftFoamConfigHandler::validateConfiguration() const
{
    bool isValid = true;

    Info<< "Validating LIFT configuration..." << endl;

    // Check numerical parameters
    isValid &= validateNumericalParameters();
    
    // Check laser parameters
    isValid &= validateLaserParameters();
    
    // Check material parameters
    isValid &= validateMaterialParameters();

    return isValid;
}

bool LiftFoamConfigHandler::validateFields() const
{
    const fvMesh& mesh = runTime_.lookupObject<fvMesh>(polyMesh::defaultRegion);
    bool isValid = true;

    Info<< "Validating field requirements..." << endl;

    // Required fields and their dimensions
    std::map<word, dimensionSet> requiredFields = {
        {"Te", dimTemperature},
        {"Tl", dimTemperature},
        {"alpha.titanium", dimless},
        {"laserSource", dimPower/dimVolume},
        {"phaseChangeRate", dimMass/dimVolume/dimTime},
        {"p", dimPressure},
        {"U", dimVelocity}
    };

    // Check field existence and dimensions
    for (const auto& field : requiredFields)
    {
        if (!mesh.foundObject<volScalarField>(field.first))
        {
            Info<< "Warning: Required field " << field.first << " not found" << endl;
            isValid = false;
            continue;
        }

        const volScalarField& vField = 
            mesh.lookupObject<volScalarField>(field.first);
        
        if (vField.dimensions() != field.second)
        {
            Info<< "Error: Field " << field.first 
                << " has incorrect dimensions" << nl
                << "Expected: " << field.second << nl
                << "Found: " << vField.dimensions() << endl;
            isValid = false;
        }
    }

    return isValid && validateFieldDimensions(mesh);
}

bool LiftFoamConfigHandler::validateLaserParameters() const
{
    bool isValid = true;

    // Check laser power
    scalar laserPower = getEntry<scalar>("laserPower", 2.0);
    if (laserPower <= 0)
    {
        Info<< "Error: Invalid laser power: " << laserPower << endl;
        isValid = false;
    }

    // Check beam parameters
    scalar beamRadius = getEntry<scalar>("beamRadius", 1.88e-3);
    if (beamRadius <= 0)
    {
        Info<< "Error: Invalid beam radius: " << beamRadius << endl;
        isValid = false;
    }

    // Check pulse parameters
    scalar pulseEnergy = getEntry<scalar>("pulseEnergy", 5e-6);
    scalar pulseDuration = getEntry<scalar>("pulseDuration", 214e-15);
    scalar pulseFrequency = getEntry<scalar>("pulseFrequency", 400e3);

    if (pulseEnergy <= 0 || pulseDuration <= 0 || pulseFrequency <= 0)
    {
        Info<< "Error: Invalid pulse parameters" << endl;
        isValid = false;
    }

    // Check absorption coefficient
    scalar absorptionCoeff = getEntry<scalar>("absorptionCoeff", 1e8);
    if (absorptionCoeff <= 0)
    {
        Info<< "Error: Invalid absorption coefficient" << endl;
        isValid = false;
    }

    // Check laser position
    vector laserPosition = getEntry<vector>("laserPosition", vector(0, 0, 0.005));
    if (!std::isfinite(mag(laserPosition)))
    {
        Info<< "Error: Invalid laser position" << endl;
        isValid = false;
    }

    return isValid;
}

bool LiftFoamConfigHandler::validateMaterialParameters() const
{
    bool isValid = true;

    // Check melting temperature
    scalar meltTemp = getEntry<scalar>("meltingTemperature", 1941);
    if (meltTemp <= 0)
    {
        Info<< "Error: Invalid melting temperature" << endl;
        isValid = false;
    }

    // Check latent heat
    scalar latentHeat = getEntry<scalar>("latentHeat", 3.65e5);
    if (latentHeat <= 0)
    {
        Info<< "Error: Invalid latent heat" << endl;
        isValid = false;
    }

    // Check undercooling coefficient
    scalar undercoolingCoeff = getEntry<scalar>("undercoolingCoeff", 1e-5);
    if (undercoolingCoeff < 0)
    {
        Info<< "Error: Invalid undercooling coefficient" << endl;
        isValid = false;
    }

    return isValid;
}

bool LiftFoamConfigHandler::validateNumericalParameters() const
{
    bool isValid = true;

    // Check phase fraction parameters
    isValid &= checkEntry<scalar>("nAlphaCorr", 1);
    isValid &= checkEntry<scalar>("nAlphaSubCycles", 2);
    isValid &= checkEntry<scalar>("cAlpha", 1.0);

    // Check time step controls
    scalar maxDeltaT = getEntry<scalar>("maxDeltaT", GREAT);
    scalar minDeltaT = getEntry<scalar>("minDeltaT", SMALL);

    if (minDeltaT <= 0 || maxDeltaT <= minDeltaT)
    {
        Info<< "Error: Invalid time step limits" << endl;
        isValid = false;
    }

    // Check electron-phonon coupling parameters
    isValid &= checkEntry<scalar>("Ce", 1000);    // Electron heat capacity
    isValid &= checkEntry<scalar>("Cl", 1000);    // Lattice heat capacity
    isValid &= checkEntry<scalar>("gamma", 1e17);  // Coupling factor

    return isValid;
}

bool LiftFoamConfigHandler::validateFieldDimensions(const fvMesh& mesh) const
{
    bool isValid = true;

    if (mesh.foundObject<volScalarField>("Te"))
    {
        const volScalarField& Te = mesh.lookupObject<volScalarField>("Te");
        if (min(Te).value() < 0)
        {
            Info<< "Error: Negative electron temperature detected" << endl;
            isValid = false;
        }
    }

    if (mesh.foundObject<volScalarField>("alpha.titanium"))
    {
        const volScalarField& alpha = 
            mesh.lookupObject<volScalarField>("alpha.titanium");
        if (min(alpha).value() < 0 || max(alpha).value() > 1)
        {
            Info<< "Error: Phase fraction out of bounds" << endl;
            isValid = false;
        }
    }

    return isValid;
}

}
