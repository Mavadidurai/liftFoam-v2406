#include "DimensionValidator.H"

namespace Foam
{

// Initialize static dimension sets
const dimensionSet DimensionValidator::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimPressure(1, -1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet DimensionValidator::dimEnergy(1, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimLaserSource(1, -1, -3, 0, 0, 0, 0);

template<class Type>
void DimensionValidator::checkFieldDimensions
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const dimensionSet& expectedDimensions,
    const word& fieldName
)
{
    if (field.dimensions() != expectedDimensions)
    {
        FatalErrorInFunction
            << "Dimension mismatch for field \"" << fieldName << "\"" << nl
            << "Expected: " << expectedDimensions << nl
            << "Actual: " << field.dimensions()
            << abort(FatalError);
    }
}
void DimensionValidator::checkLIFTFields
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& alpha1,
    const volScalarField& laserSource,
    const volScalarField& phaseChangeRate
)
{
    // Check dimensions
    checkFieldDimensions(Te, dimTemperature, "Te");
    checkFieldDimensions(Tl, dimTemperature, "Tl");
    checkFieldDimensions(p, dimPressure, "p");
    checkFieldDimensions(U, dimVelocity, "U");
    checkFieldDimensions(alpha1, dimless, "alpha1");
    checkFieldDimensions(laserSource, dimPower/dimVolume, "laserSource");
    checkFieldDimensions(phaseChangeRate, dimMass/dimVolume/dimTime, "phaseChangeRate");

    // Check bounds
    if (min(Te).value() < 0 || min(Tl).value() < 0)
    {
        FatalErrorIn("checkLIFTFields")
            << "Negative temperatures detected"
            << abort(FatalError);
    }

    if (min(alpha1).value() < -SMALL || max(alpha1).value() > 1 + SMALL)
    {
        FatalErrorIn("checkLIFTFields")
            << "Phase fraction out of bounds"
            << abort(FatalError);
    }
}
void DimensionValidator::validateLIFTFields
(
    const fvMesh& mesh,        // Add mesh parameter
    const volScalarField& Te,
    const volScalarField& Tl,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& alpha1,
    const volScalarField& laserSource
)
{
    // Validate dimensions
    checkFieldDimensions(Te, dimTemperature, "Electron Temperature");
    checkFieldDimensions(Tl, dimTemperature, "Lattice Temperature"); 
    checkFieldDimensions(p, dimPressure, "Pressure");
    checkFieldDimensions(U, dimVelocity, "Velocity");
    checkFieldDimensions(alpha1, dimless, "Phase Fraction");
    checkFieldDimensions(laserSource, dimLaserSource, "Laser Source");

    // Define and use the dimension sets
    const dimensionSet dimSurfaceTension(1, 0, -2, 0, 0, 0, 0);  // [N/m]
    const dimensionSet dimPhaseChange(1, -3, -1, 0, 0, 0, 0);    // [kg/m³/s]

    // Additional validation for surface tension and phase change fields if present
    if (mesh.foundObject<volScalarField>("surfaceTension"))
    {
        const volScalarField& sigma = 
            mesh.lookupObject<volScalarField>("surfaceTension");
        checkFieldDimensions(sigma, dimSurfaceTension, "Surface Tension");
    }

    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& pcr = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
        checkFieldDimensions(pcr, dimPhaseChange, "Phase Change Rate");
    }

    // Validate bounds
    if (min(Te).value() < 0 || min(Tl).value() < 0)
    {
        FatalErrorInFunction
            << "Negative temperature detected:" << nl
            << "Min Te: " << min(Te).value() << ", Min Tl: " << min(Tl).value()
            << abort(FatalError);
    }

    if (min(alpha1).value() < 0 || max(alpha1).value() > 1)
    {
        FatalErrorInFunction
            << "Phase fraction out of bounds [0, 1]" << nl
            << "Min alpha: " << min(alpha1).value() 
            << ", Max alpha: " << max(alpha1).value()
            << abort(FatalError);
    }

    Info << "LIFT field dimensions and bounds validated successfully." << endl;
}
void DimensionValidator::checkDimensionSet
(
    const dimensionSet& actual, 
    const dimensionSet& expected,
    const word& paramName
)  
{
    if (actual != expected)
    {
        FatalErrorInFunction
            << "Dimension mismatch for " << paramName << nl
            << "Expected: " << expected << nl
            << "Actual: " << actual 
            << abort(FatalError);
    }
}
void DimensionValidator::validateTTMParameters
(
    const dimensionedScalar& Ce,
    const dimensionedScalar& Cl,
    const dimensionedScalar& G
)
{
    if (Ce.value() <= 0 || Cl.value() <= 0 || G.value() <= 0)
    {
        FatalErrorInFunction
            << "Invalid TTM parameters: Ce=" << Ce.value() << ", Cl=" << Cl.value()
            << ", G=" << G.value() << ". All must be positive."
            << abort(FatalError);
    }

    // Changed from checkFieldDimensions to checkDimensionSet
    checkDimensionSet(Ce.dimensions(), 
                     dimEnergy/dimVolume/dimTemperature, 
                     "Electron Heat Capacity");

    checkDimensionSet(Cl.dimensions(),
                     dimEnergy/dimVolume/dimTemperature, 
                     "Lattice Heat Capacity");

    checkDimensionSet(G.dimensions(),
                     dimEnergy/dimVolume/dimTime/dimTemperature,
                     "Electron-Phonon Coupling");

    Info << "TTM parameters validated successfully." << endl;
}

void DimensionValidator::validateEnergyConservation
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volScalarField& laserSource,
    const volScalarField& phaseChangeRate,
    const dimensionedScalar& timeStep
)
{
    if (timeStep.value() <= 0)
    {
        FatalErrorInFunction
            << "Invalid time step: " << timeStep.value() << ". Must be positive."
            << abort(FatalError);
    }

    checkFieldDimensions(laserSource, dimLaserSource, "Laser Source");
    checkFieldDimensions(phaseChangeRate, dimMass/dimVolume/dimTime, "Phase Change Rate");

    Info << "Energy conservation validated successfully." << endl;
}

void DimensionValidator::validateInterfaceResolution
(
    const volScalarField& alpha1,
    const fvMesh& mesh,
    const scalar minCellsInInterface
)
{
    if (minCellsInInterface < 1)
    {
        FatalErrorInFunction
            << "Invalid minCellsInInterface value: " << minCellsInInterface
            << ". Must be >= 1."
            << abort(FatalError);
    }

    scalar interfaceCells = 0;
    forAll(mesh.C(), cellI)
    {
        if (alpha1[cellI] > 0.01 && alpha1[cellI] < 0.99)
        {
            interfaceCells++;
        }
    }

    if (interfaceCells < minCellsInInterface)
    {
        FatalErrorInFunction
            << "Insufficient interface resolution: " << interfaceCells
            << " cells (Minimum required: " << minCellsInInterface << ")"
            << abort(FatalError);
    }

    Info << "Interface resolution validated with " << interfaceCells << " cells." << endl;
}

} // End namespace Foam

