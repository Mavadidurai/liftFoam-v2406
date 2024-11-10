#include "DimensionValidator.H"

namespace Foam
{

// Initialize static dimension sets
const dimensionSet DimensionValidator::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimAcceleration(0, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimPressure(1, -1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet DimensionValidator::dimEnergy(1, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimPower(1, 2, -3, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimEnergyDensity(1, -1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimFrequency(0, 0, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimSurfaceTension(1, 0, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimThermalCond(1, 1, -3, -1, 0, 0, 0);
const dimensionSet DimensionValidator::dimHeatCapacity(0, 2, -2, -1, 0, 0, 0);
const dimensionSet DimensionValidator::dimViscosity(0, 2, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimLaserSource(1, -1, -3, 0, 0, 0, 0);

void DimensionValidator::checkFieldDimensions
(
    const dimensionSet& actualDimensions,
    const dimensionSet& expectedDimensions,
    const word& fieldName
)
{
    if (actualDimensions != expectedDimensions)
    {
        FatalErrorInFunction
            << "Dimension mismatch for field " << fieldName << nl
            << "Expected: " << expectedDimensions << nl
            << "Actual: " << actualDimensions << nl
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
    // Basic field validity checks
    if (!Te.size() || !Tl.size() || !p.size() || !U.size() 
        || !alpha1.size() || !laserSource.size() || !phaseChangeRate.size())
    {
        FatalErrorInFunction
            << "One or more fields have invalid size"
            << abort(FatalError);
    }

    // Check field dimensions
    checkFieldDimensions(Te, dimTemperature, "electron temperature");
    checkFieldDimensions(Tl, dimTemperature, "lattice temperature");
    checkFieldDimensions(p, dimPressure, "pressure");
    checkFieldDimensions(U, dimVelocity, "velocity");
    checkFieldDimensions(alpha1, dimless, "phase fraction");
    checkFieldDimensions(laserSource, dimLaserSource, "laser source");
    checkFieldDimensions
    (
        phaseChangeRate,
        dimMass/dimVolume/dimTime,
        "phase change rate"
    );

    // Check for physical bounds
    if (min(Te).value() < 0 || min(Tl).value() < 0)
    {
        FatalErrorInFunction
            << "Negative temperature detected" << nl
            << "Min Te: " << min(Te).value() << nl
            << "Min Tl: " << min(Tl).value()
            << abort(FatalError);
    }

    if (min(p).value() < 0)
    {
        FatalErrorInFunction
            << "Negative pressure detected" << nl
            << "Min p: " << min(p).value()
            << abort(FatalError);
    }

    // Check phase fraction bounds
    if (min(alpha1).value() < -SMALL || max(alpha1).value() > 1 + SMALL)
    {
        FatalErrorInFunction
            << "Phase fraction out of bounds [0,1]" << nl
            << "Min alpha: " << min(alpha1).value() << nl
            << "Max alpha: " << max(alpha1).value()
            << abort(FatalError);
    }

    // Check for NaN/Inf values using internal field
    scalar maxTeVal = max(mag(Te.internalField())).value();
    scalar maxTlVal = max(mag(Tl.internalField())).value();
    scalar maxPVal = max(mag(p.internalField())).value();
    scalar maxUVal = max(mag(U.internalField())).value();

    if (maxTeVal > GREAT || maxTlVal > GREAT)
    {
        FatalErrorInFunction
            << "Temperature field contains invalid values" << nl
            << "Max Te magnitude: " << maxTeVal << nl
            << "Max Tl magnitude: " << maxTlVal
            << abort(FatalError);
    }

    if (maxPVal > GREAT)
    {
        FatalErrorInFunction
            << "Pressure field contains invalid values" << nl
            << "Max pressure magnitude: " << maxPVal
            << abort(FatalError);
    }

    if (maxUVal > GREAT)
    {
        FatalErrorInFunction
            << "Velocity field contains invalid values" << nl
            << "Max velocity magnitude: " << maxUVal
            << abort(FatalError);
    }

    Info<< "All LIFT field dimensions and bounds validated" << endl;
}
void DimensionValidator::checkTTMDimensions
(
    const dimensionedScalar& Ce,
    const dimensionedScalar& Cl,
    const dimensionedScalar& G
)
{
    // Check for valid dimensioned scalars
    if (Ce.value() <= 0 || Cl.value() <= 0 || G.value() < 0)
    {
        FatalErrorInFunction
            << "Invalid TTM parameters:" << nl
            << "Ce = " << Ce.value() << nl
            << "Cl = " << Cl.value() << nl
            << "G = " << G.value() << nl
            << "All parameters must be positive"
            << abort(FatalError);
    }

    checkFieldDimensions
    (
        Ce.dimensions(),
        dimEnergy/dimVolume/dimTemperature,
        "electron heat capacity"
    );
    
    checkFieldDimensions
    (
        Cl.dimensions(),
        dimEnergy/dimVolume/dimTemperature,
        "lattice heat capacity"
    );
    
    checkFieldDimensions
    (
        G.dimensions(),
        dimEnergy/dimVolume/dimTime/dimTemperature,
        "electron-phonon coupling"
    );
}

void DimensionValidator::checkEnergyConservation
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volScalarField& laserSource,
    const volScalarField& phaseChangeRate,
    const dimensionedScalar& timeStep
)
{
    // Check for valid time step
    if (timeStep.value() <= 0)
    {
        FatalErrorInFunction
            << "Invalid time step value: " << timeStep.value()
            << abort(FatalError);
    }

    // Check field dimensions
    checkFieldDimensions(Te, dimTemperature, "electron temperature");
    checkFieldDimensions(Tl, dimTemperature, "lattice temperature");
    checkFieldDimensions(laserSource, dimPower/dimVolume, "laser source");
    checkFieldDimensions(phaseChangeRate, dimMass/dimVolume/dimTime, "phase change rate");

    dimensionedScalar totalEnergy = fvc::domainIntegrate
    (
        Te + Tl + timeStep*laserSource - timeStep*phaseChangeRate
    );

    if (totalEnergy.value() < 0)
    {
        FatalErrorInFunction
            << "Energy conservation violation detected" << nl
            << "Total energy: " << totalEnergy.value()
            << abort(FatalError);
    }
}

void DimensionValidator::checkInterfaceResolution
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
            << ". Must be >= 1"
            << abort(FatalError);
    }

    scalar interfaceThickness = 0.0;
    label interfaceCells = 0;

    forAll(mesh.C(), cellI)
    {
        if (alpha1[cellI] > 0.01 && alpha1[cellI] < 0.99)
        {
            interfaceCells++;
            interfaceThickness += pow(mesh.V()[cellI], 1.0/3.0);
        }
    }

    if (interfaceCells > 0)
    {
        interfaceThickness /= interfaceCells;
        
        if (interfaceCells < minCellsInInterface)
        {
            FatalErrorInFunction
                << "Insufficient interface resolution" << nl
                << "Cells in interface: " << interfaceCells << nl
                << "Minimum required: " << minCellsInInterface << nl
                << "Average interface thickness: " << interfaceThickness
                << abort(FatalError);
        }
    }

    Info<< "Interface resolution check passed with " << interfaceCells 
        << " interface cells" << endl;
}

void DimensionValidator::checkBoundaryConditions
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volScalarField& p,
    const dictionary& dict
)
{
    if (!dict.found("vacuumPressure") || !dict.found("minTemperature"))
    {
        Warning
            << "Dictionary missing vacuumPressure or minTemperature entries. "
            << "Using defaults: vacuumPressure=1e-6, minTemperature=300.0"
            << endl;
    }

    const scalar vacuumPressure = dict.lookupOrDefault<scalar>("vacuumPressure", 1e-6);
    const scalar minTemp = dict.lookupOrDefault<scalar>("minTemperature", 300.0);

    // Check boundary fields
    forAll(Te.boundaryField(), patchI)
    {
        const fvPatchScalarField& Tep = Te.boundaryField()[patchI];
        const fvPatchScalarField& Tlp = Tl.boundaryField()[patchI];
        const fvPatchScalarField& pp = p.boundaryField()[patchI];

        if (min(Tep) < minTemp || min(Tlp) < minTemp)
        {
            FatalErrorInFunction
                << "Temperature below minimum at boundary " 
                << Tep.patch().name() << nl
                << "Min Te: " << min(Tep) << nl
                << "Min Tl: " << min(Tlp) << nl
                << "Minimum allowed: " << minTemp
                << abort(FatalError);
        }

        if (Tep.patch().name() == "vacuum" && min(pp) < vacuumPressure)
        {
            FatalErrorInFunction
                << "Pressure below vacuum pressure at vacuum boundary" << nl
                << "Min pressure: " << min(pp) << nl
                << "Vacuum pressure: " << vacuumPressure
                << abort(FatalError);
        }
    }

    Info<< "Boundary conditions validated successfully" << endl;
}

void DimensionValidator::checkPhaseChangeDimensions
(
    const dimensionedScalar& latentHeat,
    const dimensionedScalar& phaseChangeRate,
    const dimensionedScalar& interfaceEnergy
)
{
    // Check for valid values
    if (latentHeat.value() <= 0)
    {
        FatalErrorInFunction
            << "Invalid latent heat: " << latentHeat.value()
            << ". Must be positive"
            << abort(FatalError);
    }

    checkFieldDimensions
    (
        latentHeat.dimensions(),
        dimEnergy/dimMass,
        "latent heat"
    );

    checkFieldDimensions
    (
        phaseChangeRate.dimensions(),
        dimMass/dimVolume/dimTime,
        "phase change rate"
    );

    checkFieldDimensions
    (
        interfaceEnergy.dimensions(),
        dimEnergy/dimArea,
        "interface energy"
    );

    Info<< "Phase change dimensions validated" << endl;
}

} // End namespace Foam
