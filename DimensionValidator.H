#ifndef DimensionValidator_H
#define DimensionValidator_H

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionSet.H"
#include "fvPatchFields.H"

namespace Foam
{

class DimensionValidator
{
public:
    // Common dimension sets
    static const dimensionSet dimVelocity;        // [m/s]
    static const dimensionSet dimAcceleration;    // [m/s²]
    static const dimensionSet dimPressure;        // [kg/m/s²]
    static const dimensionSet dimTemperature;     // [K]
    static const dimensionSet dimEnergy;          // [kg⋅m²/s²]
    static const dimensionSet dimPower;           // [kg⋅m²/s³]
    static const dimensionSet dimEnergyDensity;   // [kg/m/s²]
    static const dimensionSet dimFrequency;       // [1/s]
    static const dimensionSet dimSurfaceTension;  // [kg/s²]
    static const dimensionSet dimThermalCond;     // [kg⋅m/s³/K]
    static const dimensionSet dimHeatCapacity;    // [m²/s²/K]
    static const dimensionSet dimViscosity;       // [m²/s]
    static const dimensionSet dimLaserSource;     // [kg⋅m²/s³/m³]

    //- Check dimensions of a dimensionSet against expected dimensions
    static void checkFieldDimensions
    (
        const dimensionSet& actualDimensions,
        const dimensionSet& expectedDimensions,
        const word& fieldName
    );

    //- Check dimensions of a GeometricField against expected dimensions
    template<class Type>
    static void checkFieldDimensions
    (
        const GeometricField<Type, fvPatchField, volMesh>& field,
        const dimensionSet& expectedDimensions,
        const word& fieldName
    )
    {
        checkFieldDimensions(field.dimensions(), expectedDimensions, fieldName);
    }

    //- Check LIFT solver specific field dimensions
    static void checkLIFTFields
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volScalarField& p,
        const volVectorField& U,
        const volScalarField& alpha1,
        const volScalarField& laserSource,
        const volScalarField& phaseChangeRate
    );

    //- Check two-temperature model dimensions
    static void checkTTMDimensions
    (
        const dimensionedScalar& Ce,
        const dimensionedScalar& Cl,
        const dimensionedScalar& G
    );

    //- Check phase change model dimensions  
    static void checkPhaseChangeDimensions
    (
        const dimensionedScalar& latentHeat,
        const dimensionedScalar& phaseChangeRate,
        const dimensionedScalar& interfaceEnergy
    );

    //- Check boundary conditions
    static void checkBoundaryConditions
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volScalarField& p,
        const dictionary& dict
    );

    //- Check interface resolution
    static void checkInterfaceResolution
    (
        const volScalarField& alpha1,
        const fvMesh& mesh,
        const scalar minCellsInInterface
    );

    //- Check energy conservation
    static void checkEnergyConservation
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volScalarField& laserSource,
        const volScalarField& phaseChangeRate,
        const dimensionedScalar& timeStep
    );
};

} // End namespace Foam

#endif
