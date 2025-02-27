/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Class for validating field dimensions and physical parameters in LIFT process.
    
    Handles:
    - Field dimension validation
    - Physical parameter bounds checking
    - Energy conservation validation
    - Interface resolution checking
    
    Features:
    - Standard dimension sets for common quantities
    - Template functions for field validation
    - Interface checks for mesh quality
\*---------------------------------------------------------------------------*/

#ifndef DimensionValidator_H
#define DimensionValidator_H

#include "fvCFD.H"
#include "volFields.H"
#include "dimensionSet.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "volMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DimensionValidator Declaration
\*---------------------------------------------------------------------------*/

class DimensionValidator
{
public:
    //- Common dimension sets
    static const dimensionSet dimVelocity;        // [m/s]
    static const dimensionSet dimPressure;        // [kg⋅m⁻¹⋅s⁻²]
    static const dimensionSet dimTemperature;     // [K]
    static const dimensionSet dimEnergy;          // [kg⋅m²⋅s⁻²]
    static const dimensionSet dimLaserSource;     // [kg⋅m²⋅s⁻³⋅m⁻³]

    //- Check field dimensions
    template<class Type>
    static void checkFieldDimensions
    (
        const GeometricField<Type, fvPatchField, volMesh>& field,
        const dimensionSet& expectedDimensions,
        const word& fieldName
    );

    //- Check LIFT fields dimensions and bounds
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

    //- Validate LIFT fields
    static void validateLIFTFields
    (
        const fvMesh& mesh,
        const volScalarField& Te,
        const volScalarField& Tl,
        const volScalarField& p,
        const volVectorField& U,
        const volScalarField& alpha1,
        const volScalarField& laserSource
    );

    //- Check dimension set
    static void checkDimensionSet
    (
        const dimensionSet& actual,
        const dimensionSet& expected,
        const word& paramName
    );

    //- Validate two-temperature model parameters
    static void validateTTMParameters
    (
        const dimensionedScalar& Ce,
        const dimensionedScalar& Cl,
        const dimensionedScalar& G
    );

    //- Validate energy conservation
    static void validateEnergyConservation
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volScalarField& laserSource,
        const volScalarField& phaseChangeRate,
        const dimensionedScalar& timeStep
    );

    //- Validate interface resolution
    static void validateInterfaceResolution
    (
        const volScalarField& alpha1,
        const fvMesh& mesh,
        const scalar minCellsInInterface
    );
};

} // End namespace Foam

#endif
