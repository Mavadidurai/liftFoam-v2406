#ifndef vacuumBoundaryCondition_H
#define vacuumBoundaryCondition_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class vacuumBoundaryCondition
:
    public fixedValueFvPatchScalarField
{
    // Private data members
    scalar p0_;              // Reference pressure [Pa]
    scalar Lv_;              // Latent heat of vaporization [J/kg]
    scalar R_;               // Gas constant [J/kg/K]
    scalar Tv_;              // Vapor temperature [K]
    scalar vacuumPressure_;  // Vacuum pressure [Pa]
    scalar minDensity_;      // Minimum allowed density [kg/m³]

public:
    //- Runtime type information
    TypeName("vacuum");

    // Constructors

        //- Construct from patch and internal field
        vacuumBoundaryCondition
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        vacuumBoundaryCondition
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given vacuumBoundaryCondition onto a new patch
        vacuumBoundaryCondition
        (
            const vacuumBoundaryCondition&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        vacuumBoundaryCondition
        (
            const vacuumBoundaryCondition&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new vacuumBoundaryCondition(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new vacuumBoundaryCondition(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

}

#endif
