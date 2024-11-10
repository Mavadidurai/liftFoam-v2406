// nonEquilibriumBoundaryCondition.H
#ifndef nonEquilibriumBoundaryCondition_H
#define nonEquilibriumBoundaryCondition_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "dimensionSet.H"

namespace Foam
{

class nonEquilibriumBoundaryCondition
:
    public fixedValueFvPatchScalarField
{
private:
    // Private data
    scalar Te_;    // Electron temperature
    scalar Tl_;    // Lattice temperature

public:
    //- Runtime type information
    TypeName("nonEquilibrium");

    // Constructors

        //- Construct from patch and internal field
        nonEquilibriumBoundaryCondition
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonEquilibriumBoundaryCondition
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given boundary condition onto new patch
        nonEquilibriumBoundaryCondition
        (
            const nonEquilibriumBoundaryCondition&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        nonEquilibriumBoundaryCondition
        (
            const nonEquilibriumBoundaryCondition&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEquilibriumBoundaryCondition(*this)
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
                new nonEquilibriumBoundaryCondition(*this, iF)
            );
        }

        //- Destructor
        virtual ~nonEquilibriumBoundaryCondition() = default;

    // Member Functions

        //- Return electron temperature
        scalar Te() const
        {
            return Te_;
        }

        //- Return lattice temperature
        scalar Tl() const
        {
            return Tl_;
        }

        //- Check if boundary condition is valid
        bool valid() const
        {
            // Check if temperatures are physically meaningful
            bool tempValid = Te_ > 0 && Tl_ > 0;

            // Check if patch values are within reasonable bounds
            bool valueValid = true;
            const scalarField& values = *this;
            forAll(values, i)
            {
                if (values[i] < 0)
                {
                    valueValid = false;
                    break;
                }
            }

            // Check if dimensions are consistent
            bool dimValid = true;
            const dimensionSet dimTemp(0, 0, 0, 1, 0, 0, 0); // Temperature dimensions
            if (this->internalField().dimensions() != dimTemp)
            {
                dimValid = false;
            }

            if (!tempValid)
            {
                Info<< "NonEquilibrium BC: Invalid temperatures" << endl;
            }
            if (!valueValid)
            {
                Info<< "NonEquilibrium BC: Invalid patch values" << endl;
            }
            if (!dimValid)
            {
                Info<< "NonEquilibrium BC: Invalid dimensions" << endl;
            }

            return tempValid && valueValid && dimValid;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

}

#endif
