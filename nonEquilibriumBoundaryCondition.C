// nonEquilibriumBoundaryCondition.C
#include "nonEquilibriumBoundaryCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

defineTypeNameAndDebug(nonEquilibriumBoundaryCondition, 0);

// Add constructors to run-time selection table
addToRunTimeSelectionTable
(
    fvPatchScalarField,
    nonEquilibriumBoundaryCondition,
    dictionary
);

nonEquilibriumBoundaryCondition::nonEquilibriumBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Te_(300),    // Default electron temperature
    Tl_(300)     // Default lattice temperature
{
    if (debug)
    {
        Info<< "Creating " << type() << " boundary condition" << endl;
    }
}

nonEquilibriumBoundaryCondition::nonEquilibriumBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Te_(dict.getOrDefault<scalar>("Te", 300)),
    Tl_(dict.getOrDefault<scalar>("Tl", 300))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(Te_);
    }

    // Validate temperatures
    if (Te_ <= 0 || Tl_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid temperatures specified: Te = " << Te_
            << ", Tl = " << Tl_ << nl
            << "Both temperatures must be positive"
            << exit(FatalIOError);
    }
}

nonEquilibriumBoundaryCondition::nonEquilibriumBoundaryCondition
(
    const nonEquilibriumBoundaryCondition& bc,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(bc, p, iF, mapper),
    Te_(bc.Te_),
    Tl_(bc.Tl_)
{}

nonEquilibriumBoundaryCondition::nonEquilibriumBoundaryCondition
(
    const nonEquilibriumBoundaryCondition& bc,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(bc, iF),
    Te_(bc.Te_),
    Tl_(bc.Tl_)
{}

void nonEquilibriumBoundaryCondition::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar relaxationFactor = 0.5;
    const scalar Teq = (Te_ + Tl_) / 2.0;

    scalarField& field = *this;
    forAll(field, i)
    {
        field[i] = Teq + relaxationFactor*(field[i] - Teq);
        if (field[i] < 0)
        {
            field[i] = 0;
            Info<< "Warning: Corrected negative temperature at patch face " << i << endl;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void nonEquilibriumBoundaryCondition::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    
    os.writeEntry("type", type());
    os.writeEntry("Te", Te_);
    os.writeEntry("Tl", Tl_);
    writeEntry("value", os);
}

}
