#include "vacuumBoundaryCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

defineTypeNameAndDebug(vacuumBoundaryCondition, 0);
addToRunTimeSelectionTable
(
    fvPatchScalarField,
    vacuumBoundaryCondition,
    dictionary
);

vacuumBoundaryCondition::vacuumBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(1e5),           // Default: 1 atm
    Lv_(2.26e6),        // Default: water vaporization
    R_(461.5),          // Default: water vapor gas constant
    Tv_(373.15),        // Default: water boiling point
    vacuumPressure_(1e-6),
    minDensity_(1e-6)
{}

vacuumBoundaryCondition::vacuumBoundaryCondition
(
    const vacuumBoundaryCondition& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_),
    Lv_(ptf.Lv_),
    R_(ptf.R_),
    Tv_(ptf.Tv_),
    vacuumPressure_(ptf.vacuumPressure_),
    minDensity_(ptf.minDensity_)
{}

vacuumBoundaryCondition::vacuumBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    p0_(dict.lookupOrDefault<scalar>("p0", 1e5)),
    Lv_(dict.lookupOrDefault<scalar>("Lv", 2.26e6)),
    R_(dict.lookupOrDefault<scalar>("R", 461.5)),
    Tv_(dict.lookupOrDefault<scalar>("Tv", 373.15)),
    vacuumPressure_(dict.lookupOrDefault<scalar>("vacuumPressure", 1e-6)),
    minDensity_(dict.lookupOrDefault<scalar>("minDensity", 1e-6))
{}

vacuumBoundaryCondition::vacuumBoundaryCondition
(
    const vacuumBoundaryCondition& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wbppsf, iF),
    p0_(wbppsf.p0_),
    Lv_(wbppsf.Lv_),
    R_(wbppsf.R_),
    Tv_(wbppsf.Tv_),
    vacuumPressure_(wbppsf.vacuumPressure_),
    minDensity_(wbppsf.minDensity_)
{}

void vacuumBoundaryCondition::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get patch fields
    scalarField& pf = *this;
    
    // Get temperature field
    const scalarField& Tf = 
        patch().lookupPatchField<volScalarField, scalar>("T");
        
    // Calculate vapor pressure using saturation equation
    scalarField pv(patch().size(), 0.0);
    forAll(patch(), faceI)
    {
        scalar T = max(Tf[faceI], scalar(1.0));  // Prevent division by zero
        pv[faceI] = p0_*exp(Lv_/R_*(1.0/Tv_ - 1.0/T));
    }
    
    // Get pressure field
    const fvPatchField<scalar>& pPatch = 
        patch().lookupPatchField<volScalarField, scalar>("p");
    
    // Update pressure with vapor pressure and vacuum limit
    forAll(patch(), faceI)
    {
        pf[faceI] = max(min(pv[faceI], pPatch[faceI]), vacuumPressure_);
    }
    
    // Update density if rho field exists
    if (db().foundObject<volScalarField>("rho"))
    {
        fvPatchField<scalar>& rhoPatch = 
            const_cast<fvPatchField<scalar>&>
            (
                patch().lookupPatchField<volScalarField, scalar>("rho")
            );
            
        forAll(patch(), faceI)
        {
            scalar T = max(Tf[faceI], scalar(1.0));  // Prevent division by zero
            rhoPatch[faceI] = max(pf[faceI]/(R_*T), minDensity_);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void vacuumBoundaryCondition::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    
    os.writeEntry("p0", p0_);
    os.writeEntry("Lv", Lv_);
    os.writeEntry("R", R_);
    os.writeEntry("Tv", Tv_);
    os.writeEntry("vacuumPressure", vacuumPressure_);
    os.writeEntry("minDensity", minDensity_);
    
    writeEntry("value", os);
}

} // End namespace Foam
