#include "createAdditionalFields.H"
#include "DimensionValidator.H"

namespace Foam
{

void createMandatoryFields(fvMesh& mesh)
{
    // Essential fields that must exist before mesh operations
    if (!fieldExists(mesh, "cellLevel"))
    {
        auto cellLevelPtr = new volScalarField
        (
            IOobject
            (
                "cellLevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );
        mesh.objectRegistry::store(cellLevelPtr);
    }

    if (!fieldExists(mesh, "pointLevel")) 
    {
        auto pointLevelPtr = new volScalarField
        (
            IOobject
            (
                "pointLevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );
        mesh.objectRegistry::store(pointLevelPtr);
    }
}

bool fieldExists(const fvMesh& mesh, const word& fieldName)
{
    return mesh.foundObject<volScalarField>(fieldName);
}

bool validateField
(
    const volScalarField& field,
    const word& fieldName, 
    const dimensionSet& expectedDims
)
{
    // Change match() to operator== for dimension comparison
    if (field.dimensions() != expectedDims)  
    {
        return false;
    }

    forAll(field, cellI)
    {
        if (!std::isfinite(field[cellI]))
        {
            return false;
        }
    }
    return true;
}

void createAdditionalFields(fvMesh& mesh)
{
    // Only proceed if mandatory fields exist
    if (!fieldExists(mesh, "cellLevel") || !fieldExists(mesh, "pointLevel"))
    {
        FatalError
            << "Mandatory mesh fields missing. Run createMandatoryFields first."
            << abort(FatalError);
    }

    const word requiredFields[] = 
    {
        "U",
        "p",
        "alpha.titanium",
        "T",
        "rho"
    };

    // Create essential fields first
    for (const word& fieldName : requiredFields)
    {
        if (!fieldExists(mesh, fieldName))
        {
            auto fieldPtr = new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,  // Changed from MUST_READ_IF_PRESENT
                    IOobject::AUTO_WRITE
                ),
                mesh
            );
            mesh.objectRegistry::store(fieldPtr);
        }
    }
}

} // End namespace Foam
