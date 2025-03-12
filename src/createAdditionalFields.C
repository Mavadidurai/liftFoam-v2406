#include "createAdditionalFields.H"

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
        Info << "Created cellLevel field" << endl;
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
        Info << "Created pointLevel field" << endl;
    }
}

bool fieldExists(const fvMesh& mesh, const word& fieldName)
{
    return mesh.foundObject<regIOobject>(fieldName);
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

    const wordList requiredFields = 
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
            dimensionSet dims(0, 0, 0, 0, 0, 0, 0);
            
            // Set appropriate dimensions
            if (fieldName == "U")
                dims = dimensionSet(0, 1, -1, 0, 0, 0, 0);  // m/s
            else if (fieldName == "p")
                dims = dimensionSet(1, -1, -2, 0, 0, 0, 0); // kg/(m*s²)
            else if (fieldName == "alpha.titanium")
                dims = dimensionSet(0, 0, 0, 0, 0, 0, 0);   // dimensionless
            else if (fieldName == "T")
                dims = dimensionSet(0, 0, 0, 1, 0, 0, 0);   // K
            else if (fieldName == "rho")
                dims = dimensionSet(1, -3, 0, 0, 0, 0, 0);  // kg/m³
            
            auto fieldPtr = new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dims, 0.0)
            );
            mesh.objectRegistry::store(fieldPtr);
            Info << "Created field: " << fieldName << endl;
        }
    }
}

} // End namespace Foam
