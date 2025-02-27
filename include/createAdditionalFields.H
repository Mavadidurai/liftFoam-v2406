#ifndef createAdditionalFields_H
#define createAdditionalFields_H

#include "fvCFD.H"
#include "volFields.H"

namespace Foam
{
    // Check and create mandatory fields first
    void createMandatoryFields(fvMesh& mesh);

    // Original field validation 
    bool fieldExists(const fvMesh& mesh, const word& fieldName);
    bool validateField
    (
        const volScalarField& field,
        const word& fieldName,
        const dimensionSet& expectedDims
    );

    // Create additional fields after mandatory ones exist
    void createAdditionalFields(fvMesh& mesh);
}

#endif
