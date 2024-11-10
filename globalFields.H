#ifndef globalFields_H
#define globalFields_H

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "turbulentTransportModel.H"
#include "incompressibleTurbulenceModel.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

namespace Foam
{
    class GlobalFields
    {
    private:
        const fvMesh& mesh_;
        bool initialized_;

    public:
        // Field storage
        autoPtr<volVectorField> UPtr;
        autoPtr<surfaceScalarField> phiPtr;
        autoPtr<surfaceScalarField> rhoPhiPtr;
        autoPtr<volScalarField> pPtr;
        autoPtr<volScalarField> p_rghPtr;
        autoPtr<volScalarField> rhoPtr;
        autoPtr<volScalarField> alpha1Ptr;
        autoPtr<volScalarField> alpha2Ptr;
        autoPtr<volScalarField> TePtr;
        autoPtr<volScalarField> TlPtr;
        autoPtr<volScalarField> ghPtr;
        autoPtr<surfaceScalarField> ghfPtr;
        autoPtr<uniformDimensionedVectorField> gPtr;

        // Reference values
        label pRefCell;
        scalar pRefValue;

        // Models
        autoPtr<immiscibleIncompressibleTwoPhaseMixture> mixturePtr;
        autoPtr<incompressible::turbulenceModel> turbulencePtr;

        // Constructor
        GlobalFields(const fvMesh& mesh)
        :
            mesh_(mesh),
            initialized_(false),
            pRefCell(0),
            pRefValue(0.0)
        {}

        // Access functions
        const volVectorField& U() const { return UPtr(); }
        const surfaceScalarField& phi() const { return phiPtr(); }
        const volScalarField& p() const { return pPtr(); }
        const volScalarField& p_rgh() const { return p_rghPtr(); }
        const volScalarField& rho() const { return rhoPtr(); }
        const volScalarField& alpha1() const { return alpha1Ptr(); }
        const volScalarField& alpha2() const { return alpha2Ptr(); }
        const volScalarField& Te() const { return TePtr(); }
        const volScalarField& Tl() const { return TlPtr(); }

        // Field management
        bool initialize();
        bool validate();
        void clear();
    };

    extern GlobalFields* globalFieldsPtr;
}

#endif
