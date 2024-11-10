#ifndef advancedInterfaceCapturing_H
#define advancedInterfaceCapturing_H

#include "fvCFD.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "CMULES.H"
#include "MULES.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

namespace Foam
{

class advancedInterfaceCapturing
{
private:
    // References to mesh and fields
    const fvMesh& mesh_;
    volScalarField& alpha1_;
    const surfaceScalarField& phi_;
    const immiscibleIncompressibleTwoPhaseMixture& mixture_;

    // PLIC parameters
    scalar interfaceThickness_;
    scalar compressionCoeff_;

public:
    // Constructor
    advancedInterfaceCapturing
    (
        const fvMesh& mesh,
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const immiscibleIncompressibleTwoPhaseMixture& mixture
    );

    //- Correct the interface
    void correct();

    //- Return compression coefficient
    scalar compressionCoeff() const
    {
        return compressionCoeff_;
    }

private:
    //- Calculate interface normal
    tmp<volVectorField> interfaceNormal() const;
    
    //- Reconstruct interface using PLIC
    void reconstructInterface();
    
    //- Apply interface compression
    void applyInterfaceCompression();
};

}

#endif
