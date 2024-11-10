#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "fvm.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

namespace Foam
{

advancedInterfaceCapturing::advancedInterfaceCapturing
(
    const fvMesh& mesh,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    phi_(phi),
    mixture_(mixture),
    interfaceThickness_(mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("interfaceThickness", 1e-6)),
    compressionCoeff_(mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("compressionCoeff", 1.0))
{}

tmp<volVectorField> advancedInterfaceCapturing::interfaceNormal() const
{
    tmp<volVectorField> tInterfaceNormal
    (
        new volVectorField
        (
            IOobject
            (
                "interfaceNormal",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless, vector::zero)
        )
    );
    
    volVectorField& interfaceNormal = tInterfaceNormal.ref();

    tmp<volVectorField> gradAlpha = fvc::grad(alpha1_);
    const volVectorField& gradAlphaRef = gradAlpha();

    forAll(alpha1_, cellI)
    {
        if (alpha1_[cellI] > SMALL && alpha1_[cellI] < (1.0 - SMALL))
        {
            vector& n = interfaceNormal[cellI];
            n = gradAlphaRef[cellI];
            scalar magN = mag(n);

            if (magN > SMALL)
            {
                n /= magN;
            }
        }
    }

    return tInterfaceNormal;
}

void advancedInterfaceCapturing::reconstructInterface()
{
    tmp<volVectorField> tNormal = interfaceNormal();
    const volVectorField& normal = tNormal();

    forAll(alpha1_, cellI)
    {
        if (alpha1_[cellI] > SMALL && alpha1_[cellI] < (1.0 - SMALL))
        {
            const point& cellCenter = mesh_.C()[cellI];
            scalar signedDistance = 
                (normal[cellI] & (cellCenter - mesh_.C()[cellI]));
            alpha1_[cellI] = 0.5 * (1.0 + tanh(signedDistance/interfaceThickness_));
        }
    }
}

void advancedInterfaceCapturing::applyInterfaceCompression()
{
    // Calculate interface compression flux
    surfaceScalarField phic = mag(phi_/mesh_.magSf());
    phic = min(compressionCoeff_*phic, max(phic));
    
    // Get face normal field
    surfaceScalarField phir(phic*mixture_.nHatf());

    // Create alpha flux
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        phi_*fvc::interpolate(alpha1_)
    );

    // Add compression term
    surfaceScalarField phiCorr
    (
        IOobject
        (
            "phiCorr",
            mesh_.time().timeName(),
            mesh_
        ),
        fvc::flux(-phir, scalar(1) - alpha1_, "div(phir,alpha)")
    );

    // Create bound fields for MULES
    scalarField& allAlpha = alpha1_.primitiveFieldRef();
    scalarField allAlphaMin(allAlpha.size(), 0.0);
    scalarField allAlphaMax(allAlpha.size(), 1.0);

    // Apply interface compression
    MULES::limit
    (
        1.0/mesh_.time().deltaTValue(),
        geometricOneField(),
        alpha1_,
        phi_,
        phiCorr,
        zeroField(),
        zeroField(),
        allAlphaMin,
        allAlphaMax,
        false
    );

    // Update alpha1
    alpha1_ = max(alpha1_, scalar(0));
    alpha1_ = min(alpha1_, scalar(1));

    // Update boundary conditions
    alpha1_.correctBoundaryConditions();

    Info<< "Phase-1 volume fraction = "
        << alpha1_.weightedAverage(mesh_.V()).value()
        << "  Min(" << alpha1_.name() << ") = " << min(alpha1_).value()
        << "  Max(" << alpha1_.name() << ") = " << max(alpha1_).value()
        << endl;
}

void advancedInterfaceCapturing::correct()
{
    // Step 1: Interface reconstruction using PLIC
    reconstructInterface();

    // Step 2: Apply interface compression
    applyInterfaceCompression();

    Info<< "Interface capturing: "
        << "min(alpha1) = " << min(alpha1_).value()
        << ", max(alpha1) = " << max(alpha1_).value()
        << endl;
}

}
