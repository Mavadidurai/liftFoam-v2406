#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "wallFvPatch.H"
#include "fvPatchField.H"

namespace Foam
{

advancedInterfaceCapturing::advancedInterfaceCapturing
(
    const fvMesh& mesh,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    const volScalarField& T
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    phi_(phi),
    mixture_(mixture),
    T_(T),
    compressionCoeff_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("compressionCoeff", 1.0)
    ),
    interfaceThickness_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("interfaceThickness", 1e-6)
    ),
    baseContactAngle_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("contactAngle", 90.0)
    ),
    meltingTemp_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("meltingTemperature", 1941.0)  // Titanium melting point
    ),
    vaporTemp_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("vaporTemperature", 3560.0)  // Titanium boiling point
    ),
    latentHeat_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("latentHeat", 3.96e5)  // Titanium latent heat
    ),
    surfaceTensionCoeff_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("surfaceTension", 1.65)  // Titanium at melting point
    ),
    surfaceTensionGradient_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("surfaceTensionGradient", -0.26e-3)  // dσ/dT
    ),
    maxSurfaceTension_(2.0),  // Upper limit
    minSurfaceTension_(0.1),  // Lower limit
    interfaceNormal_
    (
        IOobject
        (
            "interfaceNormal",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    curvature_
    (
        IOobject
        (
            "curvature",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    ),
    phaseChangeRate_
    (
        IOobject
        (
            "phaseChangeRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
    ),
    ejectionVelocity_
    (
        IOobject
        (
            "ejectionVelocity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    )
{
    // Initialize interface tracking
    calculateInterfaceNormal();
    calculateCurvature();
}

void advancedInterfaceCapturing::calculateInterfaceNormal()
{
    tmp<volVectorField> gradAlpha = fvc::grad(alpha1_);
    const volScalarField magGradAlpha(mag(gradAlpha));

    forAll(mesh_.C(), cellI)
    {
        if (magGradAlpha[cellI] > SMALL)
        {
            interfaceNormal_[cellI] = gradAlpha()[cellI]/magGradAlpha[cellI];
        }
    }

    interfaceNormal_.correctBoundaryConditions();
}
void Foam::advancedInterfaceCapturing::calculateCurvature()
{
    // Calculate interface curvature using simpler formulation from curvatureModel
    volVectorField gradAlpha = fvc::grad(alpha1_);
    volScalarField magGradAlpha = mag(gradAlpha);
    
    curvature_ = -fvc::div(gradAlpha/(magGradAlpha + SMALL));
    
    // Apply bounds to curvature field
    forAll(mesh_.C(), cellI)
    {
        if (mag(curvature_[cellI]) > 1.0/interfaceThickness_)
        {
            curvature_[cellI] *= interfaceThickness_*mag(curvature_[cellI]);
        }
    }

    curvature_.correctBoundaryConditions();
}

void Foam::advancedInterfaceCapturing::applyContactAngle()
{
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))  // Changed from boundaryMesh to boundary
        {
            vectorField nf = patch.nf();
            
            forAll(patch, faceI)
            {
                vector nw = nf[faceI];
                scalar theta = baseContactAngle_*constant::mathematical::pi/180.0;
                
                nw = (interfaceNormal_.boundaryField()[patchI][faceI] -
                      cos(theta)*nf[faceI])/sin(theta);
                      
                nw /= (mag(nw) + SMALL);
                
                alpha1_.boundaryFieldRef()[patchI][faceI] = 
                    0.5*(1.0 + sin(theta));
            }
        }
    }
}

tmp<surfaceScalarField> Foam::advancedInterfaceCapturing::computeCompressionFlux() const
{
    // Calculate interface compression flux
    const surfaceScalarField deltaN
    (
        "deltaN",
        phi_/(mag(phi_) + SMALL)*
        mixture_.nHatf()
    );

    tmp<surfaceScalarField> tphiC
    (
        new surfaceScalarField
        (
            "phiC",
            deltaN*
            (
                min
                (
                    compressionCoeff_*deltaN,
                    max(deltaN)
                )
            )*
            mag(phi_)
        )
    );

    // Scale by interface indicator
    surfaceScalarField& phiC = tphiC.ref();
    
    // Get interpolated alpha1 field
    tmp<surfaceScalarField> talpha1f = fvc::interpolate(alpha1_);
    const surfaceScalarField& alpha1f = talpha1f();

    // Access internal field directly
    scalarField& phiCInt = phiC.primitiveFieldRef();
    const scalarField& deltaNInt = deltaN.primitiveField();
    const scalarField& alpha1fInt = alpha1f.primitiveField();

    forAll(phiCInt, faceI)
    {
        if (mag(deltaNInt[faceI]) > SMALL)
        {
            scalar alpha = alpha1fInt[faceI];
            if (alpha > SMALL && alpha < 1.0 - SMALL)
            {
                phiCInt[faceI] *= alpha*(1.0 - alpha);
            }
            else
            {
                phiCInt[faceI] = 0.0;
            }
        }
    }

    // Handle boundary field
    surfaceScalarField::Boundary& phiCBf = phiC.boundaryFieldRef();
    const surfaceScalarField::Boundary& deltaNBf = deltaN.boundaryField();
    const surfaceScalarField::Boundary& alpha1fBf = alpha1f.boundaryField();

    forAll(phiCBf, patchI)
    {
        fvsPatchScalarField& patchPhiC = phiCBf[patchI];
        const fvsPatchScalarField& patchDeltaN = deltaNBf[patchI];
        const fvsPatchScalarField& patchAlpha = alpha1fBf[patchI];
        
        forAll(patchPhiC, faceI)
        {
            if (mag(patchDeltaN[faceI]) > SMALL)
            {
                scalar alpha = patchAlpha[faceI];
                if (alpha > SMALL && alpha < 1.0 - SMALL)
                {
                    patchPhiC[faceI] *= alpha*(1.0 - alpha);
                }
                else
                {
                    patchPhiC[faceI] = 0.0;
                }
            }
        }
    }

    return tphiC;
}
void advancedInterfaceCapturing::calculatePhaseChange()
{
    forAll(mesh_.C(), cellI)
    {
        const scalar T = T_[cellI];
        if (T > meltingTemp_ && alpha1_[cellI] > SMALL)
        {
            // Calculate phase change rate based on temperature
            scalar overHeat = (T - meltingTemp_)/meltingTemp_;
            phaseChangeRate_[cellI] = latentHeat_ * overHeat;
            
            // Modify interface for phase change
            if (T > vaporTemp_)
            {
                alpha1_[cellI] *= max(0.0, 1.0 - overHeat);
            }
        }
    }
}

void advancedInterfaceCapturing::calculateEjectionVelocity()
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    
    forAll(mesh_.C(), cellI)
    {
        if (alpha1_[cellI] > SMALL && T_[cellI] > meltingTemp_)
        {
            // Calculate ejection velocity based on pressure and surface tension
            scalar sigma = calculateLocalSurfaceTension(T_[cellI]);
            scalar pDriving = p[cellI] + sigma*curvature_[cellI];
            
            // Simplified ejection velocity model
            ejectionVelocity_[cellI] = sqrt(2.0*max(pDriving, 0.0)/mixture_.rho1().value());
        }
    }
}

scalar advancedInterfaceCapturing::calculateLocalSurfaceTension
(
    const scalar T
) const
{
    // Temperature-dependent surface tension with limits
    scalar sigma = surfaceTensionCoeff_ + 
                  surfaceTensionGradient_*(T - meltingTemp_);
    
    return max(min(sigma, maxSurfaceTension_), minSurfaceTension_);
}

void advancedInterfaceCapturing::updateTemperatureDependentProps()
{
    forAll(mesh_.C(), cellI)
    {
        // Update contact angle based on local temperature
        if (T_[cellI] > meltingTemp_)
        {
            // Reduce contact angle when material is molten
            scalar tempRatio = min((T_[cellI] - meltingTemp_)/(vaporTemp_ - meltingTemp_), 1.0);
            // Contact angle tends toward 180° (detachment) at high temperature
            baseContactAngle_ = 90.0 + 90.0*tempRatio;
        }
    }
}

void advancedInterfaceCapturing::correct()
{
    // 1. Update temperature-dependent properties
    updateTemperatureDependentProps();

    // 2. Calculate interface geometry
    calculateInterfaceNormal();
    calculateCurvature();
    applyContactAngle();

    // 3. Handle phase change and ejection
    calculatePhaseChange();
    calculateEjectionVelocity();

    // 4. Calculate compression flux
    tmp<surfaceScalarField> tPhiC = computeCompressionFlux();
    tmp<volVectorField> tGradAlpha = fvc::grad(alpha1_);
        
        // Force cleanup when not needed
        if (mesh_.time().timeIndex() % 100 == 0)
        {
            alpha1_.clearOldTimes();
        }
    // 5. Create phase flux with phase change contribution
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        phi_*fvc::interpolate(alpha1_) + 
        fvc::interpolate(phaseChangeRate_/mixture_.rho1())
    );

    // 6. Apply MULES solver with compression
    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1_,
        phi_,
        phiAlpha,
        zeroField(),
        zeroField(),
        scalarField(alpha1_.size(), 0),
        scalarField(alpha1_.size(), 1)
    );

    // 7. Final boundedness check
    alpha1_.max(0.0);
    alpha1_.min(1.0);

    Info<< "Phase-1 volume fraction = "
        << alpha1_.weightedAverage(mesh_.V()).value()
        << "  Min(alpha1) = " << min(alpha1_).value()
        << "  Max(alpha1) = " << max(alpha1_).value()
        << "\nMax temperature = " << max(T_).value()
        << "  Max ejection velocity = " << max(ejectionVelocity_).value()
        << endl;
}

void advancedInterfaceCapturing::write() const
{
    interfaceNormal_.write();
    curvature_.write();
    phaseChangeRate_.write();
    ejectionVelocity_.write();
}


} // End namespace Foam
