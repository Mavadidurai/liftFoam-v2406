#ifndef ALPHAEQN_H
#define ALPHAEQN_H

// alphaEqn.H

{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("vanLeer");

    // Store old time alpha1
    volScalarField alpha10("alpha10", alpha1);

    // Calculate interface compression
    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(cAlpha*phic, max(phic));
    
    surfaceScalarField phir(phic*mixture.nHatf());

    // Initialize limiters
    volScalarField Sp
    (
        IOobject
        (
            "Sp",
            runTime.timeName(), 
            mesh
        ),
        mesh,
        dimensionedScalar("Sp", dimless/dimTime, 0.0)
    );

    volScalarField Su
    (
        IOobject
        (
            "Su", 
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("Su", dimless/dimTime, 0.0)
    );

    // Add phase change source if present
    if (mesh.foundObject<volScalarField>("phaseChangeRate"))
    {
        const volScalarField& pcr = 
            mesh.lookupObject<volScalarField>("phaseChangeRate");
        
        if (max(mag(pcr)).value() > SMALL)
        {
            // Create volume field
            scalarField& SpI = Sp.primitiveFieldRef();
            scalarField divPhiI = fvc::div(phi)().primitiveField();
            
            forAll(SpI, cellI)
            {
                SpI[cellI] -= divPhiI[cellI]/mesh.V()[cellI];
            }
            
            // Create source terms for Su
            scalarField& SuI = Su.primitiveFieldRef();
            const scalarField& pcrI = pcr.primitiveField();
            const scalar rho1Inv = 1.0/mixture.rho1().value();
            const scalar rho2Inv = 1.0/mixture.rho2().value();
            const scalar deltaRhoInv = rho1Inv - rho2Inv;
            
            forAll(SuI, cellI)
            {
                SuI[cellI] += deltaRhoInv*pcrI[cellI]/mesh.V()[cellI];
            }
        }
    }

    // Calculate alpha flux
    surfaceScalarField alphaPhi
    (
        "alphaPhi",
        fvc::flux(phi, alpha1, alphaScheme)
    );

    // Add interface compression
    alphaPhi += fvc::flux(-phir*fvc::snGrad(alpha1)*mesh.magSf(), alpha1, alpharScheme);

    // Apply MULES limiter with fixed bounds instead of scalar values
    scalarField psiMax(alpha1.size(), 1.0);
    scalarField psiMin(alpha1.size(), 0.0);
    
    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        phi,
        alphaPhi,
        Sp,
        Su,
        psiMax,
        psiMin
    );

    // Update alpha2
    alpha2 = scalar(1) - alpha1;

    // Apply bounds
    alpha1.max(0.0);
    alpha1.min(1.0);
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();

    // Update mixture properties
    mixture.correct();

    // Update density and mass flux
    rho = alpha1*rho1 + alpha2*rho2;
    rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

#endif
