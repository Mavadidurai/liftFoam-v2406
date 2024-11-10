scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

if (mesh.nInternalFaces())
{
    // Use interpolate to get alpha1 at faces
    surfaceScalarField alphaPhi
    (
        "alphaPhi",
        fvc::interpolate(alpha1)*phi
    );

    // Calculate Courant numbers using primitiveField()
    scalarField alphaCo
    (
        (mesh.time().deltaT().value() * mag(alphaPhi.primitiveField()/mesh.magSf().primitiveField()))

    );

    // Get max and mean values
    alphaCoNum = 0.5*gMax(alphaCo);
    meanAlphaCoNum = 0.5*(gSum(alphaCo)/gSum(mesh.magSf().primitiveField()));
}

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;
