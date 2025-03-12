// Calculate the reciprocal of the momentum equation diagonal
volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rAU", dimTime/dimDensity, 1.0)
);

// Check if there's an explicit flux correction specified in the PIMPLE dict
bool correctPhiInit = pimple.dict().lookupOrDefault<Switch>("correctPhi", false);

if (correctPhiInit)
{
    // Explicit correction of the flux
    CorrectPhi
    (
        U,
        phi,
        p_rgh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU)),
        geometricZeroField(),
        pimple
    );

    // Check continuity errors
    #include "continuityErrs.H"
}
