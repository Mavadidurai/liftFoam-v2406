#ifndef CREATE_FIELDS_H
#define CREATE_FIELDS_H

Info<< "Reading/creating fields for LIFT process\n" << endl;

// Declare reference cells and values
label pRefCell = 0;
scalar pRefValue = 0.0;

// Create velocity field
Info<< "Reading/creating velocity field U\n" << endl;
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("U", dimVelocity, vector(0, 0, 0))
);

// Create flux field
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    linearInterpolate(U) & mesh.Sf()
);

// Initialize two-phase mixture
Info<< "Setting up two-phase mixture\n" << endl;
autoPtr<immiscibleIncompressibleTwoPhaseMixture> mixturePtr
(
    new immiscibleIncompressibleTwoPhaseMixture(U, phi)
);
const immiscibleIncompressibleTwoPhaseMixture& mixture = mixturePtr();
mesh.objectRegistry::store(mixturePtr.ptr());

// Get phase fractions
volScalarField& alpha1 = const_cast<volScalarField&>(mixture.alpha1());
volScalarField& alpha2 = const_cast<volScalarField&>(mixture.alpha2());

// Create sigma field first to ensure proper registry
volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sigma", dimForce/dimLength, 1.64)  // Default surface tension
);

// Initialize phase properties
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Create density field
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Create mass flux field
surfaceScalarField rhoPhi
(
    IOobject("rhoPhi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(rho)*phi
);

// Read gravitational acceleration
#include "readGravitationalAcceleration.H"

// Create gravity-related fields
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

// Initialize pressure fields
volScalarField p_rgh
(
    IOobject("p_rgh", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("p_rgh", dimPressure, 0.0)
);

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    p_rgh + rho*gh
);

// Initialize temperature fields
volScalarField Te
(
    IOobject("Te", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Te", dimTemperature, 300.0)
);

volScalarField Tl
(
    IOobject("Tl", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Tl", dimTemperature, 300.0)
);

// Create laser and phase change fields
volScalarField laserSource
(
    IOobject("laserSource", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("laserSource", dimPower/dimVolume, 0)
);

volScalarField phaseChangeRate
(
    IOobject("phaseChangeRate", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("phaseChangeRate", dimMass/dimVolume/dimTime, 0)
);

// Create turbulence model
Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulencePtr =
    incompressible::turbulenceModel::New(U, phi, mixture);

// Set reference pressure cell
setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

#endif
