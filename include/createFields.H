#ifndef CREATE_FIELDS_H
#define CREATE_FIELDS_H

Info<< "Initializing fields\n" << endl;
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(), 
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);Info<< "Creating temperature fields\n" << endl;
volScalarField Te
(
    IOobject
    (
        "Te",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Te", dimTemperature, 300.0)
);

volScalarField Tl
(
    IOobject
    (
        "Tl", 
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Tl", dimTemperature, 300.0)
);

volScalarField laserSource
(
    IOobject
    (
        "laserSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPower/dimVolume, 0)
);

volScalarField phaseChangeRate
(
    IOobject
    (
        "phaseChangeRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
);
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

// MRF
#include "createMRF.H"

// Turbulence
turbulence = incompressible::turbulenceModel::New(U, phi, mixture);
// Initialize fvOptions
Info<< "Reading fvOptions\n" << endl;
fv::options& fvOptions(fv::options::New(mesh));

// Final validation
bool fieldsOK = 
    mesh.foundObject<volScalarField>("cellLevel") &&
    mesh.foundObject<volScalarField>("pointLevel") &&
    mesh.foundObject<volScalarField>("alpha.titanium") &&
    mesh.foundObject<volVectorField>("U") &&
    mesh.foundObject<surfaceScalarField>("phi") &&
    mesh.foundObject<volScalarField>("p") &&
    mesh.foundObject<volScalarField>("Te") &&
    mesh.foundObject<volScalarField>("Tl");

if (!fieldsOK)
{
    FatalErrorIn("createFields.H")
        << "Not all required fields were created"
        << abort(FatalError);
}

#endif
