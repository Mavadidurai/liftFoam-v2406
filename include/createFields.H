// createFields.H

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Reading gravitational acceleration
Info<< "Reading gravitational acceleration\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    uniformDimensionedVectorField
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    )
);

// Read field p_rgh
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rgh", dimPressure, 0.0)
);

// Read field U
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector::zero)
);

#include "createPhi.H"

// Create two-phase mixture model
Info<< "Creating two-phase mixture model\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

// Get references to phase fractions
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Create density field for ddt(rho, U)
Info<< "Creating field rho\n" << endl;
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
Info<< "Creating field rhoPhi\n" << endl;
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
);

// Temperature field
Info<< "Reading/creating field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T0", dimTemperature, 300.0)
);

// Femtosecond temperature fields
Info<< "Reading/creating electron temperature field Te\n" << endl;
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
    dimensionedScalar("Te0", dimTemperature, 300.0)
);

Info<< "Reading/creating lattice temperature field Tl\n" << endl;
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
    dimensionedScalar("Tl0", dimTemperature, 300.0)
);

// Laser source field
Info<< "Creating laser source field\n" << endl;
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
    dimensionedScalar("zero", dimPower/dimVolume, 0.0)
);

// Phase change rate field
Info<< "Creating phase change rate field\n" << endl;
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
    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
);

// Reference height for hydrostatic pressure
dimensionedScalar hRef
(
    "hRef",
    dimLength,
    mixture.lookupOrDefault<scalar>("hRef", 0.0)
);

// Hydrostatic pressure reference point
dimensionedScalar ghRef
(
    "ghRef",
    g.dimensions()*dimLength,
    0
);

// Hydrostatic pressure at cell centers and faces
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

// Create absolute pressure
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

// Surface tension and interface compression
dimensionedScalar sigma
(
    "sigma",
    dimForce/dimLength,
    transportProperties.lookupOrDefault<scalar>("sigma", 0.07)
);

// Create MRF
#include "createMRF.H"

// Create turbulence model
Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

// Read fvOptions
Info<< "Reading fvOptions\n" << endl;
fv::options& fvOptions(fv::options::New(mesh));

// Thermal properties
dimensionedScalar thermalDiffusivity
(
    "thermalDiffusivity",
    dimArea/dimTime,
    transportProperties.lookupOrDefault<scalar>("thermalDiffusivity", 1e-5)
);

// Interface compression coefficient
dimensionedScalar cAlphaCoeff
(
    "cAlpha",
    dimless,
    transportProperties.lookupOrDefault<scalar>("cAlpha", 1.0)
);

// Create volume-weighted flux needed in the PISO loop
surfaceScalarField phig
(
    IOobject
    (
        "phig",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    ghf*fvc::snGrad(rho)*mesh.magSf()
);

// Create volume-weighted mass flux
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh
    ),
    phi*fvc::interpolate(alpha1)
);

// Make sure all fluxes and fields are stored in the registry
mesh.objectRegistry::store(new surfaceScalarField(phi));
