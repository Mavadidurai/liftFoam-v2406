/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of the ultraFastShockWaveModel class for femtosecond LIFT
    process simulation. Contains methods for shock wave detection, propagation,
    and analysis.
\*---------------------------------------------------------------------------*/

#include "ultraFastShockWaveModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    ultraFastShockWaveModel Implementation
\*---------------------------------------------------------------------------*/

// Constructor
// Creates shock wave model with specified mesh and parameters
ultraFastShockWaveModel::ultraFastShockWaveModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    shockStrength_
    (
        IOobject("shockStrength", mesh.time().timeName(), mesh, 
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("shockStrength", dimless, 0)
    ),
    shockVelocity_
    (
        IOobject("shockVelocity", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    shockThreshold_(dict.lookupOrDefault<scalar>("shockThreshold", 1e6)),
    machNumber_(dict.lookupOrDefault<scalar>("machNumber", 1.0)),
    pressureRatio_(dict.lookupOrDefault<scalar>("pressureRatio", 2.0)),
    speedOfSound_(dict.lookupOrDefault<scalar>("speedOfSound", 340.0)),
    dampingCoeff_(dict.lookupOrDefault<scalar>("dampingCoeff", 0.1))
{
    calculateShockProperties();
}

// Calculate shock wave properties based on current flow state
void ultraFastShockWaveModel::calculateShockProperties()
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    forAll(mesh_.C(), cellI)
    {
        speedOfSound_ = calculateLocalSpeedOfSound(p[cellI], rho[cellI]);
    }
}

// Calculate local speed of sound using isentropic relation
scalar ultraFastShockWaveModel::calculateLocalSpeedOfSound
(
    const scalar p,
    const scalar rho
) const
{
    const scalar gamma = 1.4; // Specific heat ratio
    return sqrt(gamma * p / rho);
}

// Update shock wave model with current pressure and velocity fields
void ultraFastShockWaveModel::update
(
    const volScalarField& p,
    const volVectorField& U
)
{
    volScalarField gradP = mag(fvc::grad(p));
    volScalarField divU = fvc::div(U);

    forAll(mesh_.C(), cellI)
    {
        if (gradP[cellI] > shockThreshold_ && divU[cellI] < 0)
        {
            scalar pressureJump = gradP[cellI] / shockThreshold_;
            scalar localMach = mag(U[cellI])/speedOfSound_;
            
            shockStrength_[cellI] = pressureJump * 
                (1.0 + dampingCoeff_ * (localMach - machNumber_));
        }
        else
        {
            shockStrength_[cellI] = 0;
        }
    }

    updateShockVelocity(U);
}

// Update shock wave velocity based on local flow conditions
void ultraFastShockWaveModel::updateShockVelocity(const volVectorField& U)
{
    forAll(mesh_.C(), cellI)
    {
        if (shockStrength_[cellI] > SMALL)
        {
            scalar shockSpeed = machNumber_ * speedOfSound_ * 
                              sqrt(1.0 + (pressureRatio_ - 1.0)/2.0);
            
            vector flowDirection = U[cellI]/(mag(U[cellI]) + SMALL);
            shockVelocity_[cellI] = shockSpeed * flowDirection;
        }
        else
        {
            shockVelocity_[cellI] = vector::zero;
        }
    }
}

// Propagate shock wave and apply damping effects
void ultraFastShockWaveModel::propagate()
{
    const scalar deltaT = mesh_.time().deltaTValue();
    
    forAll(mesh_.C(), cellI)
    {
        if (shockStrength_[cellI] > SMALL)
        {
            shockStrength_[cellI] *= exp(-dampingCoeff_ * deltaT);
        }
    }
}

// Return shock strength field
tmp<volScalarField> ultraFastShockWaveModel::shockStrength() const
{
    return shockStrength_;
}

// Return shock velocity field
tmp<volVectorField> ultraFastShockWaveModel::velocityContribution() const
{
    return shockVelocity_;
}

// Validate model parameters and state
bool ultraFastShockWaveModel::valid() const
{
    return shockStrength_.size() > 0 && 
           machNumber_ > 1.0 && 
           pressureRatio_ > 1.0 &&
           speedOfSound_ > SMALL;
}

// Write model information to output
void ultraFastShockWaveModel::write() const
{
    Info<< "Ultra-fast shock wave model:" << nl
        << "  Threshold: " << shockThreshold_ << nl
        << "  Mach number: " << machNumber_ << nl
        << "  Pressure ratio: " << pressureRatio_ << nl
        << "  Speed of sound: " << speedOfSound_ << nl
        << "  Maximum shock strength: " << max(shockStrength_).value() << nl
        << "  Maximum shock velocity: " << max(mag(shockVelocity_)).value() 
        << endl;
}

} // End namespace Foam
