/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Class for modeling shock waves during femtosecond laser-induced forward 
    transfer (LIFT) process. Handles shock wave detection, propagation, and 
    effects in high-speed flows.

    The model:
    - Detects shock waves based on pressure gradients
    - Calculates shock strength and velocity
    - Simulates shock wave propagation
    - Accounts for local flow conditions
    
SourceFiles
    ultraFastShockWaveModel.C

\*---------------------------------------------------------------------------*/

#ifndef ultraFastShockWaveModel_H
#define ultraFastShockWaveModel_H

#include "fvCFD.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ultraFastShockWaveModel Declaration
\*---------------------------------------------------------------------------*/

class ultraFastShockWaveModel
{
private:
    // Private data
    
    //- Reference to mesh
    const fvMesh& mesh_;
    
    //- Reference to dictionary  
    const dictionary& dict_;

    //- Shock wave detection fields
    volScalarField shockStrength_;    //- Field for shock wave intensity
    volVectorField shockVelocity_;    //- Field for shock wave velocity vector
    
    //- Model parameters
    scalar shockThreshold_;    //- Minimum pressure gradient for shock detection
    scalar machNumber_;        //- Local flow Mach number
    scalar pressureRatio_;     //- Pressure ratio across shock
    scalar speedOfSound_;      //- Local speed of sound
    scalar dampingCoeff_;      //- Shock wave damping coefficient

    //- Private Member Functions
    
        //- Calculate shock wave properties
        void calculateShockProperties();
        
        //- Update shock velocity field
        void updateShockVelocity(const volVectorField& U);
        
        //- Calculate local speed of sound
        scalar calculateLocalSpeedOfSound
        (
            const scalar p,    //- Local pressure
            const scalar rho   //- Local density
        ) const;

public:
    // Constructors

        //- Construct from mesh and dictionary
        ultraFastShockWaveModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        ~ultraFastShockWaveModel() = default;

        //- Prevent copying and assignment 
        ultraFastShockWaveModel(const ultraFastShockWaveModel&) = delete;
        void operator=(const ultraFastShockWaveModel&) = delete;

    // Member Functions

        //- Update shock wave model
        void update
        (
            const volScalarField& p,  //- Pressure field
            const volVectorField& U   //- Velocity field
        );

        //- Return shock strength field
        tmp<volScalarField> shockStrength() const;

        //- Return velocity contribution from shock
        tmp<volVectorField> velocityContribution() const;

        //- Propagate shock wave
        void propagate();

        //- Check if model is valid
        bool valid() const;

        //- Write model data
        void write() const;
};

} // End namespace Foam

#endif
