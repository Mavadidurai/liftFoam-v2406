/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Two-temperature model for femtosecond laser-material interaction in LIFT process.
    Models separate electron and lattice temperatures, including:
    - Electron-phonon coupling
    - Temperature-dependent thermal properties
    - Energy conservation tracking
    - Laser energy absorption
    
    The model solves coupled equations for:
    - Electron temperature (Te)
    - Lattice temperature (Tl)
    
    With temperature-dependent properties:
    - Electron thermal conductivity (ke)
    - Lattice thermal conductivity (kl)
    - Electron-phonon coupling coefficient (G)
    
SourceFiles
    twoTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoTemperatureModel_H
#define twoTemperatureModel_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "DimensionValidator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class twoTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class twoTemperatureModel
{
private:
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to dictionary
        const dictionary& dict_;

        //- Temperature Fields
        volScalarField Te_;        //- Electron temperature [K]
        volScalarField Tl_;        //- Lattice temperature [K]

        //- Material Properties
        dimensionedScalar Ce_;     //- Electron heat capacity [J/m³/K]
        dimensionedScalar Cl_;     //- Lattice heat capacity [J/m³/K]
        dimensionedScalar G_;      //- Electron-phonon coupling [W/m³/K]

        //- Energy Conservation Tracking
        mutable dimensionedScalar lastTotalEnergy_;  //- Previous total energy [J]
        mutable bool energyInitialized_;             //- Energy tracking flag

    // Private Member Functions

        //- Validate model parameters and dimensions
        bool validateParameters() const;

        //- Validate temperature field values
        bool validateFields() const;

        //- Check energy conservation
        bool checkEnergyConservation() const;

        //- Update energy tracking variables
        void updateEnergyTracking() const;

public:
    // Constructors

        //- Construct from mesh and dictionary
        twoTemperatureModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~twoTemperatureModel();

        //- Prevent copying
        twoTemperatureModel(const twoTemperatureModel&) = delete;
        void operator=(const twoTemperatureModel&) = delete;

    // Member Functions

        //- Solve temperature equations
        void solve(const volScalarField& laserSource);

        //- Calculate electron thermal conductivity
        tmp<volScalarField> electronThermalConductivity() const;

        //- Calculate electron heat capacity
        tmp<volScalarField> electronHeatCapacity() const;

        //- Calculate electron-phonon coupling
        tmp<volScalarField> electronPhononCoupling() const;

        //- Check model validity
        bool valid() const;

        //- Write model information
        void write() const;

        //- Calculate thermal conductivities
        tmp<volScalarField> ke() const;  //- Electronic thermal conductivity
        tmp<volScalarField> kl() const;  //- Lattice thermal conductivity

        // Access Functions

        //- Return electron temperature
        const volScalarField& Te() const { return Te_; }

        //- Return lattice temperature
        const volScalarField& Tl() const { return Tl_; }

        //- Return electron heat capacity
        const dimensionedScalar& Ce() const { return Ce_; }

        //- Return lattice heat capacity
        const dimensionedScalar& Cl() const { return Cl_; }

        //- Return electron-phonon coupling
        const dimensionedScalar& G() const { return G_; }
};

} // End namespace Foam

#endif
