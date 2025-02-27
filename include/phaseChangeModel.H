/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Phase change model for LIFT process simulation, handling melting and 
    vaporization transitions during femtosecond laser heating.

    Features:
    - Equilibrium and non-equilibrium phase change
    - Temperature and pressure dependent melting
    - Residual stress calculation
    - Energy conservation tracking
    - Interface tracking
    - Source terms for momentum, energy and phase equations

    The model accounts for:
    - Melting/solidification
    - Vaporization (optional)
    - Thermal stress
    - Interface dynamics
    - Energy conservation
    
SourceFiles
    phaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "fvCFD.H"
#include "fvMatrix.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "symmTensor.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
protected:
    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to phase fraction field
        volScalarField& alpha1_;

        //- Reference to dictionary
        const dictionary& dict_;

        //- Pointer to mixture model
        const immiscibleIncompressibleTwoPhaseMixture* mixture_;

        // Material Properties

            //- Melting temperature [K]
            dimensionedScalar Tm_;

            //- Vaporization temperature [K]
            dimensionedScalar Tv_;

            //- Latent heat of melting [J/kg]
            dimensionedScalar Lm_;

            //- Latent heat of vaporization [J/kg]
            dimensionedScalar Lv_;

            //- Interface width [m]
            dimensionedScalar interfaceWidth_;

            //- Temperature-pressure coupling [K/Pa]
            dimensionedScalar dTmdp_;

            //- Reference pressure [Pa]
            dimensionedScalar pRef_;

        // Residual Stress Properties

            //- Young's modulus [Pa]
            scalar E_;

            //- Poisson's ratio [-]
            scalar nu_;

            //- Thermal expansion coefficient [1/K]
            scalar alpha_;

            //- Reference temperature [K]
            scalar Tref_;

        // Kinetic Parameters

            //- Gradual melting rate [1/s]
            dimensionedScalar gradualMeltingRate_;

            //- Gradual solidification rate [1/s]
            dimensionedScalar gradualSolidificationRate_;

            //- Undercooling coefficient [K/m]
            dimensionedScalar undercoolingCoeff_;

        // Fields

            //- Phase indicator field
            volScalarField phaseIndicator_;

            //- Density field [kg/m³]
            volScalarField rho_;

            //- Phase change rate field [J/m³/s]
            mutable tmp<volScalarField> phaseChangeRate_;

        // Energy Tracking

            //- Previous total energy [J]
            mutable dimensionedScalar lastTotalEnergy_;

            //- Energy tracking initialization flag
            mutable bool energyInitialized_;

        // Model Configuration

            //- Non-equilibrium phase change flag
            bool isNonEquilibrium_;

    // Protected Member Functions

        //- Calculate phase change rate
        virtual void calculatePhaseChangeRate() const;

        //- Update density field
        virtual void updateDensity();

        //- Validate model configuration
        virtual bool validate() const;

        //- Check energy conservation
        virtual bool checkEnergyConservation() const;

        //- Update energy tracking variables
        virtual void updateEnergyTracking() const;
   // Helper functions for numerical stability and validation
      // Numerical stability parameters
        
        //- Maximum allowed phase change rate [J/m³/s]
        dimensionedScalar maxPhaseChangeRate_;
        
        //- Temperature transition width for smooth phase change [K]
        dimensionedScalar transitionWidth_;
        
        //- Maximum allowed temperature [K]
        dimensionedScalar maxTemperature_;
        
        //- Maximum allowed stress [Pa]
        dimensionedScalar maxStress_;
        
        //- Maximum allowed temperature gradient [K/m]
        dimensionedScalar maxTemperatureGradient_;
        
        //- Energy conservation tolerance [-]
        scalar energyTolerance_;
    //- Check and constrain temperature values
    void checkAndLimitTemperature(scalar& T) const;
    
    //- Validate material property bounds
    bool checkPropertyBounds() const;
    
    //- Limit value between bounds (inline helper)
    inline scalar limitValue
    (
        const scalar value,
        const scalar minVal,
        const scalar maxVal
    ) const
    {
        return min(max(value, minVal), maxVal);
    }

    //- Helper for energy calculations
    inline scalar calculateLocalDensity(const label cellI) const;
    
    //- Helper for phase transition calculations 
    inline scalar calculateEffectiveMeltTemp
    (
        const scalar T,
        const scalar gradMag
    ) const;

public:
    TypeName("phaseChangeModel");

    // Constructors
    phaseChangeModel
    (
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict,
        bool isNonEquilibrium = false
    );

    static autoPtr<phaseChangeModel> New
    (
        const word& modelType,
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    );

    //- Destructor 
    virtual ~phaseChangeModel() = default;

    // Delete default copy/assign
    phaseChangeModel(const phaseChangeModel&) = delete;
    void operator=(const phaseChangeModel&) = delete;

    // Main functions
    virtual void correct(const volScalarField& T);
    virtual void correctPressure(const volScalarField& p);
    
    // Phase change source terms
    virtual tmp<volScalarField> source() const;
    virtual tmp<fvVectorMatrix> momentumSource() const;
    virtual tmp<volScalarField> volumetricForce() const;
    virtual tmp<volScalarField> recoilPressure() const;
    virtual tmp<volScalarField> massSource() const;
    virtual tmp<fvScalarMatrix> energySource(const volScalarField& h) const;
    virtual tmp<fvScalarMatrix> electronSource() const;
    virtual tmp<fvScalarMatrix> latticeSource() const;

    // Residual stress calculations
    tmp<volSymmTensorField> residualStress() const;

    // Helper functions
    scalar liquidFraction() const;
    scalar solidFraction() const;
    scalar vaporFraction() const;
    scalar interfaceArea() const;
    tmp<volScalarField> subcooling() const;
    dimensionedScalar totalPhaseChangeEnergy() const;

    // IO
    virtual bool valid() const;
    virtual void write() const;

    // Access functions
    const volScalarField& phaseIndicator() const { return phaseIndicator_; }
    const dimensionedScalar& meltingTemp() const { return Tm_; }
    const dimensionedScalar& vaporizationTemp() const { return Tv_; }
    bool isNonEquilibrium() const { return isNonEquilibrium_; }
    scalar youngsModulus() const { return E_; }
    scalar poissonsRatio() const { return nu_; }
    scalar thermalExpansion() const { return alpha_; }
    scalar referenceTemp() const { return Tref_; }
};

} // End namespace Foam

#endif
