/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Main solver class for Laser-Induced Forward Transfer (LIFT) process.
    
    Features:
    - Two-temperature model for electron-lattice coupling
    - Femtosecond laser pulse modeling
    - Phase change with melting and vaporization
    - Interface tracking with surface tension
    - Droplet formation and tracking
    - Adaptive mesh refinement
    - Energy conservation tracking
    
    The solver couples:
    - Thermal dynamics (electron/lattice temperatures)
    - Fluid dynamics (incompressible two-phase flow)
    - Phase change (melting, solidification, vaporization)
    - Surface tension and interface tracking
    - Laser-material interaction
    
    Reference:
    For details on the LIFT process and modeling approach, see:
    - Phipps et al. (2010) Review: Laser-Ablation Propulsion
    - Willis & Xu (2000) Transport Phenomena and Droplet Formation
    
SourceFiles
    LIFTModel.C

\*---------------------------------------------------------------------------*/

#ifndef LIFTModel_H
#define LIFTModel_H

#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOMRFZoneList.H"
#include "fvMatrix.H"
#include "radiationModel.H"
#include "basicThermo.H" 
#include "advancedInterfaceCapturing.H"
//#include "dynamicRefineFvMesh.H"

#include "twoTemperatureModel.H"  
#include "phaseChangeModel.H"     
#include "femtosecondLaserModel.H" 

// Forward declarations
namespace Foam
{
    class femtosecondLaserModel;
    class dropletModel;
}

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class LIFTModel Declaration
\*---------------------------------------------------------------------------*/

class LIFTModel
{
private:
    // References
    const fvMesh& mesh_;
    //autoPtr<dynamicRefineFvMesh> dynamicMesh_; 
    const dictionary& dict_;
    Time& runTime_;
    
    // Core fields
    volVectorField& U_;             // Velocity
    surfaceScalarField& phi_;       // Flux
    volScalarField& p_;            // Pressure
    volScalarField& alpha1_;       // Phase fraction
    volScalarField& rho_;          // Density
    volScalarField& Te_;           // Electron temperature
    volScalarField& Tl_;           // Lattice temperature
    
    // Laser and phase change fields
    volScalarField laserSource_;    // Laser energy deposition
    volScalarField phaseChangeRate_; // Phase change rate
    volScalarField phaseIndicator_; // Phase state indicator
    
    // Interface tracking fields
    volScalarField interfaceEnergy_;
    volVectorField interfaceNormal_;
    volScalarField curvature_;
    
    // Pressure-related fields
    volScalarField p_rgh_;         // Dynamic pressure
    surfaceScalarField rhoPhi_;    // Density flux
    volScalarField gh_;            // Gravity potential
    surfaceScalarField ghf_;       // Face gravity potential
    
    // Models
    autoPtr<immiscibleIncompressibleTwoPhaseMixture> mixturePtr_;
    immiscibleIncompressibleTwoPhaseMixture& mixture_;
    autoPtr<twoTemperatureModel> ttm_;
    autoPtr<femtosecondLaserModel> laser_;
    autoPtr<phaseChangeModel> phaseChange_;
    autoPtr<dropletModel> droplet_;
    autoPtr<Foam::incompressible::turbulenceModel> turbulence_;
        IOMRFZoneList MRF_;
    fv::options& fvOptions_;
    fvVectorMatrix UEqn_;
    autoPtr<advancedInterfaceCapturing> interfaceCapturing_;
    autoPtr<radiationModel> radiation_;
    autoPtr<basicThermo> thermo_;
    //autoPtr<dynamicRefineFvMesh> mesh_;

    // Physical parameters
    dimensionedScalar Tm_;          // Melting temperature
    dimensionedScalar Tv_;          // Vaporization temperature
    dimensionedScalar Lm_;          // Latent heat of melting
    dimensionedScalar Lv_;          // Latent heat of vaporization
    dimensionedScalar sigma_;       // Surface tension
    dimensionedScalar interfaceWidth_; // Interface thickness
    dimensionedScalar criticalPressure_; // Critical pressure
    dimensionedVector g_;           // Gravity
    
    // Material properties
    dimensionedScalar rho1_;        // Donor density
    dimensionedScalar rho2_;        // Ambient density
    dimensionedScalar mu1_;         // Donor viscosity
    dimensionedScalar mu2_;         // Ambient viscosity
    dimensionedScalar Cp1_;         // Donor heat capacity
    dimensionedScalar Cp2_;         // Ambient heat capacity
    dimensionedScalar k1_;          // Donor thermal conductivity
    dimensionedScalar k2_;          // Ambient thermal conductivity

    // Numerical controls
    pimpleControl& pimple_;
    scalar maxCo_;                  // Max Courant number
    scalar alphaCoNum_;             // Interface Courant number
    label pRefCell_;                // Reference pressure cell
    scalar pRefValue_;              // Reference pressure value
    scalar cumulativeContErr_;      // Cumulative continuity error

    // Private member functions
    void correctPhaseChange();
    void updateProperties();
    void calculateDropletFormation();
    dimensionedScalar calculateTotalEnergy() const;
    void checkEnergy();
    void validateFields();
        // Helper functions
    bool checkEnergyConservation
    (
        const dimensionedScalar& energyBefore,
        const dimensionedScalar& energyAfter
    ) const;

public:
    TypeName("LIFTModel");

    // Constructors
    LIFTModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Factory method
    static autoPtr<LIFTModel> New
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~LIFTModel() = default;

const pimpleControl& pimple() const { return pimple_; }
const fvVectorMatrix& UEqn() const { return UEqn_; }
const label& pRefCell() const { return pRefCell_; }
const scalar& pRefValue() const { return pRefValue_; }
//const volScalarField& phaseChangeRate() const { return phaseChangeRate_; }

    // Member Functions
    void solve();
    bool initializeFields();
    virtual void write();
    virtual bool valid() const;

    // Field access functions
    volVectorField& U() { return U_; }
    surfaceScalarField& phi() { return phi_; }
    volScalarField& p() { return p_; }
    volScalarField& alpha1() { return alpha1_; }
    volScalarField& rho() { return rho_; }
    volScalarField& Te() { return Te_; }
    volScalarField& Tl() { return Tl_; }

    // Const field access
    const volVectorField& U() const { return U_; }
    const surfaceScalarField& phi() const { return phi_; }
    const volScalarField& p() const { return p_; }
    const volScalarField& alpha1() const { return alpha1_; }
    const volScalarField& rho() const { return rho_; }
    const volScalarField& Te() const { return Te_; }
    const volScalarField& Tl() const { return Tl_; }
    const volScalarField& laserSource() const { return laserSource_; }
    const volScalarField& phaseChangeRate() const { return phaseChangeRate_; }
    const immiscibleIncompressibleTwoPhaseMixture& mixture() const 
{ 
    return mixture_; // Return reference directly instead of dereferencing pointer
}

    // Model access
    const twoTemperatureModel& ttm() const { return *ttm_; }
    const femtosecondLaserModel& laser() const { return *laser_; }
    const phaseChangeModel& phaseChange() const { return *phaseChange_; }
    const incompressible::turbulenceModel& turbulence() const { return *turbulence_; }

    // Physical parameter access  
    const dimensionedScalar& meltingTemp() const { return Tm_; }
    const dimensionedScalar& vaporizationTemp() const { return Tv_; }
    const dimensionedScalar& latentHeatMelting() const { return Lm_; }
    const dimensionedScalar& latentHeatVaporization() const { return Lv_; }
    const dimensionedScalar& surfaceTension() const { return sigma_; }
    const dimensionedScalar& interfaceWidth() const { return interfaceWidth_; }
    const dimensionedScalar& criticalPressure() const { return criticalPressure_; }

    // Material property access
    const dimensionedScalar& rho1() const { return rho1_; }
    const dimensionedScalar& rho2() const { return rho2_; }
    const dimensionedScalar& mu1() const { return mu1_; }
    const dimensionedScalar& mu2() const { return mu2_; }
    const dimensionedScalar& Cp1() const { return Cp1_; }
    const dimensionedScalar& Cp2() const { return Cp2_; }
    const dimensionedScalar& k1() const { return k1_; }
    const dimensionedScalar& k2() const { return k2_; }

    //- Solve equations
    void solvePhaseEquation();
    void solveTemperatureEquation(); 
    void solveMomentumEquation();
    void solvePressureEquation();

    // Phase state checks
    bool isLiquid(const scalar T) const
    {
        return T > Tm_.value() && T < Tv_.value();
    }

    bool isVapor(const scalar T) const
    {
        return T >= Tv_.value();
    }

    // Mesh access
    const fvMesh& mesh() const { return mesh_; }
    const dictionary& dict() const { return dict_; }
    const Time& runTime() const { return runTime_; }
    
};

} // End namespace Foam

#endif
