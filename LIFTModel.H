#ifndef LIFTModel_H
#define LIFTModel_H

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "adaptiveMeshRefinement.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "twoTemperatureModel.H"
#include "femtosecondLaserModel.H"
#include "ultraFastShockWaveModel.H"
#include "dropletModel.H"
#include "femtosecondDropletModel.H"
#include "extremeConditionMaterialProperties.H"
#include "nonEquilibriumPhaseChangeModel.H"
#include "timeResolvedAnalytics.H"
#include "multiScaleVisualization.H"
#include "adaptiveTimeStep.H"
#include "vacuumBoundaryCondition.H"
#include "twoPhaseMixtureProperties.H"
#include "residualStressModel.H"
#include "radiationModel.H"
#include "advancedInterfaceCapturing.H"
#include "curvatureModel.H"
#include "dropletMetrics.H"
#include "DimensionValidator.H"
#include "ValidationReporter.H"
#include <memory>

namespace Foam
{

class LIFTModel
{
private:
    // Base objects and references
    const fvMesh& mesh_;
    const dictionary& dict_;
    immiscibleIncompressibleTwoPhaseMixture& mixture_;
    std::unique_ptr<ValidationReporter> reporter_;

    // Primary fields
    volScalarField rho_;
    volVectorField U_;
    volScalarField p_;
    volScalarField alpha1_;
    surfaceScalarField phi_;

    // Temperature fields
    volScalarField Te_;
    volScalarField Tl_;

    // Process fields
    volScalarField laserSource_;
    volScalarField phaseChangeRate_;
    volScalarField phaseIndicator_;
    volScalarField interfaceEnergy_;
    
    // Core models
    std::unique_ptr<twoTemperatureModel> ttm_;
    std::unique_ptr<femtosecondLaserModel> laser_;
    std::unique_ptr<ultraFastShockWaveModel> shockWave_;
    std::unique_ptr<dropletModel> droplet_;
    std::unique_ptr<femtosecondDropletModel> fsDroplet_;
    std::unique_ptr<extremeConditionMaterialProperties> material_;
    std::unique_ptr<nonEquilibriumPhaseChangeModel> phaseChange_;
    std::unique_ptr<timeResolvedAnalytics> analytics_;
    std::unique_ptr<multiScaleVisualization> visualization_;
    std::unique_ptr<adaptiveTimeStep> timeStep_;
    std::unique_ptr<residualStressModel> stressModel_;
    std::unique_ptr<twoPhaseMixtureProperties> mixtureProps_;
    std::unique_ptr<advancedInterfaceCapturing> interfaceCapturing_;
    std::unique_ptr<curvatureModel> curvature_;
    std::unique_ptr<dropletMetrics> dropletMetrics_;
    std::unique_ptr<radiationModel> radiation_;
    std::unique_ptr<adaptiveMeshRefinement> meshRefinement_;
    std::unique_ptr<vacuumBoundaryCondition> vacuumBC_;

    // Physical parameters
    dimensionedScalar Tm_;        // Melting temperature
    dimensionedScalar Tv_;        // Vaporization temperature
    dimensionedScalar Lm_;        // Latent heat of melting
    dimensionedScalar Lv_;        // Latent heat of vaporization
    dimensionedScalar sigma_;     // Surface tension
    dimensionedScalar interfaceWidth_;  // Interface width

    // Energy tracking
    dimensionedScalar initialEnergy_;
    dimensionedScalar lastEnergy_;
    scalar maxEnergyChange_;

    // Solution controls
    bool solveEnergy_;
    bool solveMomentum_;
    bool adaptiveMesh_;
    bool writeFields_;

    // Private validation methods
    bool validateMixture() const;
    bool validateFields() const;
    bool validateModels() const;
    bool validateEnergy() const;
    bool validateBoundaries() const;
    bool validateMeshQuality() const;
    
    // Private initialization methods
    void initializeFields();
    void initializeLaserSource();
    void initializeMaterialRegions();
    void updatePhaseIndicator();

    // Private helper methods
    void calculateInterfaceEnergy();
    dimensionedScalar calculateTotalEnergy() const;
    bool checkEnergyConservation
    (
        const dimensionedScalar& energyBefore,
        const dimensionedScalar& energyAfter
    ) const;
    bool checkTimestepStability(const scalar deltaT) const;
    void updateEnergyTracking();
    void writeSummary() const;

public:
    // Constructors
    LIFTModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    ~LIFTModel();

    // Deleted copy constructor and assignment
    LIFTModel(const LIFTModel&) = delete;
    void operator=(const LIFTModel&) = delete;

    // Member Functions
        //- Solve the LIFT model equations
        void solve();

        //- Write the results
        void write();

        //- Calculate adaptive time step
        scalar calculateDeltaT() const;

        //- Update external fields
        void updateFields
        (
            volScalarField& alpha1,
            volVectorField& U,
            volScalarField& p,
            volScalarField& p_rgh,
            volScalarField& rho
        );

        //- Apply boundary conditions
        void applyBoundaryConditions();

        //- Validate field dimensions
        void validateDimensions() const;

        //- Check if model is valid
        bool valid() const;

        // Error reporting
        void reportError(const std::string& msg) const;
        void reportWarning(const std::string& msg) const;

    // Access functions
        //- Get electron temperature field
        const volScalarField& Te() const { return Te_; }
        
        //- Get lattice temperature field
        const volScalarField& Tl() const { return Tl_; }
        
        //- Get pressure field
        const volScalarField& p() const { return p_; }
        
        //- Get velocity field
        const volVectorField& U() const { return U_; }
        
        //- Get phase fraction field
        const volScalarField& alpha1() const { return alpha1_; }
        
        //- Get density field
        const volScalarField& rho() const { return rho_; }
        
        //- Get laser source field
        const volScalarField& laserSource() const { return laserSource_; }
        
        //- Get phase change rate field
        const volScalarField& phaseChangeRate() const { return phaseChangeRate_; }
        
        //- Get phase indicator field
        const volScalarField& phaseIndicator() const { return phaseIndicator_; }
        
        //- Get interface energy field
        const volScalarField& interfaceEnergy() const { return interfaceEnergy_; }
        
        //- Get flux field
        const surfaceScalarField& phi() const { return phi_; }
        
        //- Get mixture properties
        const immiscibleIncompressibleTwoPhaseMixture& mixture() const { return mixture_; }

    // Model access
        //- Get two-temperature model
        const twoTemperatureModel& ttm() const { return *ttm_; }
        
        //- Get laser model
        const femtosecondLaserModel& laser() const { return *laser_; }
        
        //- Get phase change model
        const nonEquilibriumPhaseChangeModel& phaseChange() const { return *phaseChange_; }
};

} // End namespace Foam

#endif
