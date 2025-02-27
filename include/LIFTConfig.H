/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Configuration class for LIFT (Laser-Induced Forward Transfer) process.
    
    Manages and validates all configuration parameters including:
    - Physical material properties
    - Laser parameters
    - Numerical control settings
    - Mesh quality criteria
    
    Features:
    - Parameter validation
    - Mesh quality checking
    - Field validation
    - Dictionary management
    - Default value handling
    
    Used to ensure consistent and valid configuration across the LIFT solver.
    
SourceFiles
    LIFTConfig.C

Reference
    Willis & Xu (2000) Transport Phenomena and Droplet Formation
    in Pulsed Laser-Induced Forward Transfer of Liquids
\*---------------------------------------------------------------------------*/

#ifndef LIFTConfig_H
#define LIFTConfig_H

#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"


namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class LIFTConfig Declaration
\*---------------------------------------------------------------------------*/

class LIFTConfig
{
private:
    // References
    const Time& runTime_;
    IOdictionary dict_;
    
    // Physical Parameters
    struct PhysicalParameters
    {
        dimensionedScalar meltTemp_;
        dimensionedScalar vaporTemp_;
        dimensionedScalar donorDensity_;
        dimensionedScalar surfaceTension_;
        dimensionedScalar latentHeatMelting_;
        dimensionedScalar latentHeatVaporization_;
        dimensionedScalar interfaceWidth_;
        dimensionedScalar criticalPressure_;
    } physicalParams_;

    // Laser Parameters
    struct LaserParameters
    {
        dimensionedScalar pulseEnergy_;
        dimensionedScalar pulseWidth_;
        dimensionedScalar spotSize_;
        dimensionedScalar wavelength_;
        dimensionedScalar absorptionCoeff_;
        vector direction_;
        point focus_;
    } laserParams_;

    // Numerical Limits
    struct NumericalLimits
    {
        scalar minTemp = 300.0;
        scalar maxTemp = 1e5;
        scalar minPressure = 1e-6;
        scalar maxPressure = 1e9;
        scalar minDensity = 1e-6;
        scalar maxCo = 0.5;
        scalar maxAlphaCo = 0.2;
        scalar energyTolerance = 1e-6;
    } limits_;

    // Mesh Quality Parameters
    struct MeshQualityParams
    {
        scalar maxSkewness = 4.0;
        scalar maxNonOrtho = 70.0;
        scalar maxAspectRatio = 1000.0;
        scalar minVolume = 1e-15;
    } meshQuality_;

    // Private Member Functions
    void readPhysicalParameters();
    void readLaserParameters();
    void readNumericalLimits();
    void readMeshQualityParams();
    bool validateConfiguration() const;
    
    // Validation Functions
    bool validateMeshQuality(const fvMesh& mesh) const;
    bool validateFields(const fvMesh& mesh) const;
    bool validateDictionaries(const fvMesh& mesh) const;
    bool validateLaserParameters() const;
    bool validatePhysicalParameters() const;

public:
    // Constructors
    LIFTConfig
    (
        const Time& runTime,
        const fileName& dictPath
    );

    //- Destructor
    virtual ~LIFTConfig() = default;

    // Delete copy constructor and assignment
    LIFTConfig(const LIFTConfig&) = delete;
    void operator=(const LIFTConfig&) = delete;

    // Public Member Functions
    bool checkInitialization(const fvMesh& mesh) const;
    autoPtr<immiscibleIncompressibleTwoPhaseMixture> createMixture
    (
        fvMesh& mesh,
        const volVectorField& U,
        const surfaceScalarField& phi
    ) const;

    // Access Functions
    const dictionary& dict() const { return dict_; }
    const PhysicalParameters& physicalParams() const { return physicalParams_; }
    const LaserParameters& laserParams() const { return laserParams_; }
    const NumericalLimits& limits() const { return limits_; }
    const MeshQualityParams& meshQuality() const { return meshQuality_; }

    // Parameter Validation Template
    template<class Type>
    Type validateParameter
    (
        const word& name,
        const dimensionSet& dims,
        const Type& value,
        const Type& minValue,
        const Type& maxValue
    ) const;
};

} // End namespace Foam

#endif
