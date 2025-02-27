#ifndef advancedInterfaceCapturing_H
#define advancedInterfaceCapturing_H

#include "fvCFD.H"
#include "surfaceFields.H"
#include "MULES.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

namespace Foam
{

class advancedInterfaceCapturing
{
private:
    // References
    const fvMesh& mesh_;
    volScalarField& alpha1_;         // Phase fraction
    const surfaceScalarField& phi_;  // Flux field
    const immiscibleIncompressibleTwoPhaseMixture& mixture_;
    const volScalarField& T_;        // Temperature field for phase change

    // LIFT-specific parameters
    scalar compressionCoeff_;        // Interface compression
    scalar interfaceThickness_;      // For phase transition region
    scalar baseContactAngle_;        // Base contact angle
    scalar meltingTemp_;            // Melting temperature
    scalar vaporTemp_;              // Vaporization temperature
    scalar latentHeat_;             // Latent heat of fusion
    
    // Temperature-dependent properties
    scalar surfaceTensionCoeff_;     // Base surface tension coefficient
    scalar surfaceTensionGradient_; // Temperature gradient of surface tension
    scalar maxSurfaceTension_;      // Maximum surface tension value
    scalar minSurfaceTension_;      // Minimum surface tension value

    // Interface tracking fields
    volVectorField interfaceNormal_;  // Interface normal vectors
    volScalarField curvature_;        // Interface curvature
    volScalarField phaseChangeRate_;  // Phase change rate
    volScalarField ejectionVelocity_; // Droplet ejection velocity

    // Private helper functions
    void calculateInterfaceNormal();
    void calculateCurvature();
    void calculatePhaseChange();
    void calculateEjectionVelocity();
    void applyContactAngle();
    void updateTemperatureDependentProps();
    tmp<surfaceScalarField> computeCompressionFlux() const;
    scalar calculateLocalSurfaceTension(const scalar T) const;

public:
    // Constructor
    advancedInterfaceCapturing
    (
        const fvMesh& mesh,
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T
    );

    //- Delete default copy constructor and assignment operator
    advancedInterfaceCapturing(const advancedInterfaceCapturing&) = delete;
    void operator=(const advancedInterfaceCapturing&) = delete;

    // Member functions
    void correct();  // Main solver routine

    // Access functions
    const volVectorField& interfaceNormal() const { return interfaceNormal_; }
    const volScalarField& curvature() const { return curvature_; }
    const volScalarField& phaseChangeRate() const { return phaseChangeRate_; }
    const volScalarField& ejectionVelocity() const { return ejectionVelocity_; }
    scalar surfaceTensionCoeff() const { return surfaceTensionCoeff_; }

    // Write functions
    void write() const;
};

} // End namespace Foam

#endif
