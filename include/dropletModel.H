#ifndef dropletModel_H 
#define dropletModel_H

#include "fvCFD.H"
#include "volFields.H"

namespace Foam
{

class dropletModel
{
protected:
    // References
    const fvMesh& mesh_;
    const volScalarField& rho_;
    const volVectorField& U_;

    // Physical properties
    scalar L_;                // Characteristic length
    scalar criticalWe_;       // Critical Weber number
    scalar viscosity_;        // Dynamic viscosity
    scalar surfaceTension_;   // Surface tension
    
    // Femtosecond properties
    scalar criticalTemperature_; // Critical temperature for phase transition
    scalar criticalPressure_;    // Critical pressure for phase transition
    bool isFemtosecond_;         // Mode indicator

    // Fields
    volScalarField dropletIndicator_;  // Droplet phase indicator field
    volVectorField ejectionVelocity_;  // Ejection velocity field (femtosecond)

    // Metrics
    dimensionedScalar velocity_;    // Droplet velocity
    dimensionedScalar diameter_;    // Equivalent spherical diameter
    scalar aspectRatio_;            // Droplet aspect ratio
    scalar circularity_;            // Droplet circularity

    // Private member functions for metrics calculation
    void calculateVelocity();
    void calculateDiameter();
    void calculateAspectRatio();
    void calculateCircularity();
    void updateMetrics();

    // Additional private helper functions
    scalar dropletVolume() const;
    scalar dropletSurfaceArea() const;
    vector dropletCenter() const;
    tmp<volSymmTensorField> dropletStress() const;
    bool validateProperties() const;
    bool validateFields() const;

public:
    //- Runtime type information
    TypeName("dropletModel");

    // Constructors
    dropletModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField& rho,
        const volVectorField& U,
        const bool isFemtosecond = false
    );

    //- Destructor
    virtual ~dropletModel() = default;

    // Delete copy constructor and assignment operator
    dropletModel(const dropletModel&) = delete;
    void operator=(const dropletModel&) = delete;

    // Member Functions
    //- Update droplet state
    virtual void update
    (
        const volScalarField& T,
        const volScalarField& p,
        const volVectorField& U
    );

    //- Check droplet formation
    bool isDropletFormed() const;

    //- Check breakup conditions
    void checkBreakup();

    //- Check coalescence conditions
    void checkCoalescence();

    //- Calculate droplet trajectory
    void calculateTrajectory();

    //- Model validation
    bool valid() const;

    //- Write function
    void write() const;

    // Dimensionless numbers
    scalar Oh() const;  // Ohnesorge number 
    scalar Re() const;  // Reynolds number

    // Access functions
    const dimensionedScalar& velocity() const { return velocity_; }
    const dimensionedScalar& diameter() const { return diameter_; }
    scalar aspectRatio() const { return aspectRatio_; }
    scalar circularity() const { return circularity_; }
    tmp<volVectorField> ejectionVelocity() const { return ejectionVelocity_; }
    const volScalarField& dropletIndicator() const { return dropletIndicator_; }
    bool isFemtosecondMode() const { return isFemtosecond_; }
    bool active() const { return isDropletFormed(); }
    tmp<fvVectorMatrix> dragForce() const;
    tmp<fvScalarMatrix> heatSource() const;

    // Helper functions
    tmp<volScalarField> We() const;
    scalar dropletMass() const;
    scalar dropletKineticEnergy() const;
    scalar dropletSurfaceEnergy() const;
};

} // End namespace Foam

#endif
