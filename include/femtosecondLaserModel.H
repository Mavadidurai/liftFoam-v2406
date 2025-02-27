// femtosecondLaserModel.H
#ifndef femtosecondLaserModel_H
#define femtosecondLaserModel_H

#include "fvCFD.H"
#include "DimensionValidator.H"
#include "mathematicalConstants.H"

namespace Foam
{

class femtosecondLaserModel
{
private:
    //- Reference to mesh
    const fvMesh& mesh_;
    
    //- Reference to dictionary
    const dictionary& dict_;

    // Laser parameters
    dimensionedScalar peakIntensity_;    // Peak intensity [W/m²]
    dimensionedScalar pulseWidth_;       // Pulse duration [s]
    dimensionedScalar wavelength_;       // Laser wavelength [m]
    dimensionedScalar absorptionCoeff_;  // Absorption coefficient [1/m]
    dimensionedScalar spotSize_;         // Laser spot size [m]
    dimensionedScalar pulseEnergy_;      // Pulse energy [J]
    
    // Geometrical parameters
    vector direction_;      // Laser direction
    point focus_;          // Focal point
    
    // Optional parameters
    scalar reflectivity_;   // Surface reflectivity
    bool gaussianProfile_;  // Use Gaussian profile
    label maxReflections_;  // Maximum number of reflections

    // Cached fields
    mutable tmp<volScalarField> tSource_;
    mutable bool sourceValid_;

    // Private member functions
    void calculateSource() const;
    bool validateParameters() const;
    bool checkPhysicalBounds() const;
    scalar calculateGaussianIntensity(const scalar R, const scalar z) const;
    bool isInBeam(const point& p) const;
    
    //- Energy conservation check
    bool checkEnergyConservation() const;

public:
    //- Runtime type information
    TypeName("femtosecondLaserModel");

    // Constructors
    femtosecondLaserModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~femtosecondLaserModel() = default;

    //- Delete copy constructor and assignment operator
    femtosecondLaserModel(const femtosecondLaserModel&) = delete;
    void operator=(const femtosecondLaserModel&) = delete;

    // Member Functions

    //- Update the laser model
    void update();

    //- Return the laser source term
    tmp<volScalarField> source() const;

    //- Write laser model data
    void write() const;
    
      virtual void correct()
    {
        sourceValid_ = false;
        calculateSource();
    }

    //- Return peak intensity
    const dimensionedScalar& peakIntensity() const
    {
        return peakIntensity_;
    }

    //- Return pulse width
    const dimensionedScalar& pulseWidth() const
    {
        return pulseWidth_;
    }

    //- Check if model is valid
    bool valid() const
    {
        return validateParameters() && checkPhysicalBounds();
    }
};

} // End namespace Foam

#endif

