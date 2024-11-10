// LIFTValidatorManager.H
#ifndef LIFTValidatorManager_H
#define LIFTValidatorManager_H

#include "fvCFD.H"
#include "DimensionValidator.H"
#include "PhysicalValidator.H"
#include "BoundaryValidator.H"

namespace Foam
{

class LIFTValidatorManager
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;
     std::unique_ptr<ValidationReporter> reporter_;

    // Validation thresholds
    scalar energyTolerance_;
    scalar maxCFL_;
    scalar maxTemperature_;
    scalar minPressure_;
    scalar maxNonOrthogonality_;
    scalar maxSkewness_;

    // State tracking
    dimensionedScalar lastTotalEnergy_;
    bool firstIteration_;

    // Private validation helper methods
    bool validateMeshQuality() const;
    bool validateEnergyConservation(const dimensionedScalar& currentEnergy) const;
    void updateEnergyTracking(const dimensionedScalar& newEnergy);

public:
    // Constructor
    LIFTValidatorManager
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    // Validation methods
    bool validatePreSolve
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volVectorField& U,
        const volScalarField& p,
        const volScalarField& alpha1
    );

    bool validatePostSolve
    (
        const volScalarField& Te,
        const volScalarField& Tl,
        const volVectorField& U,
        const volScalarField& p,
        const volScalarField& alpha1,
        const dimensionedScalar& totalEnergy
    );

    bool validateTimeStep
    (
        const volVectorField& U,
        const scalar deltaT
    ) const;

    bool validateModels
    (
        const twoTemperatureModel& ttm,
        const femtosecondLaserModel& laser,
        const phaseChangeModelBase& phaseChange,
        const advancedInterfaceCapturing& interface
    ) const;

    bool validateBoundaries() const;
};
