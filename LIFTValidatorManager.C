// LIFTValidatorManager.C
LIFTValidatorManager::LIFTValidatorManager
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    energyTolerance_(dict.lookupOrDefault<scalar>("energyTolerance", 1e-6)),
    maxCFL_(dict.lookupOrDefault<scalar>("maxCo", 0.5)),
    maxTemperature_(dict.lookupOrDefault<scalar>("maxTemperature", 1e5)),
    minPressure_(dict.lookupOrDefault<scalar>("minPressure", 1e-6)),
    maxNonOrthogonality_(dict.lookupOrDefault<scalar>("maxNonOrthogonality", 70)),
    maxSkewness_(dict.lookupOrDefault<scalar>("maxSkewness", 4)),
    lastTotalEnergy_("lastTotalEnergy", dimEnergy, 0),
    reporter_(new ValidationReporter(mesh.time(), dict)),
    firstIteration_(true)
{
    Info<< "Initializing LIFT Validator Manager" << nl
        << "  Energy tolerance: " << energyTolerance_ << nl
        << "  Max CFL: " << maxCFL_ << nl
        << "  Max temperature: " << maxTemperature_ << nl
        << "  Min pressure: " << minPressure_ << endl;
}

bool LIFTValidatorManager::validatePreSolve
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& alpha1
)
{
    bool valid = true;
  reporter_->clearResults();

    // Temperature validation
    scalar maxTe = max(Te).value();
    bool TeValid = maxTe <= maxTemperature_;
    reporter_->addResult
    (
        "electronTemperature",
        TeValid,
        maxTe,
        maxTemperature_,
        "Electron temperature exceeds maximum",
        "Te"
    );
    valid &= TeValid;

    // Pressure validation
    scalar minP = min(p).value();
    bool pValid = minP >= minPressure_;
    reporter_->addResult
    (
        "pressure",
        pValid,
        minP,
        minPressure_,
        "Pressure below minimum",
        "p"
    );
    valid &= pValid;

    // Phase fraction validation
    scalar maxAlpha = max(alpha1).value();
    scalar minAlpha = min(alpha1).value();
    bool alphaValid = maxAlpha <= 1.0 && minAlpha >= 0.0;
    reporter_->addResult
    (
        "phaseFraction",
        alphaValid,
        maxAlpha,
        1.0,
        "Phase fraction out of bounds",
        "alpha1"
    );
    valid &= alphaValid;

   
    // Dimension validation
    valid &= DimensionValidator::checkField(Te, dimTemperature, "Te");
    valid &= DimensionValidator::checkField(Tl, dimTemperature, "Tl");
    valid &= DimensionValidator::checkField(U, dimVelocity, "U");
    valid &= DimensionValidator::checkField(p, dimPressure, "p");
    valid &= DimensionValidator::checkField(alpha1, dimless, "alpha1");

    // Physical bounds validation
    valid &= PhysicalValidator::checkTemperatureBounds(Te, 0, maxTemperature_);
    valid &= PhysicalValidator::checkTemperatureBounds(Tl, 0, maxTemperature_);
    valid &= PhysicalValidator::checkPressurePositivity(p);
    valid &= PhysicalValidator::checkPhaseFractionBounds(alpha1);

    // Mesh quality validation
    valid &= validateMeshQuality();
 // Report results
    reporter_->report();
    
    return valid;
}

bool LIFTValidatorManager::validatePostSolve
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& alpha1,
    const dimensionedScalar& totalEnergy
)
{
    bool valid = true;

    // Energy conservation
    if (!firstIteration_)
    {
        valid &= validateEnergyConservation(totalEnergy);
    }
    
    updateEnergyTracking(totalEnergy);
    firstIteration_ = false;

    // Numerical stability
    valid &= PhysicalValidator::checkNumericalStability(Te);
    valid &= PhysicalValidator::checkNumericalStability(Tl);
    valid &= PhysicalValidator::checkNumericalStability(p);

    // Boundary conditions
    valid &= validateBoundaries();

    return valid;
}

bool LIFTValidatorManager::validateTimeStep
(
    const volVectorField& U,
    const scalar deltaT
) const
{
    return PhysicalValidator::checkCFL(U, mesh_, maxCFL_);
}

bool LIFTValidatorManager::validateModels
(
    const twoTemperatureModel& ttm,
    const femtosecondLaserModel& laser,
    const phaseChangeModelBase& phaseChange,
    const advancedInterfaceCapturing& interface
) const
{
    bool valid = true;

    valid &= DimensionValidator::checkTTM(ttm);
    valid &= DimensionValidator::checkLaserModel(laser);
    valid &= DimensionValidator::checkPhaseChangeModel(phaseChange);
    valid &= DimensionValidator::checkInterfaceTracking(interface);

    return valid;
}

bool LIFTValidatorManager::validateMeshQuality() const
{
    scalar maxNonOrth = 0.0;
    scalar maxSkew = 0.0;

    forAll(mesh_.cells(), cellI)
    {
        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            maxNonOrth = max(maxNonOrth, mesh_.nonOrthogonality(cellI, faceI));
            maxSkew = max(maxSkew, mesh_.skewness(cellI, faceI));
        }
    }

    if (maxNonOrth > maxNonOrthogonality_ || maxSkew > maxSkewness_)
    {
        FatalErrorInFunction
            << "Mesh quality metrics exceeded limits:" << nl
            << "Max non-orthogonality: " << maxNonOrth << nl
            << "Max skewness: " << maxSkew
            << abort(FatalError);
        return false;
    }

    return true;
}

bool LIFTValidatorManager::validateBoundaries() const
{
    return BoundaryValidator::validateAllBoundaries
    (
        mesh_,
        mesh_.lookupObject<volScalarField>("p"),
        mesh_.lookupObject<volScalarField>("T"),
        mesh_.lookupObject<volVectorField>("U"),
        dict_
    );
}

bool LIFTValidatorManager::validateEnergyConservation
(
    const dimensionedScalar& currentEnergy
) const
{
    return PhysicalValidator::checkEnergyConservation
    (
        lastTotalEnergy_,
        currentEnergy,
        energyTolerance_
    );
}

void LIFTValidatorManager::updateEnergyTracking
(
    const dimensionedScalar& newEnergy
)
{
    lastTotalEnergy_ = newEnergy;
}

} // End namespace Foam

#endif
