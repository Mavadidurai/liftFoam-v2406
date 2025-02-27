/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of LIFT configuration management class.
    
    Handles:
    - Parameter reading and validation
    - Mesh quality assessment
    - Field verification
    - Multi-phase mixture creation
    
    The implementation ensures all parameters are:
    - Physically meaningful
    - Dimensionally correct
    - Within numerical stability limits
    - Consistent with mesh quality requirements
\*---------------------------------------------------------------------------*/

#include "LIFTConfig.H"

Foam::LIFTConfig::LIFTConfig
(
    const Time& runTime,
    const fileName& dictPath
)
:
    runTime_(runTime),
    dict_
    (
        IOobject
        (
            "liftFoamDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    )
{
    readPhysicalParameters();
    readLaserParameters();
    readNumericalLimits();
    readMeshQualityParams();

    if (!validateConfiguration())
    {
        FatalErrorIn("LIFTConfig")
            << "Invalid configuration"
            << abort(FatalError);
    }
}

void Foam::LIFTConfig::readPhysicalParameters()
{
    physicalParams_.meltTemp_ = dimensionedScalar("Tm", dimTemperature, 
        dict_.getOrDefault<scalar>("meltingTemperature", 1941.0));
    physicalParams_.vaporTemp_ = dimensionedScalar("Tv", dimTemperature,
        dict_.getOrDefault<scalar>("vaporizationTemperature", 3560.0));
    physicalParams_.donorDensity_ = dimensionedScalar("rho", dimDensity,
        dict_.getOrDefault<scalar>("donorDensity", 4500.0));
    physicalParams_.surfaceTension_ = dimensionedScalar("sigma", dimForce/dimLength,
        dict_.getOrDefault<scalar>("surfaceTension", 1.65));
    physicalParams_.latentHeatMelting_ = dimensionedScalar("Lm", dimEnergy/dimMass,
        dict_.getOrDefault<scalar>("latentHeatMelting", 3.96e5));
    physicalParams_.latentHeatVaporization_ = dimensionedScalar("Lv", dimEnergy/dimMass,
        dict_.getOrDefault<scalar>("latentHeatVaporization", 8.85e6));
    physicalParams_.interfaceWidth_ = dimensionedScalar("delta", dimLength,
        dict_.getOrDefault<scalar>("interfaceWidth", 1e-6));
    physicalParams_.criticalPressure_ = dimensionedScalar("Pc", dimPressure,
        dict_.getOrDefault<scalar>("criticalPressure", 1e8));
}

void Foam::LIFTConfig::readLaserParameters()
{
    laserParams_.pulseEnergy_ = dimensionedScalar("Ep", dimEnergy,
        dict_.get<scalar>("pulseEnergy"));
    laserParams_.pulseWidth_ = dimensionedScalar("tp", dimTime,
        dict_.get<scalar>("pulseWidth"));
    laserParams_.spotSize_ = dimensionedScalar("w0", dimLength,
        dict_.get<scalar>("spotSize"));
    laserParams_.wavelength_ = dimensionedScalar("lambda", dimLength,
        dict_.get<scalar>("wavelength"));
    laserParams_.absorptionCoeff_ = dimensionedScalar("alpha", dimless/dimLength,
        dict_.get<scalar>("absorptionCoeff"));
    laserParams_.direction_ = dict_.get<vector>("laserDirection");
    laserParams_.focus_ = dict_.get<point>("focusPoint");
}

void Foam::LIFTConfig::readNumericalLimits()
{
    dict_.readIfPresent("minTemp", limits_.minTemp);
    dict_.readIfPresent("maxTemp", limits_.maxTemp);
    dict_.readIfPresent("minPressure", limits_.minPressure);
    dict_.readIfPresent("maxPressure", limits_.maxPressure);
    dict_.readIfPresent("minDensity", limits_.minDensity);
    dict_.readIfPresent("maxCo", limits_.maxCo);
    dict_.readIfPresent("maxAlphaCo", limits_.maxAlphaCo);
    dict_.readIfPresent("energyTolerance", limits_.energyTolerance);
}

void Foam::LIFTConfig::readMeshQualityParams()
{
    dict_.readIfPresent("maxSkewness", meshQuality_.maxSkewness);
    dict_.readIfPresent("maxNonOrtho", meshQuality_.maxNonOrtho);
    dict_.readIfPresent("maxAspectRatio", meshQuality_.maxAspectRatio);
    dict_.readIfPresent("minVolume", meshQuality_.minVolume);
}

bool Foam::LIFTConfig::validateMeshQuality(const fvMesh& mesh) const
{
    // Get mesh metrics
    scalar maxSkewness = 0.0;
    scalar maxNonOrtho = 0.0;
    
    // Calculate minimum volume
    scalar minVol = min(mesh.V()).value();

    // Get mesh primitive data
    const vectorField& centres = mesh.cellCentres();
    const vectorField& areas = mesh.faceAreas();
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        // Get face normal and centers
        vector fArea = areas[faceI];
        vector fNormal = fArea/(mag(fArea) + SMALL);
        
        // Get owner/neighbor cell centers
        vector d = centres[neighbour[faceI]] - centres[owner[faceI]];
        
        // Calculate non-orthogonality angle
        scalar dMag = mag(d);
        scalar cosTheta = (d & fNormal)/(dMag + SMALL);
        scalar theta = acos(min(1.0, max(-1.0, cosTheta)));
        maxNonOrtho = max(maxNonOrtho, theta*180.0/constant::mathematical::pi);

        // Calculate skewness
        vector fCenter = mesh.faceCentres()[faceI];
        vector ownerCenter = centres[owner[faceI]];
        
        vector d1 = fCenter - ownerCenter;
        scalar skewness = mag((d1 & fNormal)/(mag(d1) + SMALL));
        maxSkewness = max(maxSkewness, skewness);
    }

    // Perform validation checks
    bool volumeOK = minVol >= meshQuality_.minVolume;
    bool skewnessOK = maxSkewness <= meshQuality_.maxSkewness;
    bool orthogonalityOK = maxNonOrtho <= meshQuality_.maxNonOrtho;

    Info<< "Mesh quality metrics:" << nl
        << "  Min volume: " << minVol << nl
        << "  Max skewness: " << maxSkewness << nl
        << "  Max non-orthogonality: " << maxNonOrtho << " deg" << endl;

    return volumeOK && skewnessOK && orthogonalityOK;
}
bool Foam::LIFTConfig::validateFields(const fvMesh& mesh) const
{
    bool valid = true;
    
    // Check required fields exist
    valid &= mesh.foundObject<volScalarField>("alpha.titanium");
    valid &= mesh.foundObject<volScalarField>("T");
    valid &= mesh.foundObject<volScalarField>("p");
    valid &= mesh.foundObject<volVectorField>("U");
    
    // Check field dimensions if they exist
    if (valid)
    {
        const volScalarField& alpha = 
            mesh.lookupObject<volScalarField>("alpha.titanium");
        const volScalarField& T = 
            mesh.lookupObject<volScalarField>("T");
        
        valid &= alpha.dimensions() == dimless;
        valid &= T.dimensions() == dimTemperature;
    }
    
    return valid;
}

bool Foam::LIFTConfig::validateDictionaries(const fvMesh& mesh) const
{
    bool valid = true;
    
    // Check required dictionaries exist
    valid &= mesh.foundObject<IOdictionary>("transportProperties");
    valid &= mesh.foundObject<IOdictionary>("thermophysicalProperties");
    
    return valid;
}

bool Foam::LIFTConfig::validateLaserParameters() const
{
    bool valid = true;
    
    valid &= laserParams_.pulseEnergy_.value() > 0;
    valid &= laserParams_.pulseWidth_.value() > 0;
    valid &= laserParams_.spotSize_.value() > 0;
    valid &= laserParams_.wavelength_.value() > 0;
    valid &= laserParams_.absorptionCoeff_.value() > 0;
    valid &= mag(laserParams_.direction_) > SMALL;
    
    return valid;
}

bool Foam::LIFTConfig::validatePhysicalParameters() const
{
    bool valid = true;
    
    valid &= physicalParams_.meltTemp_.value() > 0;
    valid &= physicalParams_.vaporTemp_.value() > physicalParams_.meltTemp_.value();
    valid &= physicalParams_.donorDensity_.value() > 0;
    valid &= physicalParams_.surfaceTension_.value() > 0;
    valid &= physicalParams_.latentHeatMelting_.value() > 0;
    valid &= physicalParams_.latentHeatVaporization_.value() > 0;
    valid &= physicalParams_.interfaceWidth_.value() > 0;
    valid &= physicalParams_.criticalPressure_.value() > 0;
    
    return valid;
}

bool Foam::LIFTConfig::validateConfiguration() const
{
    bool valid = true;

    // Validate laser parameters
    valid &= validateLaserParameters();
    
    // Validate physical parameters
    valid &= validatePhysicalParameters();
    
    // Validate numerical parameters
    valid &= limits_.minTemp > 0;
    valid &= limits_.maxTemp > limits_.minTemp;
    valid &= limits_.minPressure > 0;
    valid &= limits_.maxPressure > limits_.minPressure;
    valid &= limits_.minDensity > 0;
    valid &= limits_.maxCo > 0 && limits_.maxCo < 1;
    valid &= limits_.maxAlphaCo > 0 && limits_.maxAlphaCo < 1;
    valid &= limits_.energyTolerance > 0;
    
    return valid;
}

bool Foam::LIFTConfig::checkInitialization(const fvMesh& mesh) const
{
    bool dictionariesOK = validateDictionaries(mesh);
    bool fieldsOK = validateFields(mesh);
    bool meshOK = validateMeshQuality(mesh);
    bool physicsOK = validatePhysicalParameters();
    bool laserOK = validateLaserParameters();

    return dictionariesOK && fieldsOK && meshOK && physicsOK && laserOK;
}

autoPtr<Foam::immiscibleIncompressibleTwoPhaseMixture> 
Foam::LIFTConfig::createMixture
(
    fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi
) const
{
    return autoPtr<immiscibleIncompressibleTwoPhaseMixture>
    (
        new immiscibleIncompressibleTwoPhaseMixture(U, phi)
    );
}
