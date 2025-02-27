// femtosecondLaserModel.C
#include "femtosecondLaserModel.H"

namespace Foam
{

defineTypeNameAndDebug(femtosecondLaserModel, 0);

femtosecondLaserModel::femtosecondLaserModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    peakIntensity_
    (
        "peakIntensity",
        dimPower/dimArea,
        dict.get<scalar>("peakIntensity")
    ),
    pulseWidth_
    (
        "pulseWidth",
        dimTime,
        dict.get<scalar>("pulseWidth")
    ),
    wavelength_
    (
        "wavelength",
        dimLength,
        dict.get<scalar>("wavelength")
    ),
    absorptionCoeff_
    (
        "absorptionCoeff",
        dimless/dimLength,
        dict.getOrDefault<scalar>("absorptionCoeff", 1e8)
    ),
    spotSize_
    (
        "spotSize",
        dimLength,
        dict.get<scalar>("spotSize")
    ),
    pulseEnergy_
    (
        "pulseEnergy",
        dimEnergy,
        dict.get<scalar>("pulseEnergy")
    ),
    direction_(dict.get<vector>("direction")),
    focus_(dict.get<point>("focus")),
    reflectivity_(dict.getOrDefault<scalar>("reflectivity", 0.0)),
    gaussianProfile_(dict.getOrDefault<bool>("gaussianProfile", true)),
    maxReflections_(dict.getOrDefault<label>("maxReflections", 1)),
    sourceValid_(false)
{
    // Normalize direction vector
    direction_ /= mag(direction_) + SMALL;

    if (!valid())
    {
        FatalErrorInFunction
            << "Invalid laser parameters" << nl
            << "  Peak intensity: " << peakIntensity_.value() << nl
            << "  Pulse width: " << pulseWidth_.value() << nl
            << "  Wavelength: " << wavelength_.value() << nl
            << "  Spot size: " << spotSize_.value()
            << abort(FatalError);
    }

    Info<< "Femtosecond laser model initialized:" << nl
        << "  Peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Pulse width: " << pulseWidth_.value() << " s" << nl
        << "  Wavelength: " << wavelength_.value() << " m" << nl
        << "  Spot size: " << spotSize_.value() << " m" << nl
        << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
        << "  Focus: " << focus_ << nl
        << "  Direction: " << direction_ << endl;
}

void femtosecondLaserModel::update()
{
    // Invalidate cached source field
    sourceValid_ = false;
}

bool femtosecondLaserModel::validateParameters() const
{
    bool valid = true;

    // Check dimensions
    valid &= peakIntensity_.dimensions() == dimPower/dimArea;
    valid &= pulseWidth_.dimensions() == dimTime;
    valid &= wavelength_.dimensions() == dimLength;
    valid &= absorptionCoeff_.dimensions() == dimless/dimLength;
    valid &= spotSize_.dimensions() == dimLength;
    valid &= pulseEnergy_.dimensions() == dimEnergy;

    // Check values
    if (peakIntensity_.value() <= 0)
    {
        WarningInFunction
            << "Non-positive peak intensity: " << peakIntensity_.value()
            << endl;
        valid = false;
    }

    if (pulseWidth_.value() <= 0)
    {
        WarningInFunction
            << "Non-positive pulse width: " << pulseWidth_.value()
            << endl;
        valid = false;
    }

    if (wavelength_.value() <= 0)
    {
        WarningInFunction
            << "Non-positive wavelength: " << wavelength_.value()
            << endl;
        valid = false;
    }

    return valid;
}

bool femtosecondLaserModel::checkPhysicalBounds() const
{
    bool valid = true;

    // Check femtosecond regime (100fs - 1000fs typical)
    if (pulseWidth_.value() < 1e-15 || pulseWidth_.value() > 1e-12)
    {
        WarningInFunction
            << "Pulse width outside femtosecond range: "
            << pulseWidth_.value() << " s"
            << endl;
        valid = false;
    }

    // Check wavelength (typically around 800nm for Ti:Sapphire)
    if (wavelength_.value() < 1e-7 || wavelength_.value() > 1e-6)
    {
        WarningInFunction
            << "Wavelength outside typical range: "
            << wavelength_.value() << " m"
            << endl;
        valid = false;
    }

    // Check intensity (typical range 10¹² - 10¹⁶ W/m²)
    if (peakIntensity_.value() > 1e16)
    {
        WarningInFunction
            << "Very high peak intensity may cause material damage: "
            << peakIntensity_.value() << " W/m²"
            << endl;
    }

    return valid;
}

scalar femtosecondLaserModel::calculateGaussianIntensity
(
    const scalar R,
    const scalar z
) const
{
    if (gaussianProfile_)
    {
        return exp(-2.0*sqr(R)/(wavelength_.value()*z/constant::mathematical::pi));
    }
    else
    {
        return R <= spotSize_.value() ? 1.0 : 0.0;
    }
}

bool femtosecondLaserModel::isInBeam(const point& p) const
{
    vector r = p - focus_;
    scalar z = (r & direction_);
    r -= z*direction_;
    
    return mag(r) <= 2.0*spotSize_.value() && z >= 0;
}

bool femtosecondLaserModel::checkEnergyConservation() const
{
    if (!tSource_.valid())
    {
        return true;
    }

    const dimensionedScalar totalEnergy = 
        fvc::domainIntegrate(tSource_()*mesh_.time().deltaT());

    return totalEnergy.value() <= 1.1*pulseEnergy_.value();
}

void femtosecondLaserModel::calculateSource() const
{
    if (sourceValid_)
    {
        return;
    }

    tSource_ = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "laserSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& source = tSource_.ref();
    const scalar t = mesh_.time().value();

    forAll(mesh_.C(), cellI)
    {
        const point& cellCenter = mesh_.C()[cellI];
        
        if (isInBeam(cellCenter))
        {
            vector r = cellCenter - focus_;
            scalar z = (r & direction_);
            r -= z*direction_;

            scalar R = mag(r);
            scalar temporalTerm = exp(-4.0*log(2.0)*sqr(t/pulseWidth_.value()));
            scalar spatialTerm = calculateGaussianIntensity(R, z);
            scalar absorptionTerm = exp(-absorptionCoeff_.value()*z);

            source[cellI] = peakIntensity_.value() * 
                           temporalTerm * 
                           spatialTerm * 
                           absorptionTerm * 
                           (1.0 - reflectivity_);
        }
    }

    sourceValid_ = true;

    if (!checkEnergyConservation())
    {
        WarningInFunction
            << "Energy conservation violated in laser source calculation"
            << endl;
    }
}

tmp<volScalarField> femtosecondLaserModel::source() const
{
    calculateSource();
    return tSource_;
}

void femtosecondLaserModel::write() const
{
    Info<< "Femtosecond laser model:" << nl
        << "  Peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Pulse width: " << pulseWidth_.value() << " s" << nl
        << "  Wavelength: " << wavelength_.value() << " m" << nl
        << "  Spot size: " << spotSize_.value() << " m" << nl
        << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
        << "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
        << "  Reflectivity: " << reflectivity_ << nl
        << "  Focus: " << focus_ << nl
        << "  Direction: " << direction_ << endl;

    if (tSource_.valid())
    {
        const dimensionedScalar maxIntensity = max(tSource_());
        const dimensionedScalar totalEnergy = 
            fvc::domainIntegrate(tSource_()*mesh_.time().deltaT());

        Info<< "Source statistics:" << nl
            << "  Maximum intensity: " << maxIntensity.value() << " W/m³" << nl
            << "  Total energy: " << totalEnergy.value() << " J" << endl;
    }
}

} // End namespace Foam
