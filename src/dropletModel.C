#include "dropletModel.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(dropletModel, 0);

dropletModel::dropletModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& rho,
    const volVectorField& U,
    const bool isFemtosecond
)
:
    mesh_(mesh),
    rho_(rho),
    U_(U),
    L_(dict.lookupOrDefault<scalar>("characteristicLength", 1.0)),
    criticalWe_(dict.lookupOrDefault<scalar>("criticalWeberNumber", 12.0)),
    viscosity_(dict.lookupOrDefault<scalar>("viscosity", 1.0)),
    surfaceTension_(dict.lookupOrDefault<scalar>("surfaceTension", 1.0)),
    criticalTemperature_(isFemtosecond ? dict.get<scalar>("criticalTemperature") : 0),
    criticalPressure_(isFemtosecond ? dict.get<scalar>("criticalPressure") : 0),
    isFemtosecond_(isFemtosecond),
    dropletIndicator_
    (
        IOobject("dropletIndicator", mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    ejectionVelocity_
    (
        IOobject("ejectionVelocity", mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    velocity_("velocity", dimVelocity, 0.0),
    diameter_("diameter", dimLength, 0.0),
    aspectRatio_(1.0),
    circularity_(1.0)
{
    if (!valid())
    {
        FatalErrorIn("dropletModel::dropletModel")
            << "Invalid model parameters"
            << abort(FatalError);
    }
}

void dropletModel::update
(
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U
)
{
    forAll(mesh_.C(), cellI)
    {
        if (isFemtosecond_)
        {
            // Femtosecond mode
            if (T[cellI] > criticalTemperature_ && p[cellI] > criticalPressure_)
            {
                dropletIndicator_[cellI] = 1;
                scalar ejectionSpeed = sqrt(2.0 * (T[cellI] - criticalTemperature_));
                ejectionVelocity_[cellI] = -ejectionSpeed * mesh_.C()[cellI] / 
                                          (mag(mesh_.C()[cellI]) + SMALL);
            }
            else
            {
                dropletIndicator_[cellI] = 0;
                ejectionVelocity_[cellI] = vector::zero;
            }
        }
        else
        {
            // Standard mode
            if (T[cellI] > 1000 && mag(U[cellI]) > 10)
            {
                dropletIndicator_[cellI] = 1;
            }
            else
            {
                dropletIndicator_[cellI] = 0;
            }
        }
    }

    updateMetrics();
    checkBreakup();
    checkCoalescence();
}

void dropletModel::updateMetrics()
{
    calculateVelocity();
    calculateDiameter();
    calculateAspectRatio();
    calculateCircularity();
}

void dropletModel::calculateVelocity()
{
    vector centerOld = dropletCenter();
    scalar deltaT = mesh_.time().deltaTValue();

    if (deltaT > SMALL && mesh_.time().timeIndex() > mesh_.time().startTimeIndex())
    {
        static vector centerPrev = centerOld;
        velocity_ = dimensionedScalar
        (
            "velocity",
            dimLength/dimTime,
            mag(centerOld - centerPrev) / deltaT
        );
        centerPrev = centerOld;
    }
}

void dropletModel::calculateDiameter()
{
    scalar volume = dropletVolume();
    diameter_ = dimensionedScalar
    (
        "diameter",
        dimLength,
        pow(6.0 * volume / constant::mathematical::pi, 1.0/3.0)
    );
}

void dropletModel::calculateAspectRatio()
{
    scalar xMin = GREAT, xMax = -GREAT;
    scalar yMin = GREAT, yMax = -GREAT;

    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            const point& C = mesh_.C()[cellI];
            xMin = min(xMin, C.x());
            xMax = max(xMax, C.x());
            yMin = min(yMin, C.y());
            yMax = max(yMax, C.y());
        }
    }

    scalar dx = xMax - xMin;
    scalar dy = yMax - yMin;
    
    if (dy > SMALL)
    {
        aspectRatio_ = dx/dy;
    }
    else
    {
        aspectRatio_ = GREAT;
    }
}

void dropletModel::calculateCircularity()
{
    scalar perimeter = 0.0;
    scalar area = dropletVolume();

    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            const labelList& nc = mesh_.cellCells()[cellI];
            forAll(nc, nbrI)
            {
                if (dropletIndicator_[nc[nbrI]] <= 0.5)
                {
                    perimeter += mag(mesh_.Sf().boundaryField()[nbrI][cellI]);
                }
            }
        }
    }

    if (perimeter > SMALL)
    {
        circularity_ = 4.0 * constant::mathematical::pi * area / 
                      (perimeter * perimeter);
    }
    else
    {
        circularity_ = 1.0;
    }
}

bool dropletModel::isDropletFormed() const
{
    return max(dropletIndicator_).value() > 0.5;
}

void dropletModel::checkBreakup()
{
    volScalarField We = this->We()();
    
    forAll(mesh_.C(), cellI)
    {
        if (We[cellI] > criticalWe_)
        {
            dropletIndicator_[cellI] = 0.5;
        }
    }
}

void dropletModel::checkCoalescence()
{
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0 && dropletIndicator_[cellI] < 1)
        {
            const labelList& nc = mesh_.cellCells()[cellI];
            forAll(nc, nbrI)
            {
                if (dropletIndicator_[nc[nbrI]] > 0)
                {
                    dropletIndicator_[cellI] = 1;
                    dropletIndicator_[nc[nbrI]] = 1;
                }
            }
        }
    }
}

void dropletModel::calculateTrajectory()
{
    // TODO: Implement if needed for specific application
}

bool dropletModel::validateProperties() const
{
    if (L_ <= 0 || criticalWe_ <= 0 || viscosity_ <= 0 || surfaceTension_ <= 0)
    {
        return false;
    }
    
    if (isFemtosecond_ && (criticalTemperature_ <= 0 || criticalPressure_ <= 0))
    {
        return false;
    }
    
    return true;
}

bool dropletModel::validateFields() const
{
    return dropletIndicator_.size() == mesh_.nCells() &&
           ejectionVelocity_.size() == mesh_.nCells();
}

bool dropletModel::valid() const
{
    return validateProperties() && validateFields();
}

tmp<volScalarField> dropletModel::We() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "We",
                mesh_.time().timeName(),
                mesh_
            ),
            rho_ * magSqr(U_) * L_ / surfaceTension_
        )
    );
}

scalar dropletModel::dropletMass() const
{
    scalar mass = 0.0;
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            mass += rho_[cellI] * mesh_.V()[cellI];
        }
    }
    return mass;
}

scalar dropletModel::dropletKineticEnergy() const
{
    scalar ke = 0.0;
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            ke += 0.5 * rho_[cellI] * magSqr(U_[cellI]) * mesh_.V()[cellI];
        }
    }
    return ke;
}

scalar dropletModel::dropletSurfaceEnergy() const
{
    return surfaceTension_ * dropletSurfaceArea();
}

scalar dropletModel::dropletVolume() const
{
    scalar volume = 0.0;
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            volume += mesh_.V()[cellI];
        }
    }
    return volume;
}

scalar dropletModel::dropletSurfaceArea() const
{
    scalar area = 0.0;
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            const labelList& nc = mesh_.cellCells()[cellI];
            forAll(nc, nbrI)
            {
                if (dropletIndicator_[nc[nbrI]] <= 0.5)
                {
                    area += mag(mesh_.Sf().boundaryField()[nbrI][cellI]);
                }
            }
        }
    }
    return area;
}

vector dropletModel::dropletCenter() const
{
    vector center(vector::zero);
    scalar totalVolume = 0.0;
    
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            scalar cellVolume = mesh_.V()[cellI];
            center += mesh_.C()[cellI] * cellVolume;
            totalVolume += cellVolume;
        }
    }
    
    if (totalVolume > SMALL)
    {
        center /= totalVolume;
    }
    return center;
}

tmp<volSymmTensorField> dropletModel::dropletStress() const
{
    tmp<volSymmTensorField> tStress
    (
        new volSymmTensorField
        (
            IOobject
            (
                "dropletStress",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedSymmTensor
            (
                "zero",
                dimForce/dimArea,
                symmTensor::zero
            )
        )
    );
    
    volSymmTensorField& stress = tStress.ref();
    
    // Calculate surface tension contribution
    tmp<volVectorField> tGradAlpha = fvc::grad(dropletIndicator_);
    const volVectorField& gradAlpha = tGradAlpha();
    
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            vector n = gradAlpha[cellI];
            scalar magN = mag(n);
            
            if (magN > SMALL)
            {
                n /= magN;
                stress[cellI] = surfaceTension_ * magN * symm(tensor::I - (n * n));
            }
        }
    }
    
    return tStress;
}

scalar dropletModel::Oh() const
{
    scalar charLength = diameter_.value();
    if (charLength > SMALL)
    {
        return viscosity_ / 
               sqrt(rho_[0] * surfaceTension_ * charLength);
    }
    return 0.0;
}

scalar dropletModel::Re() const
{
    scalar charVelocity = velocity_.value();
    scalar charLength = diameter_.value();
    if (viscosity_ > SMALL)
    {
        return rho_[0] * charVelocity * charLength / viscosity_;
    }
    return 0.0;
}

void dropletModel::write() const
{
    Info<< "Droplet Model Statistics:" << nl
        << "  Mode: " << (isFemtosecond_ ? "Femtosecond" : "Standard") << nl
        << "  Droplet formed: " << (isDropletFormed() ? "Yes" : "No") << nl
        << "  Properties:" << nl
        << "    Mass = " << dropletMass() << " kg" << nl
        << "    Volume = " << dropletVolume() << " m³" << nl
        << "    Surface Area = " << dropletSurfaceArea() << " m²" << nl
        << "    Diameter = " << diameter_.value() << " m" << nl
        << "    Velocity = " << velocity_.value() << " m/s" << nl
        << "    Aspect Ratio = " << aspectRatio_ << nl
        << "    Circularity = " << circularity_ << nl
        << "  Dimensionless Numbers:" << nl
        << "    Weber = " << max(We()()).value() << nl
        << "    Reynolds = " << Re() << nl
        << "    Ohnesorge = " << Oh() << nl
        << "  Energetics:" << nl
        << "    Kinetic Energy = " << dropletKineticEnergy() << " J" << nl
        << "    Surface Energy = " << dropletSurfaceEnergy() << " J" << endl;

    if (mesh_.time().writeTime())
    {
        dropletIndicator_.write();
        if (isFemtosecond_)
        {
            ejectionVelocity_.write();
        }
    }
}
// Add these implementations at the end of dropletModel.C file, before the closing namespace

tmp<fvVectorMatrix> Foam::dropletModel::dragForce() const
{
    tmp<volScalarField> tCd
    (
        new volScalarField
        (
            IOobject
            (
                "Cd",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Cd", dimless, 0.0)
        )
    );

    volScalarField& Cd = tCd.ref();
    
    // Calculate drag coefficient based on Reynolds number
    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            scalar Re = rho_[cellI]*mag(U_[cellI])*diameter_.value()/viscosity_;
            if (Re < 1)
            {
                Cd[cellI] = 24.0/Re;  // Stokes regime
            }
            else if (Re < 1000)
            {
                Cd[cellI] = 24.0/Re*(1.0 + 0.15*pow(Re, 0.687));  // Schiller-Naumann
            }
            else
            {
                Cd[cellI] = 0.44;  // Newton regime
            }
        }
    }

    return
    (
        fvm::Sp(0.5*rho_*Cd*mag(U_)*dropletIndicator_/diameter_, U_)
    );
}

tmp<fvScalarMatrix> Foam::dropletModel::heatSource() const
{
    tmp<volScalarField> tNu
    (
        new volScalarField
        (
            IOobject
            (
                "Nu",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Nu", dimless, 0.0)
        )
    );

    volScalarField& Nu = tNu.ref();

    // Calculate Nusselt number based on Reynolds and Prandtl numbers
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    scalar Pr = 0.7;  // Air Prandtl number approximation

    forAll(mesh_.C(), cellI)
    {
        if (dropletIndicator_[cellI] > 0.5)
        {
            scalar Re = rho_[cellI]*mag(U_[cellI])*diameter_.value()/viscosity_;
            Nu[cellI] = 2.0 + 0.6*sqrt(Re)*pow(Pr, 1.0/3.0);  // Ranz-Marshall correlation
        }
    }

    // Return heat transfer matrix
    return
    (
        fvm::Sp
        (
            Nu*surfaceTension_/(diameter_*diameter_)*dropletIndicator_,
            T
        )
    );
}
} // End namespace Foam
