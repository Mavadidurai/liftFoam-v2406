// residualStressModel.C
#include "residualStressModel.H"

namespace Foam
{

void residualStressModel::calculate(const volScalarField& T, const volScalarField& alpha1)
{
    Info<< "Calculating residual stresses" << endl;
    
    // Update field pointers
    T_ = &T;
    alpha1_ = &alpha1;
}

tmp<volSymmTensorField> residualStressModel::residualStress() const
{
    if (!T_ || !alpha1_)
    {
        FatalErrorIn("residualStressModel::residualStress()")
            << "Temperature or phase field not set"
            << abort(FatalError);
    }

    tmp<volSymmTensorField> tStress
    (
        new volSymmTensorField
        (
            IOobject
            (
                "residualStress",
                T_->mesh().time().timeName(),
                T_->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            T_->mesh(),
            dimensionedSymmTensor
            (
                "zero",
                dimPressure,
                symmTensor::zero
            )
        )
    );

    volSymmTensorField& stress = tStress.ref();

    // Calculate thermal strain
    volScalarField thermalStrain = alpha_ * (*T_ - dimensionedScalar("Tref", dimTemperature, Tref_));
    
    // Calculate elastic constants
    scalar lambda = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0*nu_));
    scalar mu = E_ / (2.0 * (1.0 + nu_));
    
    // Calculate stress components
    forAll(T_->mesh().C(), cellI)
    {
        scalar strain = thermalStrain[cellI];
        scalar solidFraction = (*alpha1_)[cellI];
        
        // Only calculate stress in solidified regions
        if (solidFraction > 0.5)
        {
            // Calculate stress components assuming isotropic material
            scalar sigmaxx = (2.0*mu + lambda) * strain;
            
            // Create symmetric stress tensor
            stress[cellI] = symmTensor
            (
                sigmaxx, 0, 0,
                sigmaxx, 0,
                sigmaxx
            );
        }
    }

    return tStress;
}

void residualStressModel::write() const
{
    if (T_ && alpha1_)
    {
        Info<< "Residual Stress Model:" << nl
            << "  Young's modulus: " << E_ << nl
            << "  Poisson's ratio: " << nu_ << nl
            << "  Thermal expansion: " << alpha_ << nl
            << "  Reference temp: " << Tref_ << nl
            << "  Current temp range: " << min(*T_).value()
            << " - " << max(*T_).value() << nl
            << "  Solid fraction range: " << min(*alpha1_).value()
            << " - " << max(*alpha1_).value() << endl;
    }
}

}
