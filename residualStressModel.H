#ifndef residualStressModel_H
#define residualStressModel_H

#include "fvCFD.H"
#include "volFields.H"
#include "symmTensor.H"

namespace Foam
{

class residualStressModel
{
private:
    const volScalarField* T_;
    const volScalarField* alpha1_;
    scalar E_;           // Young's modulus
    scalar nu_;          // Poisson's ratio
    scalar alpha_;       // Thermal expansion coefficient
    scalar Tref_;        // Reference temperature

public:
    //- Default constructor
    residualStressModel():
        T_(nullptr),
        alpha1_(nullptr),
        E_(0),
        nu_(0),
        alpha_(0),
        Tref_(0)
    {}

    //- Construct from fields and dictionary
    residualStressModel
    (
        const volScalarField& T,
        const volScalarField& alpha1,
        const dictionary& dict
    ):
        T_(&T),
        alpha1_(&alpha1),
        E_(dict.get<scalar>("youngModulus")),
        nu_(dict.get<scalar>("poissonRatio")),
        alpha_(dict.get<scalar>("thermalExpansionCoeff")),
        Tref_(dict.get<scalar>("referenceTemperature"))
    {}
 bool valid() const
    {
        return T_ != nullptr && alpha1_ != nullptr && E_ > 0 && nu_ > -1 && nu_ < 0.5;
    }
    //- Calculate stresses based on current fields
    void calculate(const volScalarField& T, const volScalarField& alpha1);

    //- Return the latest calculated residual stress field
    tmp<volSymmTensorField> residualStress() const;
    
    //- Write model information
    void write() const;
};

}

#endif
