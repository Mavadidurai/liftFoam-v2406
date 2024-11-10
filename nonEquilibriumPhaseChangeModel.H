#ifndef nonEquilibriumPhaseChangeModel_H
#define nonEquilibriumPhaseChangeModel_H

#include "phaseChangeModelBase.H"

namespace Foam
{

class nonEquilibriumPhaseChangeModel : public phaseChangeModelBase
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    dimensionedScalar Tm_;    // Melting temperature [K]
    dimensionedScalar Tv_;    // Vaporization temperature [K]
    dimensionedScalar Lm_;    // Latent heat of melting [J/kg]
    dimensionedScalar Lv_;    // Latent heat of vaporization [J/kg]
    volScalarField phaseIndicator_;

public:
    nonEquilibriumPhaseChangeModel
    (
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    );

    virtual ~nonEquilibriumPhaseChangeModel() = default;

    bool valid() const
    {
        return Tm_.value() > 0 && 
               Tv_.value() > 0 && 
               Lm_.value() > 0 && 
               Lv_.value() > 0 && 
               phaseIndicator_.size() > 0 &&
               Tm_.dimensions() == dimTemperature &&
               Tv_.dimensions() == dimTemperature &&
               Lm_.dimensions() == dimEnergy/dimMass &&
               Lv_.dimensions() == dimEnergy/dimMass;
    }
    
    virtual void correct(const volScalarField& T) override;
    virtual tmp<fvScalarMatrix> Sh(const volScalarField& h) const override;
    virtual void write() const override;
    virtual tmp<volScalarField> phaseChangeRate() const override;
};

}

#endif
