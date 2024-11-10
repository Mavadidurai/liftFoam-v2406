#ifndef extremeConditionMaterialProperties_H
#define extremeConditionMaterialProperties_H

#include "fvCFD.H"

namespace Foam
{

class extremeConditionMaterialProperties
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    scalar rho0_;
    scalar Cp0_;
    scalar k0_;
    scalar Te0_;
    scalar A_;
    scalar B_;

public:
    extremeConditionMaterialProperties(const fvMesh& mesh, const dictionary& dict);
bool valid() const
    {
        return rho0_ > 0 && Cp0_ > 0 && k0_ > 0 && Te0_ > 0;
    }
    void update(const volScalarField& T);
    tmp<volScalarField> rho(const volScalarField& T, const volScalarField& p) const;
    tmp<volScalarField> Cp(const volScalarField& T, const volScalarField& p) const;
    tmp<volScalarField> k(const volScalarField& T, const volScalarField& p) const;
    tmp<volScalarField> electronThermalConductivity(const volScalarField& Te) const;
    tmp<volScalarField> electronHeatCapacity(const volScalarField& Te) const;
};

}

#endif
