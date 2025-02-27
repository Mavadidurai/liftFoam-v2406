// extremeConditionMaterialProperties.H
#ifndef extremeConditionMaterialProperties_H
#define extremeConditionMaterialProperties_H

#include "fvCFD.H"
#include "volFields.H"
#include "DimensionValidator.H"

namespace Foam
{

class extremeConditionMaterialProperties
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    // Essential material properties 
    dimensionedScalar rho0_;     // Reference density
    dimensionedScalar Cp0_;      // Specific heat
    dimensionedScalar k0_;       // Thermal conductivity
    dimensionedScalar Tm_;       // Melting temperature 
    dimensionedScalar Tv_;       // Vaporization temperature
    dimensionedScalar Lm_;       // Latent heat melting
    dimensionedScalar Lv_;       // Latent heat vaporization
    dimensionedScalar sigma_;    // Surface tension

    // Temperature dependence
    scalar alphaT_;   // Thermal expansion
    scalar betaT_;    // Conductivity coefficient  
    scalar gammaT_;   // Heat capacity coefficient
    
    // Electron coupling
    scalar G0_;       // Coupling coefficient
    scalar Ce0_;      // Electron heat capacity 
    scalar Cl0_;      // Lattice heat capacity

public:
    // Constructor
    extremeConditionMaterialProperties
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    bool valid() const;

    // Property accessors
    const dimensionedScalar& rho0() const { return rho0_; }
    const dimensionedScalar& Cp0() const { return Cp0_; }
    const dimensionedScalar& k0() const { return k0_; }
    const dimensionedScalar& Tm() const { return Tm_; }
    const dimensionedScalar& Tv() const { return Tv_; }

    // Temperature dependent properties
    void update(const volScalarField& T);
    tmp<volScalarField> rho(const volScalarField& T) const;
    tmp<volScalarField> Cp(const volScalarField& T) const;
    tmp<volScalarField> k(const volScalarField& T) const;

    // Electron properties 
    tmp<volScalarField> electronThermalConductivity(const volScalarField& Te) const;
    tmp<volScalarField> electronHeatCapacity(const volScalarField& Te) const;
    tmp<volScalarField> electronPhononCoupling
    (
        const volScalarField& Te, 
        const volScalarField& Tl
    ) const;

    // Phase properties
    tmp<volScalarField> surfaceTension(const volScalarField& T) const;
    tmp<volScalarField> phaseChangeRate(const volScalarField& T) const;
};

}

#endif
