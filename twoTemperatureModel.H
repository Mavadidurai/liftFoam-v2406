#ifndef twoTemperatureModel_H
#define twoTemperatureModel_H

#include "fvCFD.H"

namespace Foam
{

class twoTemperatureModel
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    volScalarField Te_;           // Electron temperature
    volScalarField Tl_;           // Lattice temperature
    dimensionedScalar Ce_;        // Electron heat capacity [J/m³/K]
    dimensionedScalar Cl_;        // Lattice heat capacity [J/m³/K]
    dimensionedScalar G_;         // Electron-phonon coupling [W/m³/K]

public:
    twoTemperatureModel(const fvMesh& mesh, const dictionary& dict);

    //- Solve the two-temperature model equations
    void solve(const volScalarField& laserSource);

    //- Return the electron thermal conductivity field
    tmp<volScalarField> electronThermalConductivity() const;

    //- Return the electron heat capacity field
    tmp<volScalarField> electronHeatCapacity() const;

    //- Return the electron-phonon coupling field
    tmp<volScalarField> electronPhononCoupling() const;

    //- Return electron heat capacity constant
    const dimensionedScalar& Ce() const { return Ce_; }

    //- Return lattice heat capacity constant
    const dimensionedScalar& Cl() const { return Cl_; }

    //- Return electron-phonon coupling constant
    const dimensionedScalar& G() const { return G_; }

    bool valid() const
    {
        return Te_.size() > 0 && 
               Tl_.size() > 0 && 
               Ce_.value() > 0 && 
               Cl_.value() > 0 && 
               G_.value() > 0 &&
               Te_.dimensions() == dimTemperature &&
               Tl_.dimensions() == dimTemperature &&
               Ce_.dimensions() == dimEnergy/dimVolume/dimTemperature &&
               Cl_.dimensions() == dimEnergy/dimVolume/dimTemperature &&
               G_.dimensions() == dimEnergy/dimVolume/dimTime/dimTemperature;
    }

    const volScalarField& Te() const { return Te_; }
    const volScalarField& Tl() const { return Tl_; }
};

}

#endif
