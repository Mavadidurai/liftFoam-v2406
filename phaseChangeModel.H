#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "phaseChangeModelBase.H"
#include "fvCFD.H"

namespace Foam
{

class phaseChangeModel
:
    public phaseChangeModelBase
{
    // Private data
private:
    const fvMesh& mesh_;
    const volScalarField& T_;
    volScalarField& alpha1_;
    const dictionary& dict_;          
    
    // Physical parameters
    dimensionedScalar meltingTemperature_;
    dimensionedScalar latentHeat_;
    dimensionedScalar undercoolingCoeff_;
    dimensionedScalar gradualMeltingRate_;
    dimensionedScalar gradualSolidificationRate_;
    dimensionedScalar interfaceWidth_;
    dimensionedScalar dTmdp_;
    dimensionedScalar pRef_;
    
    // Fields
    volScalarField phaseIndicator_;
    volScalarField rho_;
    mutable tmp<volScalarField> phaseChangeRate_;

    // Private Member Functions
    void calculatePhaseChangeRate() const;

public:
    //- Runtime type information
    TypeName("phaseChangeModel");

    // Constructors
    phaseChangeModel
    (
        const fvMesh& mesh,
        const volScalarField& T,
        volScalarField& alpha1,
        const dictionary& dict
    );

    //- Destructor
    virtual ~phaseChangeModel() = default;

    // Member Functions
    virtual void correct(const volScalarField& T);
    virtual tmp<fvScalarMatrix> Sh(const volScalarField& h) const;
    virtual void write() const;
    virtual tmp<volScalarField> phaseChangeRate() const;
};

} // End namespace Foam

#endif
