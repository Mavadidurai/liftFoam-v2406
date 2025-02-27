/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Model for handling radiation effects in LIFT process simulation.
    
    Features:
    - Temperature-dependent absorptivity and emissivity
    - Radiative heat transfer calculation
    - Stefan-Boltzmann radiation
    - Property bounds enforcement
    - Automatic boundary condition handling
    
    The model calculates:
    - Surface absorption of laser energy
    - Radiative heat losses
    - Temperature-dependent material properties
    
SourceFiles
    radiationModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModel_H
#define radiationModel_H

#include "fvCFD.H"
#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class radiationModel Declaration
\*---------------------------------------------------------------------------*/

class radiationModel
{
protected:
    // Protected Data
    
        //- Reference to mesh
        const fvMesh& mesh_;
        
        //- Reference to dictionary
        dictionary dict_;
        
        //- Material absorptivity field [-]
        mutable volScalarField absorptivity_;
        
        //- Material emissivity field [-]
        mutable volScalarField emissivity_;
        
        // Temperature Dependencies
        
            //- Temperature coefficient [1/K]
            scalar tempCoeff_;
            
            //- Reference temperature [K]
            scalar referenceTemp_;

    // Protected Member Functions

        //- Calculate temperature-dependent properties
        void calculateProperties(const volScalarField& T) const;
        
        //- Calculate radiative heat transfer
        tmp<volScalarField> calculateRadiativeHeatTransfer
        (
            const volScalarField& T
        ) const;

public:
    // Constructors

        //- Construct from mesh and dictionary
        radiationModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );
        
        //- Destructor
        virtual ~radiationModel() {}

    // Member Functions

        //- Update radiation properties
        virtual void correct(const volScalarField& T);
        
        //- Return radiative source term
        virtual tmp<volScalarField> Ru() const;
        
        //- Update temperature-dependent properties
        virtual void updateProperties() const;
        
        //- Write radiation model data
        virtual void write() const;

    // Access Functions

        //- Return absorptivity field
        const volScalarField& absorptivity() const { return absorptivity_; }
        
        //- Return emissivity field
        const volScalarField& emissivity() const { return emissivity_; }
};

} // End namespace Foam

#endif
