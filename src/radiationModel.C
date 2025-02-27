/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of radiation model for LIFT process simulation.
    
    Handles:
    - Temperature-dependent property updates
    - Radiative heat transfer calculations
    - Property bounds enforcement
    - Interface tracking
    
    Uses Stefan-Boltzmann law for radiation calculations and includes
    temperature-dependent material properties.
\*---------------------------------------------------------------------------*/

#include "radiationModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        radiationModel Implementation
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

radiationModel::radiationModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    absorptivity_
    (
        IOobject
        (
            "absorptivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "absorptivity",
            dimless,
            dict.lookupOrDefault<scalar>("absorptivity", 0.5)
        )
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "emissivity",
            dimless,
            dict.lookupOrDefault<scalar>("emissivity", 0.5)
        )
    ),
    tempCoeff_(dict.lookupOrDefault<scalar>("temperatureCoefficient", 0.001)),
    referenceTemp_(dict.lookupOrDefault<scalar>("referenceTemperature", 300.0))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void radiationModel::correct(const volScalarField& T)
{
    calculateProperties(T);
}

void radiationModel::calculateProperties(const volScalarField& T) const
{
    // Update temperature-dependent properties
    forAll(mesh_.C(), cellI)
    {
        scalar tempDiff = T[cellI] - referenceTemp_;
        scalar tempFactor = 1.0 + tempCoeff_ * tempDiff;
        
        // Apply temperature correction
        absorptivity_[cellI] *= tempFactor;
        emissivity_[cellI] *= tempFactor;
        
        // Enforce physical bounds
        absorptivity_[cellI] = min(max(absorptivity_[cellI], 0.0), 1.0);
        emissivity_[cellI] = min(max(emissivity_[cellI], 0.0), 1.0);
    }

    // Update boundary conditions
    absorptivity_.correctBoundaryConditions();
    emissivity_.correctBoundaryConditions();
}


tmp<volScalarField> radiationModel::calculateRadiativeHeatTransfer
(
    const volScalarField& T
) const
{
    const scalar sigmaSB = 5.67e-8; // Stefan-Boltzmann constant

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "qrad",
                mesh_.time().timeName(),
                mesh_
            ),
            emissivity_ * sigmaSB * pow4(T)
        )
    );
}

tmp<volScalarField> radiationModel::Ru() const
{
    if (!mesh_.foundObject<volScalarField>("T"))
    {
        FatalErrorIn("radiationModel::Ru()")
            << "Temperature field 'T' not found in the mesh."
            << exit(FatalError);
    }

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    tmp<volScalarField> tRu = calculateRadiativeHeatTransfer(T);
    
    return tRu;
}

void radiationModel::updateProperties() const
{
    if (mesh_.foundObject<volScalarField>("T"))
    {
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        calculateProperties(T);
    }
}

void radiationModel::write() const
{
    Info<< "Radiation Model:" << nl
        << "  Average absorptivity: " << gAverage(absorptivity_) << nl
        << "  Average emissivity: " << gAverage(emissivity_) << nl
        << "  Temperature coefficient: " << tempCoeff_ << nl
        << "  Reference temperature: " << referenceTemp_ << " K" << endl;
}

} // End namespace Foam
