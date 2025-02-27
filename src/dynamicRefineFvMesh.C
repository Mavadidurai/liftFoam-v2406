/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Dynamic mesh refinement implementation for femtosecond LIFT process.
\*---------------------------------------------------------------------------*/

#include "dynamicRefineFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "meshTools.H"
#include "fvCFD.H"
#include "autoPtr.H"

namespace Foam 
{
    
    defineTypeNameAndDebug(dynamicRefineFvMesh, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Modified constructor in dynamicRefineFvMesh.C
dynamicRefineFvMesh::dynamicRefineFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshDict_
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    alpha1Ptr_(nullptr),
    maxRefinementLevel_(dynamicMeshDict_.lookupOrDefault<scalar>("maxRefinementLevel", 2)),
    refinementThreshold_(dynamicMeshDict_.lookupOrDefault<scalar>("refinementThreshold", 0.1)),
    unrefinementThreshold_(dynamicMeshDict_.lookupOrDefault<scalar>("unrefinementThreshold", 0.05)),
    maxSkewness_(dynamicMeshDict_.lookupOrDefault<scalar>("maxSkewness", 4.0)),
    maxNonOrtho_(dynamicMeshDict_.lookupOrDefault<scalar>("maxNonOrtho", 70.0)), 
    minVolume_(dynamicMeshDict_.lookupOrDefault<scalar>("minVolume", 1e-13)),
    nBufferLayers_(dynamicMeshDict_.lookupOrDefault<label>("nBufferLayers", 1)),
    nRefinedCells_(0),
    nUnrefinedCells_(0)
{
    // Initialize fields after base class construction
    initializeFields();
}

// Add new helper function 
void dynamicRefineFvMesh::initializeFields()
{
    // Create cellLevel field
    cellLevelPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0)
        )
    );

    // Create refinementCriteria field  
    refinementCriteriaPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "refinementCriteria",
                time().timeName(), 
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0)
        )
    );

    // Store in registry
    if (cellLevelPtr_.valid())
    {
        objectRegistry::store(cellLevelPtr_.ptr());
    }
    if (refinementCriteriaPtr_.valid()) 
    {
        objectRegistry::store(refinementCriteriaPtr_.ptr());
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void dynamicRefineFvMesh::updateRefinementCriteria()
{
    if (!alpha1Ptr_)
    {
        return;
    }

    const volScalarField& alpha1 = *alpha1Ptr_;
    volScalarField& criteria = refinementCriteriaPtr_();

    // Interface gradient-based refinement
    tmp<volVectorField> tgradAlpha = fvc::grad(alpha1);
    const volVectorField& gradAlpha = tgradAlpha();
    const volScalarField magGradAlpha(mag(gradAlpha));

    // Calculate interface width
    scalar interfaceWidth = 0.0;
    label nInterfaceCells = 0;
    forAll(alpha1, cellI)
    {
        if (alpha1[cellI] > 0.01 && alpha1[cellI] < 0.99)
        {
            interfaceWidth += 1.0/magGradAlpha[cellI];
            nInterfaceCells++;
        }
    }
    reduce(interfaceWidth, sumOp<scalar>());
    reduce(nInterfaceCells, sumOp<label>());
    if (nInterfaceCells > 0)
    {
        interfaceWidth /= nInterfaceCells;
    }

    // Set refinement criteria
    forAll(criteria, cellI)
    {
        scalar alpha = alpha1[cellI];
        
        if (alpha > 0.01 && alpha < 0.99)
        {
            // Interface region
            criteria[cellI] = magGradAlpha[cellI]*interfaceWidth;
        }
        else
        {
            // Bulk regions
            criteria[cellI] = 0.0;
        }
    }
}

void dynamicRefineFvMesh::updateMesh(const mapPolyMesh& map)
{
    // Update base mesh first
    dynamicFvMesh::updateMesh(map);

    // Map fields
    if (cellLevelPtr_.valid())
    {
        cellLevelPtr_().primitiveFieldRef().map
        (
            cellLevelPtr_().primitiveField(),
            map.cellMap()
        );
        cellLevelPtr_().correctBoundaryConditions();
    }

    if (refinementCriteriaPtr_.valid())
    {
        refinementCriteriaPtr_().primitiveFieldRef().map
        (
            refinementCriteriaPtr_().primitiveField(),
            map.cellMap()
        );
        refinementCriteriaPtr_().correctBoundaryConditions();
    }
}

bool dynamicRefineFvMesh::checkMeshQuality() const
{
    if (debug)
    {
        Info<< "Checking mesh quality..." << endl;
    }

    const scalar maxSkew = maxSkewness();
    const scalar maxNonOrth = maxNonOrthogonality();
    const scalar minVol = gMin(V());

    bool meshOK = true;

    if (maxSkew > maxSkewness_)
    {
        meshOK = false;
        Info<< "High mesh skewness: " << maxSkew << endl;
    }

    if (maxNonOrth > maxNonOrtho_)
    {
        meshOK = false;
        Info<< "High non-orthogonality: " << maxNonOrth << endl;
    }

    if (minVol < minVolume_)
    {
        meshOK = false;
        Info<< "Small cell volume: " << minVol << endl;
    }

    return meshOK;
}
scalar dynamicRefineFvMesh::maxSkewness() const
{
    scalar maxSkew = 0.0;

    const pointField& points = this->points();
    const faceList& faces = this->faces();
    const cellList& cells = this->cells();
    const vectorField& cellCentres = this->cellCentres();
    const vectorField& faceCentres = this->faceCentres();

    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        forAll(c, faceI)
        {
            label faceID = c[faceI];
            
            if (this->isInternalFace(faceID))
            {
                vector d = faceCentres[faceID] - cellCentres[cellI];
                // Use areaNormal() instead of deprecated normal()
                vector s = faces[faceID].areaNormal(points);
                s /= mag(s) + SMALL;

                // Calculate face skewness
                scalar skewness = mag((d & s)/(mag(d) + SMALL));
                maxSkew = max(maxSkew, skewness);
            }
        }
    }

    reduce(maxSkew, maxOp<scalar>());
    return maxSkew;
}

scalar dynamicRefineFvMesh::maxNonOrthogonality() const
{
    scalar maxNonOrtho = 0.0;

    const pointField& points = this->points();
    const faceList& faces = this->faces();
    const vectorField& cellCentres = this->cellCentres();
    const labelList& owner = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();

    // Process faces up to nInternalFaces
    const label nIntFaces = this->nInternalFaces();

    for(label faceI = 0; faceI < nIntFaces; faceI++)
    {
        // Get owner and neighbour cell centers
        vector d = cellCentres[neighbour[faceI]] - cellCentres[owner[faceI]];
        
        // Use areaNormal() instead of deprecated normal()
        vector s = faces[faceI].areaNormal(points);
        s /= mag(s) + SMALL;

        // Calculate non-orthogonality angle
        scalar cosPhi = (d & s)/(mag(d) + SMALL);
        cosPhi = min(1.0, max(-1.0, cosPhi));
        scalar phi = acos(cosPhi);

        maxNonOrtho = max(maxNonOrtho, phi);
    }

    maxNonOrtho *= 180.0/constant::mathematical::pi;
    reduce(maxNonOrtho, maxOp<scalar>());
    
    return maxNonOrtho;
}
void dynamicRefineFvMesh::selectRefinementCells()
{
    labelList cellsToRefine;
    labelList cellsToUnrefine;
    
    volScalarField& criteria = refinementCriteriaPtr_();
    volScalarField& cellLevel = cellLevelPtr_();

    // Select cells for refinement/unrefinement
    forAll(criteria, cellI)
    {
        if (criteria[cellI] > refinementThreshold_ && cellLevel[cellI] < maxRefinementLevel_)
        {
            cellsToRefine.append(cellI);
        }
        else if (criteria[cellI] < unrefinementThreshold_ && cellLevel[cellI] > 0)
        {
            cellsToUnrefine.append(cellI);
        }
    }

    // Protect buffer layers
    protectBufferLayers(cellsToRefine);

    // Perform refinement/unrefinement
    if (cellsToRefine.size() > 0)
    {
        refine(cellsToRefine);
        nRefinedCells_ += cellsToRefine.size();
    }

    if (cellsToUnrefine.size() > 0)
    {
        unrefine(cellsToUnrefine);
        nUnrefinedCells_ += cellsToUnrefine.size();
    }
}

void dynamicRefineFvMesh::protectBufferLayers(labelList& cellsToRefine) const
{
    for (label n = 0; n < nBufferLayers_; n++)
    {
        labelList currentLayer = cellsToRefine;

        forAll(currentLayer, i)
        {
            const labelList& nbrCells = cellCells()[currentLayer[i]];
            
            forAll(nbrCells, nbrI)
            {
                if (!cellsToRefine.found(nbrCells[nbrI]))
                {
                    cellsToRefine.append(nbrCells[nbrI]);
                }
            }
        }
    }
}

bool dynamicRefineFvMesh::update()
{
    // Only refine if phase field is set
    if (!alpha1Ptr_)
    {
        return false;
    }

    // Update refinement criteria
    updateRefinementCriteria();

    // Check mesh quality
    if (!checkMeshQuality())
    {
        return false;
    }

    // Select and refine cells
    selectRefinementCells();

    if (debug)
    {
        Info<< "Refined cells: " << nRefinedCells_
            << " Unrefined cells: " << nUnrefinedCells_ << endl;
    }

    return true;
}

bool dynamicRefineFvMesh::refine(const labelList& cells)
{
    // Implementation for cell refinement
    // This is a placeholder - actual refinement implementation needed
    return true;
}

bool dynamicRefineFvMesh::unrefine(const labelList& cells)
{
    // Implementation for cell unrefinement
    // This is a placeholder - actual unrefinement implementation needed
    return true;
}

dynamicRefineFvMesh::~dynamicRefineFvMesh()
{
    if (cellLevelPtr_.valid())
    {
        objectRegistry::erase("cellLevel");
    }
    if (refinementCriteriaPtr_.valid())
    {
        objectRegistry::erase("refinementCriteria");
    }
}
bool dynamicRefineFvMesh::validateCellZones() const
{
    // Check cell zone existence and orientation
    const cellZoneMesh& czm = this->cellZones();
    
    forAll(czm, zoneI)
    {
        const cellZone& cz = czm[zoneI];
        
        // Verify zone cells have valid owners
        const labelList& cells = cz;
        forAll(cells, i)
        {
            if (cells[i] < 0 || cells[i] >= nCells())
            {
                return false;
            }
        }
    }
    return true;
}
} // End namespace Foam
