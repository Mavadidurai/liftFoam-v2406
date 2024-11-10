#include "dynamicRefineFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "polyTopoChange.H"
#include "cellSet.H"
#include "syncTools.H"
#include "meshTools.H"

namespace Foam {

defineTypeNameAndDebug(dynamicRefineFvMesh, 0);

dynamicRefineFvMesh::dynamicRefineFvMesh(const Time& runTime)
: 
    dynamicFvMesh(runTime),
    dynamicMeshDict_
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    )
{}

autoPtr<mapPolyMesh> dynamicRefineFvMesh::publicRefine(const labelList& cellsToRefine)
{
    return refine(cellsToRefine);
}

autoPtr<mapPolyMesh> dynamicRefineFvMesh::refine(const labelList& cellsToRefine)
{
    // Create empty topology modifier
    polyTopoChange meshMod(*this);

    // Get cells that can be refined
    labelList refinableCells(cellsToRefine);
    
    // Mark cells for refinement
    cellSet selectedCells(*this, "refineCells", refinableCells);
    
    // Refine marked cells
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    
    // Update points and topology
    if (map().hasMotionPoints())
    {
        // Get new points
        const pointField& newPoints = map().preMotionPoints();
        
        // Move mesh to new point position
        movePoints(newPoints);
    }

    // Update topology
    updateMesh(map());
    
    return map;
}

bool dynamicRefineFvMesh::update()
{
    // Read refinement controls from the dictionary
    if (!dynamicMeshDict_.found("refinementControls"))
    {
        return false;
    }

    const dictionary& controls = dynamicMeshDict_.subDict("refinementControls");
    
    // Read control parameters
    const scalar maxRefinementLevel = 
        controls.getOrDefault<scalar>("maxRefinementLevel", 2);
    const scalar refinementThreshold = 
        controls.getOrDefault<scalar>("refinementThreshold", 0.1);
    const scalar unrefineThreshold = 
        controls.getOrDefault<scalar>("unrefineThreshold", 0.05);
    
    // No mesh changes by default
    bool meshChanged = false;

    if (debug)
    {
        Info<< "Mesh refinement controls:" << nl
            << "  maxRefinementLevel: " << maxRefinementLevel << nl
            << "  refinementThreshold: " << refinementThreshold << nl
            << "  unrefineThreshold: " << unrefineThreshold << endl;
    }
    
    return meshChanged;
}

} // End namespace Foam
