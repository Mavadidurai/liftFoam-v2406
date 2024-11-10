#include "adaptiveMeshRefinement.H"
#include "dynamicRefineFvMesh.H"
#include "mapPolyMesh.H"
#include "fvc.H"
#include "cellSet.H"

namespace Foam {

adaptiveMeshRefinement::adaptiveMeshRefinement
(
    dynamicRefineFvMesh& mesh,
    const volScalarField& alpha1
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    refinementThreshold_(0.1),
    maxRefinementLevel_(2),
    unrefinementThreshold_(0.05)
{
    readControls();
}

void adaptiveMeshRefinement::readControls()
{
    const dictionary& dict = mesh_.dynamicMeshDict();
    if (dict.found("refinementControls"))
    {
        const dictionary& controls = dict.subDict("refinementControls");
        
        refinementThreshold_ = controls.getOrDefault<scalar>
        (
            "refinementThreshold",
            0.1
        );
        
        maxRefinementLevel_ = controls.getOrDefault<scalar>
        (
            "maxRefinementLevel",
            2
        );
        
        unrefinementThreshold_ = controls.getOrDefault<scalar>
        (
            "unrefinementThreshold",
            0.05
        );
    }
}

bool adaptiveMeshRefinement::checkRefinementLevels(const labelList& cells) const
{
    // Get current refinement levels if stored
    const volScalarField* refinementLevelPtr =
        mesh_.findObject<volScalarField>("refinementLevel");
        
    if (!refinementLevelPtr)
    {
        // If no refinement level field exists, assume all cells are at base level
        return true;
    }
    
    const volScalarField& refinementLevel = *refinementLevelPtr;
    
    // Check each cell's refinement level
    forAll(cells, i)
    {
        if (refinementLevel[cells[i]] >= maxRefinementLevel_)
        {
            return false;
        }
    }
    
    return true;
}

void adaptiveMeshRefinement::refine()
{
    // Calculate refinement criteria
    volScalarField gradAlphaMag = mag(fvc::grad(alpha1_));
    
    // Create list of cells to refine
    DynamicList<label> cellsToRefine(mesh_.nCells()/10);
    
    forAll(mesh_.cells(), cellI)
    {
        if (gradAlphaMag[cellI] > refinementThreshold_)
        {
            cellsToRefine.append(cellI);
        }
    }
    
    // Only refine if cells marked and refinement levels okay
    if (cellsToRefine.size() > 0 && checkRefinementLevels(cellsToRefine))
    {
        // Perform refinement
        autoPtr<mapPolyMesh> map = mesh_.publicRefine(cellsToRefine);

        if (map.valid())
        {
            // Update fields using the mesh map
            mesh_.updateMesh(map());
            
            Info<< "Refined " << cellsToRefine.size() 
                << " cells." << endl;
        }
    }
}

bool adaptiveMeshRefinement::valid() const
{
    if (refinementThreshold_ <= 0 || maxRefinementLevel_ < 1)
    {
        return false;
    }
    
    if (unrefinementThreshold_ >= refinementThreshold_)
    {
        return false;
    }
    
    return true;
}

} // End namespace Foam
