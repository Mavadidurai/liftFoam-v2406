#ifndef dynamicRefineFvMesh_H
#define dynamicRefineFvMesh_H

#include "dynamicFvMesh.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "cellSet.H"
#include "dictionary.H"
#include "IOdictionary.H"

namespace Foam {

class dynamicRefineFvMesh : public dynamicFvMesh
{
private:
    //- Reference to refinement dictionary
    mutable IOdictionary dynamicMeshDict_;

public:
    TypeName("dynamicRefineFvMesh");

    // Constructors
    dynamicRefineFvMesh(const Time& runTime);
    
    //- Destructor
    virtual ~dynamicRefineFvMesh() {}

    // Public Member Functions
    
    //- Return const access to refinement controls dictionary
    const dictionary& dynamicMeshDict() const
    {
        return dynamicMeshDict_;
    }
    
    //- Wrapper to access refine method
    autoPtr<mapPolyMesh> publicRefine(const labelList& cellsToRefine);

protected:
    //- Refine cells according to cellsToRefine list
    virtual autoPtr<mapPolyMesh> refine(const labelList& cellsToRefine);

    //- Update the mesh topology
    virtual bool update();

private:
    //- Disallow default bitwise copy construct
    dynamicRefineFvMesh(const dynamicRefineFvMesh&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const dynamicRefineFvMesh&) = delete;
};

} // End namespace Foam

#endif
