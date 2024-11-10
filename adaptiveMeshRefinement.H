#ifndef adaptiveMeshRefinement_H
#define adaptiveMeshRefinement_H

#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "volFields.H"

namespace Foam
{

class adaptiveMeshRefinement
{
private:
    //- Reference to mesh
    dynamicRefineFvMesh& mesh_;
    
    //- Reference to phase fraction field
    const volScalarField& alpha1_;
    
    // Refinement parameters
    scalar refinementThreshold_;
    scalar maxRefinementLevel_;
    scalar unrefinementThreshold_;
    
    //- Dictionary for refinement controls
    dictionary refineDict_;

    //- Private member functions
    void readControls();
    bool checkRefinementLevels(const labelList& cells) const;

public:
    // Constructors
    adaptiveMeshRefinement
    (
        dynamicRefineFvMesh& mesh,
        const volScalarField& alpha1
    );

    //- Destructor
    virtual ~adaptiveMeshRefinement() = default;

    // Member Functions
    
    //- Update the mesh refinement
    void refine();
    
    //- Check if configuration is valid
    bool valid() const;
    
    //- Access functions
    scalar refinementThreshold() const { return refinementThreshold_; }
    scalar maxRefinementLevel() const { return maxRefinementLevel_; }
    scalar unrefinementThreshold() const { return unrefinementThreshold_; }
};

} // End namespace Foam

#endif
