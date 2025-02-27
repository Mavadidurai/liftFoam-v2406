/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Dynamic mesh refinement class for femtosecond LIFT process simulation.
    Handles adaptive mesh refinement based on phase interface tracking and
    temperature gradients. Includes mesh quality controls and multi-level
    refinement/unrefinement capabilities.
\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineFvMesh_H
#define dynamicRefineFvMesh_H

#include "dynamicFvMesh.H"
#include "fvMeshMapper.H"
#include "mapPolyMesh.H"
#include "cellSet.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "labelList.H"
#include "autoPtr.H"
#include "syncTools.H"
#include <cmath>

namespace Foam
{

class dynamicRefineFvMesh
:   public dynamicFvMesh
{
    // Private data
    private:
        //- Mesh control dictionary
        mutable IOdictionary dynamicMeshDict_;

        //- Refinement tracking fields
        autoPtr<volScalarField> cellLevelPtr_;
        autoPtr<volScalarField> refinementCriteriaPtr_;

        //- Reference to phase field for interface tracking
        const volScalarField* alpha1Ptr_;

        //- Refinement parameters
        scalar maxRefinementLevel_;
        scalar refinementThreshold_;
        scalar unrefinementThreshold_;
        scalar maxSkewness_;
        scalar maxNonOrtho_;
        scalar minVolume_;

        //- Buffer layer control
        label nBufferLayers_;

        //- Statistics tracking
        label nRefinedCells_;
        label nUnrefinedCells_;

    // Private Member Functions
    private:
        //- Initialize refinement fields
        void initializeFields();

        //- No copy construct
        dynamicRefineFvMesh(const dynamicRefineFvMesh&) = delete;

        //- No copy assignment
        void operator=(const dynamicRefineFvMesh&) = delete;

        //- Update refinement criteria based on phase field
        void updateRefinementCriteria();

        //- Check mesh quality
        bool checkMeshQuality() const;

        //- Select cells for refinement/unrefinement
        void selectRefinementCells();

        //- Protect buffer layers around refined regions
        void protectBufferLayers(labelList& cellsToRefine) const;

        //- Calculate mesh quality metrics
        scalar maxSkewness() const;
        scalar maxNonOrthogonality() const;

        //- Refine and unrefine cells
        bool refine(const labelList& cells);
        bool unrefine(const labelList& cells);

        //- Check mesh zone validity
        bool validateCellZones() const;

    public:
        //- Runtime type information
        TypeName("dynamicRefineFvMesh");

        //- Constructor
        explicit dynamicRefineFvMesh(const IOobject& io);

        //- Destructor
        virtual ~dynamicRefineFvMesh();

    // Member Functions
        //- Update mesh refinement
        virtual bool update();

        //- Update mesh after topology change
        virtual void updateMesh(const mapPolyMesh&);

        //- Set phase field for interface tracking
        void setPhaseField(const volScalarField& alpha1)
        {
            alpha1Ptr_ = &alpha1;
        }

        //- Access functions
        const volScalarField& cellLevel() const
        {
            return cellLevelPtr_();
        }

        const volScalarField& refinementCriteria() const
        {
            return refinementCriteriaPtr_();
        }

        //- Statistics access
        label nRefinedCells() const
        {
            return nRefinedCells_;
        }

        label nUnrefinedCells() const
        {
            return nUnrefinedCells_;
        }

        //- Parameters access
        scalar maxRefinementLevel() const
        {
            return maxRefinementLevel_;
        }

        scalar refinementThreshold() const
        {
            return refinementThreshold_;
        }

        //- Mesh quality access
        scalar maxAllowedSkewness() const
        {
            return maxSkewness_;
        }

        scalar maxAllowedNonOrtho() const
        {
            return maxNonOrtho_;
        }

        //- Debug information
        bool writeDebug() const
        {
            return debug;
        }
};

} // End namespace Foam

#endif
