// dropletMetrics.C
#include "dropletMetrics.H"

namespace Foam
{

dropletMetrics::dropletMetrics(const fvMesh& mesh, const volScalarField& alpha1)
:
    mesh_(mesh),
    alpha1_(alpha1)
{}

scalar dropletMetrics::aspectRatio() const
{
    scalar xMin = GREAT, xMax = -GREAT, yMin = GREAT, yMax = -GREAT;

    forAll(mesh_.C(), cellI)
    {
        if (alpha1_[cellI] > 0.5)
        {
            const point& C = mesh_.C()[cellI];
            xMin = min(xMin, C.x());
            xMax = max(xMax, C.x());
            yMin = min(yMin, C.y());
            yMax = max(yMax, C.y());
        }
    }

    return (xMax - xMin) / (yMax - yMin);
}

scalar dropletMetrics::circularity() const
{
    scalar perimeter = 0;
    scalar area = 0;

    forAll(mesh_.C(), cellI)
    {
        if (alpha1_[cellI] > 0.5)
        {
            area += mesh_.V()[cellI];
            
            labelList neighborCells = mesh_.cellCells()[cellI];
            forAll(neighborCells, neighborI)
            {
                if (alpha1_[neighborCells[neighborI]] <= 0.5)
                {
                    perimeter += mag(mesh_.Sf().boundaryField()[neighborI][cellI]);
                }
            }
        }
    }

    return 4 * M_PI * area / (perimeter * perimeter);
}

}
