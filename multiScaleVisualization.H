#ifndef multiScaleVisualization_H
#define multiScaleVisualization_H

#include "fvCFD.H"

namespace Foam
{

class multiScaleVisualization
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;
    
    volScalarField macroScale_;
    volScalarField mesoScale_;
    volScalarField microScale_;

public:
    multiScaleVisualization(const fvMesh& mesh, const dictionary& dict);

    void update(const volScalarField& T, const volScalarField& alpha1);
    void write();
};

}

#endif
