#ifndef femtosecondLaserModel_H
#define femtosecondLaserModel_H

#include "fvCFD.H"

namespace Foam
{

class femtosecondLaserModel
{
private:
    const fvMesh& mesh_;
    const dictionary& dict_;

    scalar peakIntensity_;
    scalar pulseWidth_;
    scalar wavelength_;
    vector direction_;
    point focus_;

public:
    femtosecondLaserModel(const fvMesh& mesh, const dictionary& dict);
 bool valid() const
    {
        return peakIntensity_ > 0 && pulseWidth_ > 0 && wavelength_ > 0;
    }
    void update();
    tmp<volScalarField> source() const;
};

}

#endif
