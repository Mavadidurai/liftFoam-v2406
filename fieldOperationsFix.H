#ifndef fieldOperationsFix_H
#define fieldOperationsFix_H

#include "fvCFD.H"

namespace Foam
{
    //- Helper function to convert tmp<volScalarField> to scalarField
    tmp<scalarField> convertToScalarField
    (
        const tmp<volScalarField>& tmpVsf
    )
    {
        return tmp<scalarField>(new scalarField(tmpVsf().internalField()));
    }

    //- Overloaded operator for scalarField + tmp<volScalarField>
    tmp<scalarField> operator+
    (
        const scalarField& sf,
        const tmp<volScalarField>& tmpVsf
    )
    {
        return sf + convertToScalarField(tmpVsf);
    }

    //- Overloaded operator for scalarField - tmp<volScalarField>
    tmp<scalarField> operator-
    (
        const scalarField& sf,
        const tmp<volScalarField>& tmpVsf
    )
    {
        return sf - convertToScalarField(tmpVsf);
    }
}

#endif
