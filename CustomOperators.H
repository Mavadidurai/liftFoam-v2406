#ifndef CustomOperators_H
#define CustomOperators_H

#include "fvCFD.H"

namespace Foam
{

// Operator overload for subtracting a vector from a volVectorField
template<class Cmpt>
tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> operator-
(
    const GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& vf,
    const Vector<Cmpt>& v
)
{
    tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> tRes
    (
        new GeometricField<Vector<Cmpt>, fvPatchField, volMesh>
        (
            IOobject
            (
                "(" + vf.name() + " - " + Foam::name(v) + ")",
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            vf.dimensions()
        )
    );

    GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& res = tRes.ref();

    forAll(res, cellI)
    {
        res[cellI] = vf[cellI] - v;
    }

    return tRes;
}

// Operator overload for adding a vector to a volVectorField
template<class Cmpt>
tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> operator+
(
    const GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& vf,
    const Vector<Cmpt>& v
)
{
    tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> tRes
    (
        new GeometricField<Vector<Cmpt>, fvPatchField, volMesh>
        (
            IOobject
            (
                "(" + vf.name() + " + " + Foam::name(v) + ")",
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            vf.dimensions()
        )
    );

    GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& res = tRes.ref();

    forAll(res, cellI)
    {
        res[cellI] = vf[cellI] + v;
    }

    return tRes;
}

// Operator overload for dividing a volVectorField by a scalar
template<class Cmpt>
tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> operator/
(
    const GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& vf,
    const scalar& s
)
{
    tmp<GeometricField<Vector<Cmpt>, fvPatchField, volMesh>> tRes
    (
        new GeometricField<Vector<Cmpt>, fvPatchField, volMesh>
        (
            IOobject
            (
                "(" + vf.name() + " / " + Foam::name(s) + ")",
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            vf.dimensions()/dimless
        )
    );

    GeometricField<Vector<Cmpt>, fvPatchField, volMesh>& res = tRes.ref();

    forAll(res, cellI)
    {
        res[cellI] = vf[cellI] / s;
    }

    return tRes;
}

// Add more custom operators as needed for your LIFT process simulation

}

#endif
