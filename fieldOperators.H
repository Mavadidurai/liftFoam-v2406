#ifndef fieldOperators_H
#define fieldOperators_H

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "Field.H"
#include "tmp.H"
#include "geometricOneField.H"
#include "MULES.H"

namespace Foam
{

// Field conversion utilities
template<class Type>
tmp<Field<Type>> convertToField
(
    const GeometricField<Type, fvPatchField, volMesh>& gf
)
{
    return tmp<Field<Type>>(new Field<Type>(gf.internalField()));
}

// Scalar field conversion
inline tmp<scalarField> convertToScalarField(const volScalarField& vsf)
{
    return tmp<scalarField>(new scalarField(vsf.internalField()));
}

// Phase fraction limiters
inline void limitPhaseFraction(volScalarField& alpha)
{
    alpha.primitiveFieldRef() = max(min(alpha.primitiveField(), 1.0), 0.0);
    alpha.correctBoundaryConditions();
}

// MULES-specific helpers 
inline tmp<surfaceScalarField> computePhaseFlux
(
    const surfaceScalarField& phi,
    const volScalarField& alpha,
    const dimensionedScalar& phaseChangeRate,
    const dimensionedScalar& rho
)
{
    return phi*fvc::interpolate(alpha) + 
           fvc::interpolate(phaseChangeRate/rho);
}

// Temperature field helpers
inline void limitTemperature
(
    volScalarField& T,
    const dimensionedScalar& minTemp,
    const dimensionedScalar& maxTemp
)
{
    T.primitiveFieldRef() = max(min(T.primitiveField(), maxTemp), minTemp);
    T.correctBoundaryConditions();
}

// Vector operators needed for laser and shock modeling
template<class Type>
tmp<GeometricField<Vector<Type>, fvPatchField, volMesh>> operator*
(
    const dimensionedScalar& s,
    const GeometricField<Vector<Type>, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Vector<Type>, fvPatchField, volMesh>> tRes
    (
        new GeometricField<Vector<Type>, fvPatchField, volMesh>
        (
            IOobject
            (
                "(" + s.name() + " * " + vf.name() + ")",
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            s.dimensions() * vf.dimensions()
        )
    );

    tRes.ref().primitiveFieldRef() = s.value() * vf.primitiveField();
    tRes.ref().correctBoundaryConditions();
    return tRes;
}

// Energy conservation helpers
inline dimensionedScalar calculateTotalEnergy
(
    const volScalarField& Te,
    const volScalarField& Tl,
    const dimensionedScalar& Ce,
    const dimensionedScalar& Cl,
    const volScalarField& rho
)
{
    return fvc::domainIntegrate(rho*(Ce*Te + Cl*Tl));
}

} // End namespace Foam

#endif
