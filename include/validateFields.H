#ifndef validateFields_H
#define validateFields_H

Info<< "Validating fields\n" << endl;

// Check presence of required fields
if (!mesh.foundObject<surfaceScalarField>("phi"))
{
    FatalErrorIn("validateFields.H")
        << "Required field phi not found in registry"
        << abort(FatalError);
}

if (!mesh.foundObject<volScalarField>("Te"))
{
    FatalErrorIn("validateFields.H")
        << "Required field Te not found in registry"
        << abort(FatalError);
}

// Check phase fraction bounds (only if alpha1 exists)
if (mesh.foundObject<volScalarField>("alpha1"))
{
    if (min(alpha1).value() < -SMALL || max(alpha1).value() > 1 + SMALL)
    {
        FatalErrorIn("validateFields.H")
            << "Phase fraction out of bounds"
            << abort(FatalError);
    }
}
else
{
    FatalErrorIn("validateFields.H")
        << "Phase fraction field alpha1 not found in registry"
        << abort(FatalError);
}

// Check temperature positivity (only if Te and Tl exist)
if (mesh.foundObject<volScalarField>("Te") && mesh.foundObject<volScalarField>("Tl"))
{
    if (min(Te).value() < 0 || min(Tl).value() < 0)
    {
        FatalErrorIn("validateFields.H")
            << "Negative temperature detected"
            << abort(FatalError);
    }
}
else
{
    FatalErrorIn("validateFields.H")
        << "Temperature fields Te or Tl not found in registry"
        << abort(FatalError);
}

Info<< "Field validation completed\n" << endl;

#endif

