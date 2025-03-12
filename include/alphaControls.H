#ifndef ALPHACONTROLS_H
#define ALPHACONTROLS_H

// alphaControls.H

// Get alpha equation control settings from the solver dictionary
const dictionary& alphaControls = mesh.solverDict(alpha1.name());

// Interface compression coefficient
const scalar cAlpha = alphaControls.getOrDefault<scalar>("cAlpha", 1.0);

// Number of alpha corrector steps
const label nAlphaCorr = alphaControls.getOrDefault<label>("nAlphaCorr", 1);

// Number of alpha sub-cycles
const label nAlphaSubCycles = alphaControls.getOrDefault<label>("nAlphaSubCycles", 1);
#endif
