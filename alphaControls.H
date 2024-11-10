#ifndef alphaControls_H
#define alphaControls_H

// Get alpha controls
const dictionary& alphaControls = mesh.solverDict(alpha1.name());

// Read core parameters
const label nAlphaSubCycles(alphaControls.getOrDefault<label>("nAlphaSubCycles", 2));
const label nAlphaCorr(alphaControls.getOrDefault<label>("nAlphaCorr", 2));
const scalar cAlpha(alphaControls.getOrDefault<scalar>("cAlpha", 1.0));

// Alpha equation schemes
const word alphaScheme("div(phi,alpha)");
const word alpharScheme("div(phirb,alpha)");

#endif
