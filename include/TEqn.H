#ifndef TEqn_H  
#define TEqn_H

#include "fvCFD.H"

void TEqn(Foam::fvMesh &mesh)
{
    Foam::Info << "Solving Energy equation..." << Foam::endl;

    // Get material properties from two-temperature model
    const Foam::scalar Ce = ttm_->Ce().value();
    const Foam::scalar Cl = ttm_->Cl().value();
    const Foam::scalar G = ttm_->G().value();

    // Get thermal conductivities
    Foam::tmp<Foam::volScalarField> tke = ttm_->ke();
    Foam::tmp<Foam::volScalarField> tkl = ttm_->kl();
    const Foam::volScalarField& ke = tke();
    const Foam::volScalarField& kl = tkl();

    // Store initial energy
    const Foam::dimensionedScalar initialEnergy = 
        Foam::fvc::domainIntegrate(Ce * rho_ * Te_ + Cl * rho_ * Tl_);

    // Update laser source term
    laser_->correct();

    // Electron temperature equation
    Foam::fvScalarMatrix TeEqn
    (
        Foam::fvm::ddt(Ce * rho_, Te_)
      - Foam::fvm::laplacian(ke, Te_)
      + Foam::fvm::Sp(G * rho_, Te_)
      - G * rho_ * Tl_
     ==
        laser_->source() 
      + radiation_->Ru()
    );

    // Add phase change contribution
    if (phaseChange_)
    {
        TeEqn += phaseChange_->electronSource();
    }

    // Add droplet heat transfer
    if (droplet_ && droplet_->active())
    {
        TeEqn += droplet_->heatSource();
    }

    TeEqn.relax();
    TeEqn.solve();

    // Lattice temperature equation
    Foam::fvScalarMatrix TlEqn
    (
        Foam::fvm::ddt(Cl * rho_, Tl_)
      + Foam::fvm::div(phi_, Tl_)
      - Foam::fvm::laplacian(kl, Tl_)
      - G * rho_ * Te_
      + Foam::fvm::Sp(G * rho_, Tl_)
     ==
        Foam::dimensionedScalar("zero", Foam::dimPower / Foam::dimVolume, 0.0)
    );

    // Add phase change contribution
    if (phaseChange_)
    {
        TlEqn += phaseChange_->latticeSource();
    }

    TlEqn.relax();
    TlEqn.solve();

    // Apply physical bounds
    const Foam::dimensionedScalar Tmin("Tmin", Foam::dimTemperature, 300.0);
    const Foam::dimensionedScalar Tmax("Tmax", Foam::dimTemperature, 1e4);

    Te_ = Foam::max(Foam::min(Te_, Tmax), Tmin);
    Tl_ = Foam::max(Foam::min(Tl_, Tmax), Tmin);

    // Update material properties
    thermo_->correct();
    mixture_.correct();

    // Update phase state and radiation
    phaseChange_->correct(Tl_);
    radiation_->correct(Tl_);

    // Track energy conservation
    const Foam::dimensionedScalar finalEnergy = 
        Foam::fvc::domainIntegrate(Ce * rho_ * Te_ + Cl * rho_ * Tl_);

    const Foam::scalar energyError = Foam::mag
    (
        (finalEnergy.value() - initialEnergy.value()) /
        (Foam::mag(initialEnergy.value()) + Foam::SMALL)
    );

    if (mesh.time().writeTime())  // Fixed `mesh_` issue
    {
        Foam::Info << "Temperature solve:" << Foam::nl
            << "  Te range: " << Foam::gMin(Te_) << " - " << Foam::gMax(Te_) << " K" << Foam::nl
            << "  Tl range: " << Foam::gMin(Tl_) << " - " << Foam::gMax(Tl_) << " K" << Foam::nl
            << "  Max electron conductivity: " << Foam::max(ke).value() << Foam::nl
            << "  Max lattice conductivity: " << Foam::max(kl).value() << Foam::nl
            << "  Energy error: " << energyError * 100 << " %" << Foam::endl;

        if (phaseChange_)
        {
            Foam::Info << "Phase fractions:" << Foam::nl
                << "  Liquid: " << phaseChange_->liquidFraction() << Foam::nl
                << "  Solid: " << phaseChange_->solidFraction() << Foam::nl 
                << "  Vapor: " << phaseChange_->vaporFraction() << Foam::endl;
        }
    }
}

#endif

