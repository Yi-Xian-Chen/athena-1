//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gi_disk_3d_eos_v2.cpp
//! \brief 3D rotating GI disk initial conditions with weak sources for Helmholtz EOS.

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

namespace {

constexpr Real c_cgs = 2.99792458e10;
constexpr Real msun_cgs = 1.98847e33;
constexpr Real sqrt_two_pi = 2.5066282746310005024;

Real gm_cgs, rg_cgs;
Real r_in, r_out, r_ref;
Real sigma_ref, sigma_slope;
Real poverrho_ref, poverrho_slope;
Real q_ref;
Real ye_init, ye_vacuum, t_vacuum;
Real rho_floor;
Real vr_init;
int ye_index;
bool ye_source_enabled;
Real ye_eq, ye_tau;
bool enable_qw_source, pressure_corrected_rotation;
Real source_radius, source_radius_inv2;
Real l_nu, l_nubar, l_mu, eps_nu, eps_nubar, eps_mu;
Real qdot_scale, yedot_scale, qw_rho_min, qw_energy_fraction_limit, ye_floor, ye_ceil;
int inner_bc_mode, outer_bc_mode;

Real CylRadius(Real r, Real theta);
Real SigmaAtR(Real r_cyl);
Real PoverRhoAtR(Real r_cyl);
Real OmegaK(Real r_cyl);
Real ScaleHeight(Real r_cyl);
Real RhoMid(Real r_cyl);
Real DiskDensity(Real r, Real theta);
bool IsAtmosphere(Real r, Real theta);
Real VphiSquared(Real r, Real theta);
void FillCellState(MeshBlock *pmb, int k, int j, int i);
Real FermiApprox(int n, Real eta);
void LuminosityData(AthenaArray<Real> &out);
Real QWEta(Real rho, Real temp, Real ye);
void FreeNucleonFractions(EquationOfState *peos, Real rho, Real temp, Real ye,
                          Real *ypfree, Real *ynfree);
Real QdotQW(Real temp_mev, Real rho, Real x, Real ypfree, Real ynfree,
            Real etaqw, AthenaArray<Real> &le);
Real YeSource(Real temp_mev, Real ye, Real ypfree, Real ynfree,
              Real x, Real rho, AthenaArray<Real> &le);
void DiskSourceTerms(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &prim_scalar,
                     const AthenaArray<Real> &bcc,
                     AthenaArray<Real> &cons,
                     AthenaArray<Real> &cons_scalar);
void DiskInnerX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                      FaceField &b, Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                      FaceField &b, Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh);

} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gi_disk_3d_eos_v2" << std::endl
        << "This problem generator requires --coord=spherical_polar.";
    ATHENA_ERROR(msg);
  }
  if (!GENERAL_EOS) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gi_disk_3d_eos_v2" << std::endl
        << "This problem generator requires --eos=general/helmholtz.";
    ATHENA_ERROR(msg);
  }
  if (NSCALARS <= 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gi_disk_3d_eos_v2" << std::endl
        << "This problem generator requires at least one passive scalar for Ye.";
    ATHENA_ERROR(msg);
  }

  const Real m_msun = pin->GetOrAddReal("problem", "M_central_msun", 3.0);
  const Real gm_default = 6.67430e-8*m_msun*msun_cgs;
  gm_cgs = pin->GetOrAddReal("problem", "GM", gm_default);
  rg_cgs = gm_cgs/SQR(c_cgs);

  const Real r_in_rg = pin->GetOrAddReal("problem", "r_in_rg", 50.0);
  const Real r_out_rg = pin->GetOrAddReal("problem", "r_out_rg", 150.0);
  const Real r_ref_rg = pin->GetOrAddReal("problem", "r_ref_rg", 100.0);
  r_in = r_in_rg*rg_cgs;
  r_out = r_out_rg*rg_cgs;
  r_ref = r_ref_rg*rg_cgs;

  q_ref = pin->GetOrAddReal("problem", "Q_ref", 1.0);
  const Real cs_frac = pin->GetOrAddReal("problem", "cs_frac", 0.03);
  poverrho_ref = pin->GetOrAddReal("problem", "pressure_over_rho_ref",
                                   SQR(cs_frac*c_cgs));
  poverrho_slope = pin->GetOrAddReal("problem", "pressure_over_rho_slope", 0.0);
  sigma_slope = pin->GetOrAddReal("problem", "Sigma_slope", -1.5);
  sigma_ref = pin->GetOrAddReal("problem", "Sigma_ref",
                                std::sqrt(poverrho_ref)*OmegaK(r_ref)/(PI*6.67430e-8*q_ref));

  ye_index = pin->GetOrAddInteger("hydro", "helm_ye_index", 0);
  ye_init = pin->GetOrAddReal("problem", "Ye0", 0.4);
  ye_vacuum = pin->GetOrAddReal("problem", "Ye_vacuum", ye_init);
  rho_floor = pin->GetOrAddReal("problem", "rho_floor", 1.0e6);
  t_vacuum = pin->GetOrAddReal("problem", "T_vacuum", -1.0);
  vr_init = pin->GetOrAddReal("problem", "vr", 0.0);

  ye_source_enabled = pin->GetOrAddBoolean("problem", "ye_source_enabled", false);
  ye_eq = pin->GetOrAddReal("problem", "ye_eq", ye_init);
  ye_tau = pin->GetOrAddReal("problem", "ye_tau", 1.0e30);
  pressure_corrected_rotation =
      pin->GetOrAddBoolean("problem", "pressure_corrected_rotation", false);
  enable_qw_source = pin->GetOrAddBoolean("problem", "enable_qw_source", false);
  source_radius = pin->GetOrAddReal("problem", "source_radius", pin->GetReal("mesh", "x1min"));
  source_radius_inv2 = 1.0/SQR(source_radius);
  qdot_scale = pin->GetOrAddReal("problem", "qdot_scale", 1.0);
  yedot_scale = pin->GetOrAddReal("problem", "yedot_scale", 1.0);
  qw_rho_min = pin->GetOrAddReal("problem", "qw_rho_min", 10.0*rho_floor);
  qw_energy_fraction_limit =
      pin->GetOrAddReal("problem", "qw_energy_fraction_limit", 0.1);
  ye_floor = pin->GetOrAddReal("problem", "ye_floor", 1.0e-6);
  ye_ceil = pin->GetOrAddReal("problem", "ye_ceil", 1.0 - 1.0e-6);

  l_nu = pin->GetOrAddReal("problem", "L_nu", 0.0);
  l_nubar = pin->GetOrAddReal("problem", "L_nubar", 0.0);
  l_mu = pin->GetOrAddReal("problem", "L_mu", 0.0);
  eps_nu = pin->GetOrAddReal("problem", "eps_nu", 10.0);
  eps_nubar = pin->GetOrAddReal("problem", "eps_nubar", 10.0);
  eps_mu = pin->GetOrAddReal("problem", "eps_mu", 10.0);

  inner_bc_mode = pin->GetOrAddInteger("problem", "inner_bc_mode", 2);
  outer_bc_mode = pin->GetOrAddInteger("problem", "outer_bc_mode", 2);

  if (r_in <= 0.0 || r_out <= r_in || r_ref <= 0.0 || q_ref <= 0.0
      || poverrho_ref <= 0.0 || sigma_ref <= 0.0 || rho_floor <= 0.0
      || ye_init <= 0.0 || ye_init >= 1.0 || ye_vacuum <= 0.0 || ye_vacuum >= 1.0
      || ye_index < 0 || ye_index >= NSCALARS
      || ye_tau <= 0.0 || source_radius <= 0.0 || eps_nu <= 0.0 || eps_nubar <= 0.0
      || eps_mu <= 0.0 || qdot_scale < 0.0 || yedot_scale < 0.0
      || qw_rho_min < 0.0 || qw_energy_fraction_limit < 0.0
      || ye_floor < 0.0 || ye_ceil > 1.0 || ye_floor >= ye_ceil
      || t_vacuum == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gi_disk_3d_eos_v2" << std::endl
        << "Invalid disk parameters.";
    ATHENA_ERROR(msg);
  }

  if (ye_source_enabled || enable_qw_source) EnrollUserExplicitSourceFunction(DiskSourceTerms);
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    if (inner_bc_mode == 2) EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1Diode);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    if (outer_bc_mode == 2) EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1Diode);
  }

  if (Globals::my_rank == 0) {
    std::cout << "gi_disk_3d_eos_v2: GM=" << gm_cgs << ", rg=" << rg_cgs
              << " cm, Q_ref=" << q_ref << ", Sigma_ref=" << sigma_ref
              << ", P/rho_ref=" << poverrho_ref << ", Ye0=" << ye_init
              << ", Ye_vacuum=" << ye_vacuum
              << ", QW=" << enable_qw_source
              << std::endl;
  }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(8);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "xalpha");
  SetUserOutputVariableName(2, "rho_mid");
  SetUserOutputVariableName(3, "h_over_r");
  SetUserOutputVariableName(4, "qdot");
  SetUserOutputVariableName(5, "yedot");
  SetUserOutputVariableName(6, "ye");
  SetUserOutputVariableName(7, "atm");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const int il = is - NGHOST;
  const int iu = ie + NGHOST;
  const int jl = js - NGHOST;
  const int ju = je + NGHOST;
  const int kl = (block_size.nx3 > 1) ? ks - NGHOST : ks;
  const int ku = (block_size.nx3 > 1) ? ke + NGHOST : ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        FillCellState(this, k, j, i);
      }
    }
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real r_cell[NSCALARS];
        for (int n=0; n<NSCALARS; ++n) r_cell[n] = pscalars->r(n,k,j,i);

        const Real rho = phydro->w(IDN,k,j,i);
        const Real pres = phydro->w(IPR,k,j,i);
        const Real temp = peos->TFromRhoP(rho, pres, r_cell);
        const Real ye = r_cell[ye_index];
        const Real r_cyl = std::max(CylRadius(pcoord->x1v(i), pcoord->x2v(j)), r_in);
        user_out_var(0,k,j,i) = temp;
        user_out_var(1,k,j,i) = peos->XalphaFromRhoTYe(rho, temp, ye);
        user_out_var(2,k,j,i) = RhoMid(r_cyl);
        user_out_var(3,k,j,i) = ScaleHeight(r_cyl)/r_cyl;
        user_out_var(4,k,j,i) = 0.0;
        user_out_var(5,k,j,i) = 0.0;
        user_out_var(6,k,j,i) = ye;
        user_out_var(7,k,j,i) = IsAtmosphere(pcoord->x1v(i), pcoord->x2v(j)) ? 1.0 : 0.0;
        if (enable_qw_source) {
          AthenaArray<Real> le;
          le.NewAthenaArray(12);
          LuminosityData(le);
          Real ypfree, ynfree;
          FreeNucleonFractions(peos, rho, temp, ye, &ypfree, &ynfree);
          const Real x = std::sqrt(std::max(0.0, 1.0 - SQR(source_radius/pcoord->x1v(i))));
          const Real temp_mev = temp*8.617333262e-11;
          const Real eta = QWEta(rho, temp, ye);
          Real qdot_erg = qdot_scale*QdotQW(temp_mev, rho, x, ypfree, ynfree, eta, le)
              *1.602176634e-6;
          if (rho < qw_rho_min) qdot_erg = std::max(0.0, qdot_erg);
          user_out_var(4,k,j,i) = qdot_erg;
          user_out_var(5,k,j,i) = (rho < qw_rho_min && qdot_erg <= 0.0) ? 0.0
              : yedot_scale*YeSource(temp_mev, ye, ypfree, ynfree, x, rho, le);
        }
      }
    }
  }
}

namespace {

Real CylRadius(Real r, Real theta) {
  return r*std::sin(theta);
}

Real SigmaAtR(Real r_cyl) {
  const Real rc = std::min(std::max(r_cyl, r_in), r_out);
  return sigma_ref*std::pow(rc/r_ref, sigma_slope);
}

Real PoverRhoAtR(Real r_cyl) {
  const Real rc = std::min(std::max(r_cyl, r_in), r_out);
  return poverrho_ref*std::pow(rc/r_ref, poverrho_slope);
}

Real OmegaK(Real r_cyl) {
  return std::sqrt(gm_cgs/(r_cyl*r_cyl*r_cyl));
}

Real ScaleHeight(Real r_cyl) {
  return std::sqrt(PoverRhoAtR(r_cyl))/OmegaK(r_cyl);
}

Real RhoMid(Real r_cyl) {
  return SigmaAtR(r_cyl)/(sqrt_two_pi*ScaleHeight(r_cyl));
}

Real DiskDensity(Real r, Real theta) {
  if (r < r_in || r > r_out) return rho_floor;
  const Real r_cyl = CylRadius(r, theta);

  const Real z = r*std::cos(theta);
  const Real h = ScaleHeight(r_cyl);
  return std::max(rho_floor, RhoMid(r_cyl)*std::exp(-0.5*SQR(z/h)));
}

bool IsAtmosphere(Real r, Real theta) {
  if (r < r_in || r > r_out) return true;
  const Real r_cyl = CylRadius(r, theta);

  const Real z = r*std::cos(theta);
  const Real h = ScaleHeight(r_cyl);
  const Real rho_raw = RhoMid(r_cyl)*std::exp(-0.5*SQR(z/h));
  return rho_raw <= rho_floor;
}

Real VphiSquared(Real r, Real theta) {
  const Real r_cyl = CylRadius(r, theta);
  if (r_cyl <= 0.0) return 0.0;
  if (!pressure_corrected_rotation || r_cyl < r_in || r_cyl > r_out) return gm_cgs/r_cyl;

  const Real h = ScaleHeight(r_cyl);
  const Real z = r*std::cos(theta);
  const Real dln_h = 0.5*poverrho_slope + 1.5;
  const Real dln_rho = sigma_slope - dln_h + SQR(z/h)*dln_h;
  const Real pressure_term = PoverRhoAtR(r_cyl)*(dln_rho + poverrho_slope);
  return gm_cgs/r_cyl + pressure_term;
}

void FillCellState(MeshBlock *pmb, int k, int j, int i) {
  const Real r = pmb->pcoord->x1v(i);
  const Real theta = pmb->pcoord->x2v(j);
  const Real r_cyl = std::max(CylRadius(r, theta), r_in);
  const bool is_atmosphere = IsAtmosphere(r, theta);
  const Real rho = DiskDensity(r, theta);
  const Real poverrho = PoverRhoAtR(r_cyl);
  const Real vphi2 = std::max(VphiSquared(r, theta), 0.0);
  const Real vphi = std::sqrt(vphi2);

  Real r_cell[NSCALARS];
  for (int n=0; n<NSCALARS; ++n) r_cell[n] = 0.0;
  r_cell[ye_index] = is_atmosphere ? ye_vacuum : ye_init;
  const Real pres = (is_atmosphere && t_vacuum > 0.0)
      ? pmb->peos->PresFromRhoT(rho, t_vacuum, r_cell) : rho*poverrho;
  const Real egas = (is_atmosphere && t_vacuum > 0.0)
      ? pmb->peos->EgasFromRhoT(rho, t_vacuum, r_cell)
      : pmb->peos->EgasFromRhoP(rho, pres, r_cell);

  pmb->phydro->u(IDN,k,j,i) = rho;
  pmb->phydro->u(IM1,k,j,i) = rho*vr_init;
  pmb->phydro->u(IM2,k,j,i) = 0.0;
  pmb->phydro->u(IM3,k,j,i) = rho*vphi;
  pmb->phydro->u(IEN,k,j,i) = egas + 0.5*rho*(SQR(vr_init) + vphi2);

  pmb->phydro->w(IDN,k,j,i) = rho;
  pmb->phydro->w(IVX,k,j,i) = vr_init;
  pmb->phydro->w(IVY,k,j,i) = 0.0;
  pmb->phydro->w(IVZ,k,j,i) = vphi;
  pmb->phydro->w(IPR,k,j,i) = pres;

  for (int n=0; n<NSCALARS; ++n) {
    pmb->pscalars->s(n,k,j,i) = rho*r_cell[n];
    pmb->pscalars->r(n,k,j,i) = r_cell[n];
  }
}

void DiskSourceTerms(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &prim_scalar,
                     const AthenaArray<Real> &bcc,
                     AthenaArray<Real> &cons,
                     AthenaArray<Real> &cons_scalar) {
  AthenaArray<Real> le;
  if (enable_qw_source) {
    le.NewAthenaArray(12);
    LuminosityData(le);
  }
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real rho = prim(IDN,k,j,i);
        const Real ye = std::max(1.0e-12, std::min(1.0 - 1.0e-12,
                                                   prim_scalar(ye_index,k,j,i)));
        if (enable_qw_source) {
          Real r_cell[NSCALARS];
          for (int n=0; n<NSCALARS; ++n) r_cell[n] = prim_scalar(n,k,j,i);
          r_cell[ye_index] = ye;

          const Real temp = pmb->peos->TFromRhoP(rho, prim(IPR,k,j,i), r_cell);
          const Real temp_mev = temp*8.617333262e-11;
          const Real x = std::sqrt(std::max(0.0, 1.0 - SQR(source_radius/pmb->pcoord->x1v(i))));
          Real ypfree, ynfree;
          FreeNucleonFractions(pmb->peos, rho, temp, ye, &ypfree, &ynfree);
          const Real eta = QWEta(rho, temp, ye);
          Real qdot_erg = qdot_scale
              *QdotQW(temp_mev, rho, x, ypfree, ynfree, eta, le)*1.602176634e-6;
          if (rho < qw_rho_min) qdot_erg = std::max(0.0, qdot_erg);
          Real de = dt*rho*qdot_erg;
          const Real ekin = 0.5*(SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i))
                                 + SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
          const Real eint = std::max(cons(IEN,k,j,i) - ekin, 0.0);
          const Real de_lim = qw_energy_fraction_limit*eint;
          de = std::max(-de_lim, std::min(de_lim, de));
          cons(IEN,k,j,i) += de;

          const Real yedot = (rho < qw_rho_min && qdot_erg <= 0.0) ? 0.0
              : yedot_scale*YeSource(temp_mev, ye, ypfree, ynfree, x, rho, le);
          const Real dYe = dt*yedot;
          const Real ye_new = std::max(ye_floor, std::min(ye_ceil, ye + dYe));
          cons_scalar(ye_index,k,j,i) += rho*(ye_new - ye);
        }
        if (ye_source_enabled) {
          const Real dye = (ye_eq - ye)*(1.0 - std::exp(-dt/ye_tau));
          cons_scalar(ye_index,k,j,i) += cons(IDN,k,j,i)*dye;
        }
      }
    }
  }
}

Real FermiApprox(int n, Real eta) {
  if (n == 2) {
    Real a = std::exp(-std::abs(eta));
    Real s = std::pow(eta,3.0)/3.0 + 3.2898681*eta;
    Real ff = 2.0*(a - std::pow(a,2.0)/8.0 + std::pow(a,3.0)/27.0
                   - std::pow(a,4.0)/64.0 + std::pow(a,5.0)/125.0
                   - std::pow(a,6.0)/216.0);
    return eta < 0.0 ? ff : s + ff;
  }
  if (n == 3) {
    Real a = std::exp(-std::abs(eta));
    Real s = std::pow(eta,4.0)/4.0 + 4.9348022*SQR(eta) + 11.351273;
    Real ff = 6.0*(a - std::pow(a,2.0)/16.0 + std::pow(a,3.0)/81.0
                   - std::pow(a,4.0)/256.0);
    return eta < 0.0 ? ff : s - ff;
  }
  if (n == 4) {
    Real a = std::exp(-std::abs(eta));
    Real s = std::pow(eta,5.0)/5.0 + 6.5797363*std::pow(eta,3.0) + 45.457576*eta;
    Real ff = 24.0*(a - std::pow(a,2.0)/32.0 + std::pow(a,3.0)/243.0);
    return eta < 0.0 ? ff : s + ff;
  }
  if (n == 5) {
    Real a = std::exp(-std::abs(eta));
    Real s = std::pow(eta,6.0)/6.0 + 8.2246703*std::pow(eta,4.0)
             + 113.64394*SQR(eta) + 236.53226;
    Real ff = 120.0*(a - std::pow(a,2.0)/64.0 + std::pow(a,3.0)/729.0);
    return eta < 0.0 ? ff : s - ff;
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in FermiApprox" << std::endl
      << "Argument different from 2, 3, 4 or 5 passed." << std::endl;
  ATHENA_ERROR(msg);
  return 0.0;
}

void LuminosityData(AthenaArray<Real> &out) {
  const Real ferm5 = 118.266;
  const Real ferm4 = 23.3309;
  const Real ferm3 = 5.6822;
  const Real ferm2 = 1.80309;
  const Real ktemp_nu = eps_nu*ferm2/ferm3;
  const Real ktemp_nubar = eps_nubar*ferm2/ferm3;
  out(0) = l_nu;
  out(1) = l_nubar;
  out(2) = l_mu;
  out(3) = std::sqrt(SQR(ktemp_nu)*ferm5/ferm3);
  out(4) = std::sqrt(SQR(ktemp_nubar)*ferm5/ferm3);
  out(5) = eps_mu;
  out(6) = eps_nu;
  out(7) = SQR(ktemp_nu)*ferm4/ferm2;
  out(8) = eps_nubar;
  out(9) = SQR(ktemp_nubar)*ferm4/ferm2;
  out(10) = l_nu*SQR(ktemp_nu)*ferm5/ferm3;
  out(11) = l_nubar*SQR(ktemp_nubar)*ferm5/ferm3;
}

Real QWEta(Real rho, Real temp, Real ye) {
  const Real third = 1.0/3.0;
  const Real kb = 1.380649e-16;
  const Real mn = 1.6726e-24;
  const Real hbar = 6.62607015e-27/(2.0*PI);
  const Real c3 = std::pow(kb/(hbar*c_cgs), 3.0);
  const Real eta_den_const = std::pow(6.0, 2.0*third);
  const Real eta_den_a = std::pow(2.0, third)/eta_den_const;
  const Real eta_den_b = 2.0*std::pow(3.0, third)/eta_den_const;
  const Real vol = mn/rho;
  const Real a = c3*std::pow(temp,3.0)*vol*(PI/3.0);
  const Real a2 = SQR(a);
  const Real a4 = SQR(a2);
  const Real a6 = a2*a4;
  const Real b = std::sqrt(4.0*a6 + 27.0*a4*SQR(ye));
  const Real term = std::pow(9.0*a2*ye + std::sqrt(3.0)*b, third);
  return PI*(eta_den_a*term/a - eta_den_b*a/term);
}

void FreeNucleonFractions(EquationOfState *peos, Real rho, Real temp, Real ye,
                          Real *ypfree, Real *ynfree) {
  const Real xalpha = peos->XalphaFromRhoTYe(rho, temp, ye);
  *ypfree = std::max(0.0, std::min(1.0, ye - 0.5*xalpha));
  *ynfree = std::max(0.0, std::min(1.0, 1.0 - ye - 0.5*xalpha));
}

Real QdotQW(Real temp, Real rho, Real x, Real ypfree, Real ynfree,
            Real etaqw, AthenaArray<Real> &le) {
  const Real Na = 6.02214076e23;
  const Real ferm6 = 714.668;
  const Real ferm5 = 118.266;
  const Real ferm4 = 23.3309;
  const Real ferm3 = 5.6822;
  const Real ferm2 = 1.80309;
  const Real lnu_y = le(0);
  const Real lnubar_y = le(1);
  const Real lmu_y = le(2);
  const Real epsmu_y = le(5);
  const Real eavg1_nu = le(6);
  const Real eavg1_nubar = le(8);
  const Real coeff_nu = le(10);
  const Real coeff_nubar = le(11);

  Real h1 = 1e12*9.65*Na*(ynfree*coeff_nu + ypfree*coeff_nubar)
            *(1.0 - x)*source_radius_inv2;
  const Real f2e = FermiApprox(2, etaqw);
  const Real f3e = FermiApprox(3, etaqw);
  const Real f4e = FermiApprox(4, etaqw);
  const Real f2m = FermiApprox(2, -etaqw);
  const Real f3m = FermiApprox(3, -etaqw);
  const Real f4m = FermiApprox(4, -etaqw);
  const Real f1 = 1e8*(1.511/rho)*Na*std::pow(temp,4.0);
  Real h2_nubar = f1*lnubar_y*(0.926*f2e*(eavg1_nubar*ferm2/ferm3 - 4.0*temp*f3e/f4e)
                   + 2.20858*f2m*(eavg1_nubar*ferm2/ferm3 - 4.0*temp*f3m/f4m));
  Real h2_nu = f1*lnu_y*(2.20858*f2e*(eavg1_nu*ferm2/ferm3 - 4.0*temp*f3e/f4e)
                + 0.926*f2m*(eavg1_nu*ferm2/ferm3 - 4.0*temp*f3m/f4m));
  Real h2_mu = 2.0*f1*lmu_y*(f2e*(0.360592*epsmu_y*ferm2/ferm3 - 4.0*temp*f3e/f4e)
                 + 0.31*f2m*(epsmu_y*ferm2/ferm3 - 4.0*temp*f3m/f4m));
  Real h2_mubar = 2.0*f1*lmu_y*(f2e*(0.31*epsmu_y*ferm2/ferm3 - 4.0*temp*f3e/f4e)
                    + 0.360592*f2m*(epsmu_y*ferm2/ferm3 - 4.0*temp*f3m/f4m));
  Real h2 = 1e12*(h2_nubar + h2_nu + h2_mu + h2_mubar)
            *(1.0 - x)*source_radius_inv2;
  const Real mn = 1.6726e-24;
  Real fact = (lnubar_y*SQR(eavg1_nubar)*(eavg1_nubar*ferm2/ferm3 - 6.0*temp*ferm5/ferm6)
              + lnu_y*SQR(eavg1_nu)*(eavg1_nu*ferm2/ferm3 - 6.0*temp*ferm5/ferm6)
              + 4.0*lmu_y*SQR(epsmu_y)*(epsmu_y*ferm2/ferm3 - 6.0*temp*ferm5/ferm6))/mn;
  Real h3 = 6.2415e5*(5.2225e4*ynfree*fact + 4.3215e4*ypfree*fact)
            *(1.0 - x)*source_radius_inv2;
  Real psi = (x*x + 4.0*x + 5.0)*std::pow(1.0 - x, 4.0);
  Real h4 = (1.6e19*6.2415e5*1e32*psi/(rho*std::pow(source_radius,4.0)))
            *(lnu_y*lnubar_y*(eavg1_nubar + eavg1_nu) + (6.0/7.0)*lmu_y*lmu_y*epsmu_y);
  const Real heating = (std::isfinite(h1) ? h1 : 0.0)
                     + (std::isfinite(h2) ? h2 : 0.0)
                     + (std::isfinite(h3) ? h3 : 0.0)
                     + (std::isfinite(h4) ? h4 : 0.0);

  const Real cf = ypfree*(FermiApprox(5, etaqw)/ferm5)
                + ynfree*(FermiApprox(5, -etaqw)/ferm5);
  Real c1 = 2.073*Na*cf*std::pow(temp,6.0);
  Real c2 = (0.145*Na*1e8*std::pow(temp,9.0)/rho)
            *((f4e*f3m + f4m*f3e)/(2.0*ferm4*ferm3));
  const Real cooling = (std::isfinite(c1) ? c1 : 0.0) + (std::isfinite(c2) ? c2 : 0.0);
  return heating - cooling;
}

Real YeSource(Real temp, Real ye, Real ypfree, Real ynfree,
              Real x, Real rho, AthenaArray<Real> &le) {
  const Real alpha = 1.254;
  const Real gf = 1.16637e-11;
  const Real delta = 1.2935;
  const Real hbar = 6.582119569e-22;
  const Real erg2mev = 6.24151e5;
  const Real ferm4 = 23.3309;
  const Real lnu_y = le(0);
  const Real lnubar_y = le(1);
  const Real eavg1_nu = le(6);
  const Real eavg2_nu = le(7);
  const Real eavg1_nubar = le(8);
  const Real eavg2_nubar = le(9);
  Real lambda_nue_n = SQR(hbar*c_cgs)*((1.0 + 3.0*SQR(alpha))/(2.0*PI*PI))
      *SQR(gf)*lnu_y*(1e51*erg2mev*source_radius_inv2)
      *(eavg2_nu/eavg1_nu + 2.0*delta + SQR(delta)/eavg1_nu)*(1.0 - x);
  Real lambda_nuebar_p = SQR(hbar*c_cgs)*((1.0 + 3.0*SQR(alpha))/(2.0*PI*PI))
      *SQR(gf)*lnubar_y*(1e51*erg2mev*source_radius_inv2)
      *(eavg2_nubar/eavg1_nubar - 2.0*delta + SQR(delta)/eavg1_nubar)*(1.0 - x);
  if (!std::isfinite(lambda_nue_n)) lambda_nue_n = 0.0;
  if (!std::isfinite(lambda_nuebar_p)) lambda_nuebar_p = 0.0;

  const Real eta = QWEta(rho, temp/8.617333262e-11, ye);
  const Real cc = 0.448*std::pow(temp,5.0)/ferm4;
  Real lambda_ele_p = cc*FermiApprox(4, eta);
  Real lambda_pos_n = cc*FermiApprox(4, -eta);
  if (!std::isfinite(lambda_ele_p)) lambda_ele_p = 0.0;
  if (!std::isfinite(lambda_pos_n)) lambda_pos_n = 0.0;
  return (lambda_nue_n + lambda_pos_n)*ynfree
         - (lambda_nuebar_p + lambda_ele_p)*ypfree;
}

void DiskInnerX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                      FaceField &b, Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        const Real r1 = pco->x1v(il-i);
        const Real r2 = pco->x1v(il+i-1);
        prim(IDN,k,j,il-i) = prim(IDN,k,j,il+i-1);
        prim(IVX,k,j,il-i) = std::min(0.0, prim(IVX,k,j,il+i-1)*SQR(r2/r1));
        prim(IVY,k,j,il-i) = prim(IVY,k,j,il+i-1);
        prim(IVZ,k,j,il-i) = prim(IVZ,k,j,il+i-1)*r1/r2;
        prim(IPR,k,j,il-i) = prim(IPR,k,j,il+i-1);
        for (int n=0; n<NSCALARS; ++n) {
          pmb->pscalars->r(n,k,j,il-i) = pmb->pscalars->r(n,k,j,il+i-1);
        }
      }
    }
  }
}

void DiskOuterX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                      FaceField &b, Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        const Real r1 = pco->x1v(iu+i);
        const Real r2 = pco->x1v(iu-i+1);
        prim(IDN,k,j,iu+i) = prim(IDN,k,j,iu-i+1);
        prim(IVX,k,j,iu+i) = std::max(0.0, prim(IVX,k,j,iu-i+1)*SQR(r2/r1));
        prim(IVY,k,j,iu+i) = prim(IVY,k,j,iu-i+1);
        prim(IVZ,k,j,iu+i) = prim(IVZ,k,j,iu-i+1)*r1/r2;
        prim(IPR,k,j,iu+i) = prim(IPR,k,j,iu-i+1);
        for (int n=0; n<NSCALARS; ++n) {
          pmb->pscalars->r(n,k,j,iu+i) = pmb->pscalars->r(n,k,j,iu-i+1);
        }
      }
    }
  }
}

} // namespace
