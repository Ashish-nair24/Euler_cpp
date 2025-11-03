
#pragma once
#include <vector>
#include <algorithm>
#include "euler.hpp"
#include "recon.hpp"
#include "flux_hll.hpp"

namespace cfd {

inline double residual_1d(const std::vector<Cons>& U, std::vector<Cons>& R, double dx) {
  const int n = static_cast<int>(U.size());
  std::fill(R.begin(), R.end(), Cons{0,0,0,0});
  double amax = 0.0;
  for (int i = 0; i < n-1; ++i) {
    Cons UL, UR;
    reconstruct_pc(U.data(), i, UL, UR);
    auto [F, a] = flux_hll(UL, UR);
    amax = std::max(amax, a);
    R[i].rho   -= F.rho  / dx;
    R[i].rhou  -= F.rhou / dx;
    R[i].rhov  -= F.rhov / dx;
    R[i].E     -= F.E    / dx;

    R[i+1].rho   += F.rho  / dx;
    R[i+1].rhou  += F.rhou / dx;
    R[i+1].rhov  += F.rhov / dx;
    R[i+1].E     += F.E    / dx;
  }
  return amax;
}

inline void advance_rk3_1d(std::vector<Cons>& U, double dt, double dx) {
  const int n = static_cast<int>(U.size());
  std::vector<Cons> R(n), U1(n), U2(n);

  residual_1d(U, R, dx);
  for (int i = 0; i < n; ++i) {
    U1[i].rho  = U[i].rho  + dt * R[i].rho;
    U1[i].rhou = U[i].rhou + dt * R[i].rhou;
    U1[i].rhov = U[i].rhov + dt * R[i].rhov;
    U1[i].E    = U[i].E    + dt * R[i].E;
  }

  residual_1d(U1, R, dx);
  for (int i = 0; i < n; ++i) {
    U2[i].rho  = 0.75*U[i].rho  + 0.25*(U1[i].rho  + dt * R[i].rho);
    U2[i].rhou = 0.75*U[i].rhou + 0.25*(U1[i].rhou + dt * R[i].rhou);
    U2[i].rhov = 0.75*U[i].rhov + 0.25*(U1[i].rhov + dt * R[i].rhov);
    U2[i].E    = 0.75*U[i].E    + 0.25*(U1[i].E    + dt * R[i].E);
  }

  residual_1d(U2, R, dx);
  for (int i = 0; i < n; ++i) {
    U[i].rho  = (1.0/3.0)*U[i].rho  + (2.0/3.0)*(U2[i].rho  + dt * R[i].rho);
    U[i].rhou = (1.0/3.0)*U[i].rhou + (2.0/3.0)*(U2[i].rhou + dt * R[i].rhou);
    U[i].rhov = (1.0/3.0)*U[i].rhov + (2.0/3.0)*(U2[i].rhov + dt * R[i].rhov);
    U[i].E    = (1.0/3.0)*U[i].E    + (2.0/3.0)*(U2[i].E    + dt * R[i].E);
  }
}

} // namespace cfd
