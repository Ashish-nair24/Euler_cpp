
#pragma once
#include <cmath>
#include <algorithm>

namespace cfd {

struct Cons { double rho, rhou, rhov, E; };
struct Prim { double rho, u, v, p; };

inline constexpr double gamma_gas = 1.4;

inline double pressure_from_cons(const Cons& U) {
  const double u = (U.rho > 0.0) ? (U.rhou / U.rho) : 0.0;
  const double v = (U.rho > 0.0) ? (U.rhov / U.rho) : 0.0;
  const double kinetic = 0.5 * U.rho * (u*u + v*v);
  const double p = (gamma_gas - 1.0) * std::max(0.0, U.E - kinetic);
  return p;
}

inline Prim cons2prim(const Cons& U) {
  Prim W;
  W.rho = U.rho;
  W.u   = (U.rho > 0.0) ? (U.rhou / U.rho) : 0.0;
  W.v   = (U.rho > 0.0) ? (U.rhov / U.rho) : 0.0;
  W.p   = pressure_from_cons(U);
  return W;
}

inline Cons prim2cons(const Prim& W) {
  Cons U;
  U.rho  = W.rho;
  U.rhou = W.rho * W.u;
  U.rhov = W.rho * W.v;
  const double e = W.p / (gamma_gas - 1.0) + 0.5 * W.rho * (W.u*W.u + W.v*W.v);
  U.E = e;
  return U;
}

inline Cons flux_x(const Cons& U) {
  const Prim W = cons2prim(U);
  Cons F;
  F.rho  = W.rho * W.u;
  F.rhou = W.rho * W.u * W.u + W.p;
  F.rhov = W.rho * W.u * W.v;
  F.E    = (U.E + W.p) * W.u;
  return F;
}

inline double a_sound(const Prim& W) {
  return std::sqrt(gamma_gas * std::max(1e-12, W.p) / std::max(1e-12, W.rho));
}

} // namespace cfd
