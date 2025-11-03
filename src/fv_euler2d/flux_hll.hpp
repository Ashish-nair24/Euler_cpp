
#pragma once
#include <utility>
#include <algorithm>
#include "euler.hpp"

namespace cfd {

inline std::pair<Cons,double> flux_hll(const Cons& UL, const Cons& UR) {
  const Prim WL = cons2prim(UL);
  const Prim WR = cons2prim(UR);
  const double aL = a_sound(WL);
  const double aR = a_sound(WR);

  const double sL = std::min(WL.u - aL, WR.u - aR);
  const double sR = std::max(WL.u + aL, WR.u + aR);

  const Cons FL = flux_x(UL);
  const Cons FR = flux_x(UR);

  Cons FH{};
  double amax = 0.0;
  if (sL >= 0.0) {
    FH = FL;
    amax = std::max(std::abs(sL), std::abs(sR));
  } else if (sR <= 0.0) {
    FH = FR;
    amax = std::max(std::abs(sL), std::abs(sR));
  } else {
    const double inv = 1.0 / (sR - sL);
    FH.rho  = (sR*FL.rho  - sL*FR.rho  + sR*sL*(UR.rho  - UL.rho )) * inv;
    FH.rhou = (sR*FL.rhou - sL*FR.rhou + sR*sL*(UR.rhou - UL.rhou)) * inv;
    FH.rhov = (sR*FL.rhov - sL*FR.rhov + sR*sL*(UR.rhov - UL.rhov)) * inv;
    FH.E    = (sR*FL.E    - sL*FR.E    + sR*sL*(UR.E    - UL.E   )) * inv;
    amax = std::max(std::abs(sL), std::abs(sR));
  }
  return {FH, amax};
}

} // namespace cfd
