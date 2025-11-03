
#pragma once
#include "euler.hpp"
namespace cfd {
// Piecewise-constant reconstruction
inline void reconstruct_pc(const Cons* U, int i, Cons& UL, Cons& UR) {
  UL = U[i];
  UR = U[i+1];
}
} // namespace cfd
