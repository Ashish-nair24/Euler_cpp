
#pragma once
#include <vector>
#include "grid.hpp"

namespace cfd {

struct FieldSoA {
  std::vector<double> rho, rhou, rhov, E;
  Grid2D g;
  explicit FieldSoA(const Grid2D& grid) : g(grid) {
    const auto n = g.nx_tot() * g.ny_tot();
    rho.resize(n, 0.0);
    rhou.resize(n, 0.0);
    rhov.resize(n, 0.0);
    E.resize(n, 0.0);
  }
  inline std::size_t idx(std::size_t i, std::size_t j) const { return g.idx(i, j); }
};

} // namespace cfd
