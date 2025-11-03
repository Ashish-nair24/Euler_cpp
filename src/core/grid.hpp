
#pragma once
#include <cstddef>
#include <cassert>

namespace cfd {

struct Grid2D {
  std::size_t nx{}, ny{};
  double dx{1.0}, dy{1.0};
  std::size_t ng{1};
  std::size_t nx_tot() const { return nx + 2*ng; }
  std::size_t ny_tot() const { return ny + 2*ng; }
  std::size_t idx(std::size_t i, std::size_t j) const {
    assert(i < nx_tot() && j < ny_tot());
    return j * nx_tot() + i;
  }
};

} // namespace cfd
