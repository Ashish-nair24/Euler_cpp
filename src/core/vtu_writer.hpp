
#pragma once
#include <cstdio>
#include <string>
#include "grid.hpp"

namespace cfd {

inline void write_csv_stub(const std::string& path, const Grid2D& g) {
  FILE* f = std::fopen(path.c_str(), "w");
  if (!f) return;
  std::fprintf(f, "i,j,x,y\n");
  for (std::size_t j=0;j<g.ny_tot();++j) {
    for (std::size_t i=0;i<g.nx_tot();++i) {
      const double x = (static_cast<double>(i)-g.ng+0.5)*g.dx;
      const double y = (static_cast<double>(j)-g.ng+0.5)*g.dy;
      std::fprintf(f, "%zu,%zu,%.6f,%.6f\n", i,j,x,y);
    }
  }
  std::fclose(f);
}

} // namespace cfd
