#pragma once
#include <cstdio>
#include <vector>
#include <string>
#include "fv_euler2d/euler.hpp"

namespace cfd {

inline void write_sod_csv(const std::string& path,
                          const std::vector<Cons>& U,
                          int ng, double dx) {
  FILE* f = std::fopen(path.c_str(), "w");
  if (!f) return;
  std::fprintf(f, "x,rho,u,p\n");
  const int n = static_cast<int>(U.size());
  for (int i = ng; i < n - ng; ++i) {
    const Prim W = cons2prim(U[i]);
    const double x = (i - ng + 0.5) * dx;
    std::fprintf(f, "%.8f,%.8f,%.8f,%.8f\n", x, W.rho, W.u, W.p);
  }
  std::fclose(f);
}

} // namespace cfd

