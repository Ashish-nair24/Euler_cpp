
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#ifdef USE_MPI
  #include <mpi.h>
#endif
#include "core/grid.hpp"
#include "core/field_soa.hpp"
#include "core/vtu_writer.hpp"
#include "fv_euler2d/euler.hpp"
#include "fv_euler2d/rk3.hpp"
#include "io/csv_writer.hpp"

static int to_int(const char* s, int defv) {
  if (!s) return defv;
  try { return std::stoi(s); } catch (...) { return defv; }
}
static double to_double(const char* s, double defv) {
  if (!s) return defv;
  try { return std::stod(s); } catch (...) { return defv; }
}

int main(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  (void)size;  
#else
  int rank = 0;
#endif

  std::string problem = "sod1d";
  int nx = 400, ng = 2;
  double tend = 0.2, cfl = 0.5;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--help" || a == "-h") {
      if (rank == 0) {
        std::cout << "euler2d (beginner starter)\n"
                  << "  --problem sod1d\n"
                  << "  --nx N     (default 400)\n"
                  << "  --tend T   (default 0.2)\n"
                  << "  --cfl  C   (default 0.5)\n";
      }
#ifdef USE_MPI
      MPI_Finalize();
#endif
      return 0;
    } else if (a == "--problem" && i+1 < argc) {
      problem = argv[++i];
    } else if (a == "--nx" && i+1 < argc) {
      nx = to_int(argv[++i], nx);
    } else if (a == "--tend" && i+1 < argc) {
      tend = to_double(argv[++i], tend);
    } else if (a == "--cfl" && i+1 < argc) {
      cfl = to_double(argv[++i], cfl);
    }
  }

  if (problem != "sod1d") {
    if (rank == 0) std::cerr << "Only --problem sod1d is implemented.\n";
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  const int n_tot = nx + 2*ng;
  const double dx = 1.0 / nx;
  std::vector<cfd::Cons> U(n_tot);

  // Sod initial condition
  const int mid = nx / 2;
  for (int i = 0; i < n_tot; ++i) {
    const int icell = i - ng;
    cfd::Prim W;
    if (icell < mid) { W = {1.0, 0.0, 0.0, 1.0}; }
    else             { W = {0.125, 0.0, 0.0, 0.1}; }
    U[i] = cfd::prim2cons(W);
  }

  auto apply_bc = [&](std::vector<cfd::Cons>& A) {
    A[0] = A[ng];
    A[1] = A[ng];
    A[n_tot-2] = A[n_tot-ng-1];
    A[n_tot-1] = A[n_tot-ng-1];
  };

  double t = 0.0; int iter = 0;
  while (t < tend && iter < 200000) {
    std::vector<cfd::Cons> R(n_tot);
    double amax = cfd::residual_1d(U, R, dx);
    const double dt = cfl * dx / (amax + 1e-12);
    apply_bc(U);
    cfd::advance_rk3_1d(U, dt, dx);
    t += dt; ++iter;
  }

  if (rank == 0) {
    cfd::write_sod_csv("sod_density.csv", U, ng, dx);
    std::cout << "Done. Wrote sod_density.csv (x,rho,u,p)\n";
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
