# cfd-lite (beginner-friendly starter)

This is a **minimal C++20 CFD portfolio** starter with a working **Sod 1D** shock tube using HLL flux and TVD RK3, plus CSV output you can plot.

## Build & run
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
./euler2d --problem sod1d --nx 400 --tend 0.2 --cfl 0.5
# Output: sod_density.csv  (columns: x,rho,u,p)
```
Open `sod_density.csv` in Python/Excel to visualize.

## What to change next
- Swap piecewise-constant reconstruction for MUSCL (limiters).
- Add 2D grid and x/y-face loops.
- Introduce MPI row-wise halos later using `src/parallel/`.
