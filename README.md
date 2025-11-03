# CFD From Scratch (C++20)

This is an ongoing personal project where I'm building **computational fluid dynamics (CFD) solvers from the ground up** in modern **C++20**.

The goal is to strengthen my understanding of **numerical methods**, **HPC**, and **software design for scientific computing** by implementing CFD solvers incrementally — from simple 1D problems to multi-dimensional, parallel, and GPU-accelerated systems.

---

## 🚀 Current Progress

✅ **Sod Shock Tube (1D Euler equations)**  
- Finite-Volume formulation with HLL flux  
- 3rd-order TVD Runge–Kutta time integration  
- Simple outflow boundary conditions  
- CSV output (`sod_density.csv`) for visualization  

Run it:
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
./euler2d --problem sod1d --nx 400 --tend 0.2 --cfl 0.5
```

🧩 Next Steps
In the coming updates, this project will expand toward a full 2D compressible flow solver:
MUSCL-type reconstruction with slope limiters (Van Leer / Venkatakrishnan)
HLLC flux for sharper contact resolution
Local time stepping and adaptive CFL control

💻 Parallel & GPU Acceleration
Row-wise domain decomposition with MPI halo exchange
Shared-memory parallelism using OpenMP
CUDA backend for flux evaluation on GPUs

🧠 Learning Focus
Finite-Volume methods for hyperbolic PDEs
Stability, conservation, and flux consistency
Modern C++ (RAII, templates, header-only design)
Parallel computing and performance tuning (MPI, CUDA, OpenMP)
