# Sparse Test: 2D Poisson with CSR, Gauss–Seidel/Jacobi

This project solves 2D Poisson problems on structured grids using sparse matrices:
- L-shaped domain with triangle cut-outs (non-convex)
- Full rectangle (naive) with Dirichlet boundary

Core features:
- CSR-like sparse storage and operations in C (Spmat.c/Spmat.h)
- Assembly of Laplacian on active grid points (mesh.hpp/cpp)
- Iterative solvers: Gauss–Seidel (GS) and Jacobi (OpenMP-enabled)
- Verification against an analytic solution; plotting scripts; LaTeX report

## Project structure

- main.cpp                Entry, examples (L-shape + naive rectangle)
- lib\Spmat.h             Sparse matrix C API (triplets, solvers, utils)
- lib\mesh.hpp            Mesh class interface
- lib\data.hpp            Data class (wraps mesh + RHS + solution)
- source\mesh.cpp         Mesh2D implementation (build lists, assemble)
- source\data.cpp         Data methods (Solve_Poisson, IO, etc.)
- source\Spmat.c          Sparse ops + GS/Jacobi (C, OpenMP)
- results\                Outputs (solution.txt, naive_solution.txt, ...)
- plot.py, naive.py       Matplotlib visualization scripts
- report\report1.tex      LaTeX report

## Requirements

- Windows + MSYS2 UCRT64 (recommended) or MinGW-w64
- C++17 compiler + OpenMP (e.g., g++ with -fopenmp)
- Python 3 with numpy, matplotlib (for plotting)
- LaTeX (optional, to build report)
- Doxygen (optional, to build code docs)

## Build (Windows, MSYS2 UCRT64)

1) Open “MSYS2 UCRT64” shell, or ensure g++.exe on PATH.
2) Create folders:
   - mkdir -p build results
3) Compile:
   D:\env\msys2\ucrt64\bin\g++.exe -std=c++17 -O2 -fopenmp ^
     main.cpp ^
     source\mesh.cpp source\data.cpp source\Spmat.c ^
     -I lib ^
     -o build\poisson.exe

Notes:
- You can also use g++ from MSYS2 UCRT64 without absolute path if on PATH.
- Add -g for debugging, or -static only if your toolchain supports it.

## Run

From the workspace root:
- build\poisson.exe

Behavior:
- Logs run info to output.log and echoes it to terminal at the end.
- L-shaped example saves results\solution.txt (x y u).
- Naive rectangle example saves:
  - results\naive_solution.txt (x y u)
  - results\naive_rhs.txt      (x y f)

Adjust resolution:
- Edit k in main.cpp. The code sets nx, ny from k (e.g., nx=2k+1, ny=4k+1).
- Larger k increases DOFs and nnz; GS iterations will take longer.

## Plot

- L-shaped:
  python plot.py
- Rectangle (naive):
  python naive.py

Tips:
- We triangulate and mask triangles that lie in cut-out regions to avoid convex-hull “repair”.
- Ensure results\*.txt exist before plotting.

## Analytic solution (default)

- u(x,y) = exp(x*y - x) + y
- Δu(x,y) = (x^2 + (y-1)^2) * exp(x*y - x)

RHS sign:
- Keep RHS consistent with your assembly. Current main.cpp sets:
  - boundary rows: b = u_exact
  - interior rows: b = -Δu_exact  (see comments in main.cpp)

## Notes on solvers

- Gauss–Seidel (GS): single pass over nonzeros per sweep; in-place updates.
- Jacobi: single pass accumulate row sums, parallelized with OpenMP.
- Missing diagonal entries are replaced by a small epsilon with a warning; prefer assembling a valid diagonal.

## Documentation

- Report: report\report1.tex (Chinese), includes pt_type=2 near-slope special stencil.
- Doxygen (optional): create a Doxyfile and run `doxygen Doxyfile` to index lib/*.hpp and lib/Spmat.h.

## Troubleshooting

- Link errors about Guass_Seidel vs Gauss_Seidel: a compatibility wrapper is provided; keep names consistent.
- Large nx, ny: ensure no out-of-bounds in mesh assembly and that MaxNum in Spmat.h is large enough for nnz.
- OpenMP threads: controlled by OMP_NUM_THREADS env var or defaults; code reports detected threads.

## License and authors

- Author: Chen Bojin
- Educational use for numerical PDEs and sparse linear algebra
