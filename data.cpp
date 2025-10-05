#include "data.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

// If Spmat.h doesn't declare Guass_Seidel, declare it here
extern "C" {
    int Guass_Seidel(SpMat* A, double* b, double* x, int max_iter, double tol);
}

Data::Data() : sfunc(nullptr), solution(nullptr) {}
Data::~Data() { delete[] sfunc; delete[] solution; }

void Data::init(int nx, int ny, double xlim, double ylim) {
    mesh.init(nx, ny, xlim, ylim);

    mesh.pt_type = new int*[mesh.nx];
    for (int i = 0; i < mesh.nx; ++i) {
        mesh.pt_type[i] = new int[mesh.ny];
        for (int j = 0; j < mesh.ny; ++j) {
            if (i == 0 || j == 0 || i == mesh.nx - 1 || j == mesh.ny - 1)
                mesh.pt_type[i][j] = 1; // boundary
            else
                mesh.pt_type[i][j] = 3; // inner
        }
    }

    mesh.build_list(&mesh);
    mesh.build_stiff(&mesh);

    int size = mesh.allocated_size;
    sfunc = new double[size];
    solution = new double[size];
    for (int k = 0; k < size; ++k) { sfunc[k] = 0.0; solution[k] = 0.0; }

    std::printf("Data initialized with mesh size %dx%d, %d active points\n", nx, ny, size);
}

void Data::set_boundary_conditions() {
    int size = mesh.allocated_size;
    for (int lin = 0; lin < size; ++lin) {
        int i = mesh.xlist[lin];
        int j = mesh.ylist[lin];
        if (mesh.pt_type[i][j] == 1) sfunc[lin] = 0.0;
        else sfunc[lin] = 1.0;
    }
    std::printf("Boundary conditions and source function set\n");
}

void Data::Solve_Poisson() {
    if (!sfunc || !solution) { std::printf("Error: Data not initialized\n"); return; }
    int n = mesh.allocated_size;
    if (n <= 0) { std::printf("Error: mesh has no active points\n"); return; }

    int max_iter = 100000;
    double tol = 1e-10;
    std::printf("Solving with Gauss-Seidel: n=%d, max_iter=%d, tol=%g\n", n, max_iter, tol);

    Gauss_Seidel(&mesh.stiff, sfunc, solution, max_iter, tol);

    double res_norm = 0.0;
    for (int idx = 0; idx < n; ++idx) {
        double Ax_i = 0.0;
        for (int k = 0; k < mesh.stiff.tu; ++k) {
            if (mesh.stiff.elem[k].i - 1 == idx) {
                int col = mesh.stiff.elem[k].j - 1;
                Ax_i += mesh.stiff.elem[k].data * solution[col];
            }
        }
        double r = sfunc[idx] - Ax_i;
        res_norm += r * r;
    }
    res_norm = std::sqrt(res_norm);
    std::printf("Gauss-Seidel finished, residual norm = %.6e\n", res_norm);
}

void Data::print_solution() {
    if (!solution) { std::printf("Error: No solution to print\n"); return; }
    std::printf("Solution:\n");
    for (int lin = 0; lin < mesh.allocated_size; ++lin) {
        int i = mesh.xlist[lin];
        int j = mesh.ylist[lin];
        std::printf("Point (%d,%d): u = %.6f\n", i, j, solution[lin]);
    }
}

void Data::save_solution_to_file(const char* filename) {
    if (!solution || mesh.allocated_size <= 0) {
        std::printf("Error: no solution to save\n"); return;
    }
    FILE* file = fopen(filename, "w");
    if (!file) { std::printf("Error: Could not open file %s for writing\n", filename); return; }
    for (int lin = 0; lin < mesh.allocated_size; ++lin) {
        int i = mesh.xlist[lin];
        int j = mesh.ylist[lin];
        fprintf(file, "%d %d %.12e\n", i, j, solution[lin]);
    }
    fclose(file);
    std::printf("Solution saved to %s\n", filename);
}

double Data::get_solution_at(int i, int j) {
    if (i < 0 || i >= mesh.nx || j < 0 || j >= mesh.ny) return 0.0;
    if (mesh.linlist[i][j] == -1) return 0.0;
    int lin = mesh.linlist[i][j];
    return solution[lin];
}