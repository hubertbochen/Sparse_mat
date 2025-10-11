#include "lib\data.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm> // std::max
#include <omp.h>     // OpenMP for omp_get_wtime

// exact solution and its Laplacian for u = exp(x*y)
static inline double exact_u(double x, double y) {
    return std::exp(x*y-x) + y;
}
static inline double laplacian_exact(double x, double y) {
    return (x*x + (y-1)*(y-1)) * std::exp(x*y-x);
}

int main() {
     //record all terminal output to a log file
     //both output in terminal and saved to file
     //to disable, comment out the following line

    freopen("output.log", "w", stdout);
    std::printf("2D Poisson solver on L-shaped domain using FIM and sparse matrix\n");
    std::printf("C++17, OpenMP %d threads\n", omp_get_max_threads());
    std::printf("Assemble laplacian matrix and solve using Gauss-Seidel iteration\n");
    std::printf("Exact solution: u = exp(x*y/4) + y, Laplacian = (x^2 + y^2)*exp(x*y/4)/16\n");
    std::printf("Domain: L-shape with vertices (0,-2),(0,2),(2,2),(2,1),(1,1),(1,-2)\n");
    std::printf("Mesh: uniform grid with nx by ny points, spacing hx, hy\n");
    std::printf("Boundary condition: Dirichlet u=0 on outer boundary\n");
    std::printf("Inner boundary (the cut-out triangle): u = exact solution\n");
    std::printf("------------------------------------------------------------\n");
    // Set up mesh and data


    // choose k and enforce nx = ny = 4*k + 1
    const int k = 64;               // change k for resolution
    const int nx = 2 * k + 1;
    const int ny = 4 * k + 1;

    const double x0 = 0.0, x1 = 2.0;      // x in [0,2]
    const double y0 = -2.0, y1 = 2.0;     // y in [-2,2]
    printf("Mesh size: nx = %d, ny = %d\n", nx, ny);
    printf("Domain: x in [%.2f, %.2f], y in [%.2f, %.2f]\n", x0, x1, y0, y1);
    Data data;
    data.mesh.init(nx, ny, x1 - x0, y1 - y0);

    data.mesh.pt_type = new int*[nx];
    for (int i = 0; i < nx; ++i) data.mesh.pt_type[i] = new int[ny];

    double hx = data.mesh.hx;
    double hy = data.mesh.hy;

    // Use index-based detection since nx = 4*k + 1
    // index mapping: x=0 -> i=0, x=1 -> i=2k, x=2 -> i=4k
    //                y=-2 -> j=0, y=-1 -> j=k, y=0 -> j=2k, y=1 -> j=3k, y=2 -> j=4k

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            // default inner
            data.mesh.pt_type[i][j] = 3;

            // Outer rectangle boundary
            if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
                data.mesh.pt_type[i][j] = 1; // boundary
            }

            // Triangle 1: vertices indices (0,3k), (0,4k), (2k,4k)
            // condition (integer form): j >= 3k and 2*j >= 6k + i
            if ( (j > 3*k) && (j > 3*k + i) ) {
                data.mesh.pt_type[i][j] = 0; // removed
                continue;
            }

             if ( (j >= 3*k) && (j == 3*k + i) ) {
                data.mesh.pt_type[i][j] = 1; 
                continue;
            }

            // Triangle 2: vertices indices (2k,4k), (4k,4k), (4k,3k)
            // condition: i >= 2k and 2*j >= 10k - i  (from j >= 5k - 0.5*i)
            if ( (i > k) && (j > 5*k - i) ) {
                data.mesh.pt_type[i][j] = 0; // removed
                continue;
            }

          if ( (i >= k) && (j == 5*k - i) ) {
                data.mesh.pt_type[i][j] = 1; 
                continue;
            }


            // Triangle 3: vertices indices (2k,k), (4k,3k), (4k,0)
            // region defined by i >= 2k and 4k - i <= 2*j <= 2*i - 2k
            if (i > k) {
                int lhs = 4*k - 2*i;         // lower bound *2
                int rhs = 4*i - 2*k;      // upper bound *2
                int twoj = 2*j;
                if (twoj > lhs && twoj < rhs) {
                    data.mesh.pt_type[i][j] = 0; // removed
                    continue;
                }
            }

            if (i >= k) {
                int lhs = 4*k - 2*i;         // lower bound *2
                int rhs = 4*i - 2*k;      // upper bound *2
                int twoj = 2*j;
                if (twoj == lhs || twoj == rhs) {
                    data.mesh.pt_type[i][j] = 1; // removed
                    continue;
                }
            }
        }
    }

    // promote inner points adjacent to removed to near-boundary (2)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (data.mesh.pt_type[i][j] == 3) {
                bool adj_nonexist = false;
                const int di[1] = {1};
                const int dj[1] = {0};
                for (int m = 0; m < 1; ++m) {
                    int ni = i + di[m], nj = j + dj[m];
                    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny) {
                        if (data.mesh.pt_type[ni][nj] == 0) { adj_nonexist = true; break; }
                    }
                }
                if (adj_nonexist) data.mesh.pt_type[i][j] = 2;
            }
        }
    }
    // display pt_type for verification
    if (k <= 16) {
        std::printf("Point types (0:non-exist, 1:boundary, 2:near-boundary, 3:inner):\n");
        for (int j = ny - 1; j >= 0; --j)
        {
            for (int i = 0; i < nx; ++i) {
                std::printf("%d  ", data.mesh.pt_type[i][j]);
            }
            std::printf("\n");
        }
    }
    data.mesh.build_list(&data.mesh);
    data.mesh.build_stiff(&data.mesh);

    int size = data.mesh.allocated_size;
    data.sfunc = new double[size];
    data.solution = new double[size];
    for (int kidx = 0; kidx < size; ++kidx) { data.sfunc[kidx] = 0.0; data.solution[kidx] = 0.0; }

    // set RHS using exact solution: solver uses positive Laplacian Δu = f
    for (int lin = 0; lin < size; ++lin) {
        int ii = data.mesh.xlist[lin];
        int jj = data.mesh.ylist[lin];
        double x = x0 + ii * hx;
        double y = y0 + jj * hy;
        if (data.mesh.pt_type[ii][jj] == 1) {
            data.sfunc[lin] = exact_u(x,y); // Dirichlet
        
        } else {
            data.sfunc[lin] = -laplacian_exact(x,y); // Δu
        }
    }
    // record calculation time
    double start_time = omp_get_wtime();
    data.Solve_Poisson();
    double end_time = omp_get_wtime();
    std::printf("Calculation time: %.6f seconds\n", end_time - start_time);
    //data.print_solution();

    // compute error vs exact
    double err_max = 0.0;
    double err_l2 = 0.0;
    for (int lin = 0; lin < size; ++lin) {
        int i = data.mesh.xlist[lin];
        int j = data.mesh.ylist[lin];
        double x = x0 + i * hx;
        double y = y0 + j * hy;
        double u_ex = exact_u(x,y);
        double u_num = data.solution[lin];
        double e = u_num - u_ex;
        err_max = std::max(err_max, fabs(e));
        err_l2 += e*e;
    }
    err_l2 = std::sqrt(err_l2 / size);
    std::printf("Error: max = %.6e, L2 = %.6e\n", err_max, err_l2);

    data.save_solution_to_file("results/solution.txt");
    // free memory

  
    



      // change stdout to terminal 
    freopen("CON", "w", stdout); // Use "CON" for Windows, "/dev/tty" for Unix
    std::fflush(stdout);
    std::printf("All done.\n");
    //print the file output.log to terminal
    std::FILE* logFile = std::fopen("output.log", "r");
    if (logFile) {
        char buffer[256];
        while (std::fgets(buffer, sizeof(buffer), logFile)) {
            std::printf("%s", buffer);
        }
        std::fclose(logFile);
    } else {
        std::printf("Error: Could not open log file for reading\n");
    }

   

    return 0;


}