/*Class Data declaration and methods for managing mesh and solving Poisson equation
variables:
    mesh: Mesh2D object representing the computational mesh
    sfunc: array representing the source function (right-hand side)
    solution: array representing the solution vector
methods:
    init: initialize the Data object with mesh dimensions and allocate arrays
    set_boundary_conditions: set boundary conditions and source function values
    Solve_Poisson: solve the Poisson equation using Gauss-Seidel method
    print_solution: print the solution array to standard output
    save_solution_to_file: save the solution array to a text file
    get_solution_at: get the solution value at a specific (i,j) grid point

*/

#ifndef DATA_H
#define DATA_H

#include "mesh.hpp"
#include "Spmat.h"

#ifdef __cplusplus
class Data {
public:
    Mesh2D mesh;
    double* sfunc;    // Source function (right-hand side)
    double* solution; // Solution vector

    // Constructor and destructor
    Data();
    ~Data();

    // Methods
    void init(int nx, int ny, double xlim, double ylim);
    void set_boundary_conditions();
    void Solve_Poisson();
    void print_solution();
	void save_solution_to_file(const char* filename);
    double get_solution_at(int i, int j);
};
#endif

#endif // DATA_H
