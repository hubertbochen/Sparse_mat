/**
 * @brief Data class for managing mesh and solving Poisson equation 
 * @file data.hpp
 * @author Chen Bojin
 * @date 2025-10-13
 * @note Uses C++17 standard
 * @note Depends on mesh.hpp and Spmat.h
 * @version 0.1
 * @note This class encapsulates a 2D mesh, source function, and solution vector. It provides methods to initialize the mesh,
 *       set boundary conditions, and solve the Poisson equation.
 *      The Poisson equation is solved using the Gauss-Seidel iterative method implemented in Spmat.c.
 *      The mesh is represented by the Mesh2D class, and the sparse matrix operations are handled by the SpMat structure.
 *      Memory management is handled in the constructor and destructor.
 */






#include "mesh.hpp"


// Expose to Doxygen even if __cplusplus isn’t defined

/**
 * @brief Manages mesh, RHS and solution; solves Poisson using Gauss–Seidel.
 */
class Data {
public:
    Mesh2D mesh;
    double* sfunc;/**< source function */
    double* solution;/**< solution vector */

    Data();
    ~Data();

    void init(int nx, int ny, double xlim, double ylim);/**< initialize mesh and vectors */
    void set_boundary_conditions();/**< set Dirichlet boundary conditions */
    void Solve_Poisson();/**< solve Poisson equation */
    void print_solution();/**< print solution vector */
    void save_solution_to_file(const char* filename);/**< save solution to a text file */
    double get_solution_at(int i, int j);/**< get solution at mesh point (i,j) */
};



