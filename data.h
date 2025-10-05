#ifndef DATA_H
#define DATA_H

#include "mesh.h"
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
