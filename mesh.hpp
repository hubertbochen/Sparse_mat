/*
Mesh2D class declaration and C interface for mesh operations
variables:
    nx, ny: number of grid points in x and y directions
    allocated_size: number of active points (non-zero type)
    xlim, ylim: physical dimensions of the mesh
    hx, hy: grid spacing in x and y directions
    pt_type: 2D array indicating point types (0: non_exist, 1: boundary, 2: nearboundary, 3: inner)
    xlist, ylist: arrays mapping linear index to (i,j) coordinates
    linlist: 2D array mapping (i,j) coordinates to linear index
    stiff: sparse matrix representing the Laplacian operator
methods:
    init: initialize the mesh with given dimensions and spacing
    build_list: build the mapping arrays (xlist, ylist, linlist)
    build_stiff: assemble the sparse laplacian matrix using finite difference stencils
*/

#ifndef MESH_H
#define MESH_H

#include "Spmat.h"

#ifdef __cplusplus
extern "C" {
#endif

// C declarations for mesh (add C structs/functions here if needed)

#ifdef __cplusplus
}

// C++ Mesh2D class declaration
class Mesh2D {
public:
    int nx;
    int ny;
    int allocated_size;
    double xlim;
    double ylim;
    double hx;
    double hy;
    int** pt_type;    // 2D array. non_exist: 0, boundary: 1, nearboundary: 2, inner: 3
    int* xlist;       // xlist[lin] to find x index
    int* ylist;       // ylist[lin] to find y index
    int** linlist;    // linlist[i][j] to find the lin index
    SpMat stiff;      // Laplacian matrix

    // Constructor and destructor
    Mesh2D();
    ~Mesh2D();

    // Methods
    void init(int nx_, int ny_, double xlim_, double ylim_);
    void build_list(Mesh2D* mesh);
    void build_stiff(Mesh2D* mesh);
};
#endif

#endif // MESH_H
