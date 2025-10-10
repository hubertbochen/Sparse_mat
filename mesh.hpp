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
    int** pt_type;    // 2D array: non_exist: 0, boundary: 1, nearboundary: 2, inner: 3
    int* xlist;       // xlist[lin] to find x index
    int* ylist;       // ylist[lin] to find y index
    int** linlist;    // linlist[i][j] to find the lin index
    SpMat stiff;      // Stiffness matrix

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
