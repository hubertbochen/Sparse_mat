/**
 * @brief Mesh2D class declaration and C interface for mesh operations
 * @file mesh.hpp
 * @author Chen Bojin
 * @date 2025-10-13
 * @note Uses C++17 standard
 * @note Depends on Spmat.h
 * @version 0.1
 * @note This file declares the Mesh2D class for managing a tensor 2D mesh and its associated sparse matrix representation of the Laplacian operator.
 *      The class provides methods to initialize the mesh, build the list of active points, and construct the sparse Laplacian matrix.
 *     The mesh is represented using a structured grid with point types indicating boundary and inner points.
 *    Memory management is handled in the constructor and destructor.
 *    The file also includes C declarations for potential C functions related to mesh operations.
 *    The Mesh2D class is designed to be used in conjunction with the Data class for solving PDEs on the mesh.
 *    The sparse matrix operations are handled by the SpMat structure defined in Spmat.h.
 * 
 */




#include "Spmat.h"

#ifdef __cplusplus
extern "C" {
#endif

// C declarations for mesh (add C structs/functions here if needed)

#ifdef __cplusplus
} // extern "C"
#endif



/**
 * @brief Tensor 2D mesh and its sparse Laplacian matrix.
 *
 * Manages grid geometry, active-point lists, and builds the Laplacian (SpMat).
 * pt_type: 0=non_exist, 1=boundary, 2=near-boundary, 3=inner.
 */
class Mesh2D {
public:
    int nx;/**< number of grid points in x direction */
    int ny;/**< number of grid points in y direction */
    int allocated_size;/**< allocated size of the mesh */
    double xlim;/**< x-axis limits */
    double ylim;/**< y-axis limits */
    double hx;/**< grid spacing in x direction */
    double hy;/**< grid spacing in y direction */
    int** pt_type;/**< point types: 0=non_exist, 1=boundary, 2=near-boundary, 3=inner */
    int* xlist;/**< list of active x coordinates */
    int* ylist;/**< list of active y coordinates */
    int** linlist;/**< list of active point indices */
    SpMat stiff;/**< sparse Laplacian matrix */

    Mesh2D();
    ~Mesh2D();

    void init(int nx_, int ny_, double xlim_, double ylim_);/**< initialize mesh and its sparse matrix */
    void build_list(Mesh2D* mesh);/**< build list of active points */
    void build_stiff(Mesh2D* mesh);/**< build sparse Laplacian matrix */
};



