#ifndef MESH_H
#define MESH_H
#include "Spmat.h"
#ifdef __cplusplus
extern "C" {
#endif

// C declarations for mesh (add C structs/functions here)

#ifdef __cplusplus
}

// C++ Mesh2D class declaration
class Mesh2D {
public:
	int nx;
	int ny;
	int xlim;
	int ylim;
	int* xlist;
	int* ylist;
	int* linlist;
	struct SpMat stiff;
     //struct Dense pt_type;// non_exist: 0 , boundary:1 , inner :2
     Mesh2D();
	~Mesh2D();
	void build_list();
	void build_stiff();
};
#endif

#endif // MESH_H
