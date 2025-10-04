
#ifndef DATA_H
#define DATA_H

#include "mesh.h"
#include "Spmat.h"
#include <cstddef>

#ifdef __cplusplus
class Data {
public:
	Mesh2D mesh;
	int* sfunc;
	int* solution;

	Data();
	~Data();

	void Solve_Poisson(Data* data);
	// Add methods for data management as needed
};
#endif

#endif // DATA_H
