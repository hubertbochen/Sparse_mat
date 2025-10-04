#include "mesh.h"
#include "Spmat.h"
#include <cstdio>

class Data {
public:
    Mesh2D mesh;
    int* sfunc;
    int* solution;

    Data() : sfunc(nullptr), solution(nullptr) {}

    ~Data() {
        delete[] sfunc;
        delete[] solution;
    }


    void Solve_Poisson(Data* data) {
        // Implement the Poisson solver using mesh and SpMat
        
    }
    // Add methods for data management as needed
};

int main() {
    printf("Data class defined.\n");
    return 0;
}