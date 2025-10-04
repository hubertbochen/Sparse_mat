#include <stdio.h>
#include "dense.h"
#include "Spmat.h"

class Mesh2D {
public:
    int nx;
    int ny;
    int allocated_size;
    int xlim;
    int ylim;
    int hx;// = (xlim)/(nx-1);
    int hy;// = (ylim)/(ny-1);
    int** pt_type; // non_exist: 0, boundary: 1, nearboundary :2 ,inner: 3
    int* xlist;    // xlist[lin] to find x index
    int* ylist;    // ylist[lin] to find y index
    int** linlist;  // linlist[i][j] to find the lin index
    SpMat stiff;

    Mesh2D()
        : nx(0), ny(0), allocated_size(0), xlim(0), ylim(0), hx(0), hy(0), 
          pt_type(nullptr), xlist(nullptr), ylist(nullptr), linlist(nullptr), stiff() {}

    void init(int nx_, int ny_, double xlim_, double ylim_) {
        nx = nx_;
        ny = ny_;
        xlim = static_cast<int>(xlim_);
        ylim = static_cast<int>(ylim_);
        if (nx > 1 && ny > 1) {
            hx = xlim / (nx - 1);
            hy = ylim / (ny - 1);
        } else {
            fprintf(stderr, "Error: nx and ny must be greater than 1.\n");
            hx = 0;
            hy = 0;
        }
    }

    ~Mesh2D() {
        delete[] xlist;
        delete[] ylist;
        if (pt_type) {
            for (int i = 0; i < nx; i++) {
                delete[] pt_type[i];
            }
            delete[] pt_type;
        }
        if (linlist) {
            for (int i = 0; i < nx; i++) {
                delete[] linlist[i];
            }
            delete[] linlist;
        }
    }

    // Build index lists based on pt_type
    void build_list(Mesh2D* mesh) {
        pt_type = mesh->pt_type;
        nx = mesh->nx;
        ny = mesh->ny;
        xlim = mesh->xlim;
        ylim = mesh->ylim;
        allocated_size = nx * ny;
        xlist = new int[allocated_size];
        ylist = new int[allocated_size];
        linlist = new int*[nx];
        for (int i = 0; i < nx; i++) {
            linlist[i] = new int[ny];
        }
        
        // Initialize all to -1 first
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                linlist[i][j] = -1;
            }
        }
        
        // Then set valid indices
        int count = 0;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                if (pt_type[i][j] != 0) {
                    xlist[count] = i;
                    ylist[count] = j;
                    linlist[i][j] = count;
                    count++;
                }
            }
        }
        allocated_size = count;
    }

    // Build stiffness matrix
    void build_stiff(Mesh2D* mesh) {
        stiff.m = allocated_size;
        stiff.n = allocated_size;
        stiff.tu = 0;
        hx = mesh->hx;
        hy = mesh->hy;

        for (int lin = 0; lin < allocated_size; lin++) {
            int i = xlist[lin];
            int j = ylist[lin];
            
            if (stiff.tu >= MaxNum) break;
            
            double temp = 0.0;
            switch (static_cast<int>(mesh->pt_type[i][j])) {
                case 0: // non_exist
                    break;

                case 1: // boundary
                    stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].j = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].data = 1.0;
                    stiff.tu++;
                    break;

                case 2: // nearboundary
                    temp = 0.0;
                    if (i > 0 && mesh->pt_type[i - 1][j] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i - 1][j] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hx;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (i < nx - 1 && mesh->pt_type[i + 1][j] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i + 1][j] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hx;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (j > 0 && mesh->pt_type[i][j - 1] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i][j - 1] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hy;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (j < ny - 1 && mesh->pt_type[i][j + 1] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i][j + 1] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hy;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }

                    stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].j = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].data = -temp;
                    stiff.tu++;
                    break;

                case 3: // inner
                    temp = 0.0;
                    if (i > 0 && mesh->pt_type[i - 1][j] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i - 1][j] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hx;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (i < nx - 1 && mesh->pt_type[i + 1][j] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i + 1][j] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hx;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (j > 0 && mesh->pt_type[i][j - 1] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i][j - 1] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hy;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }
                    if (j < ny - 1 && mesh->pt_type[i][j + 1] >= 1) {
                        stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                        stiff.elem[stiff.tu].j = linlist[i][j + 1] + 1;
                        stiff.elem[stiff.tu].data = -1.0/hy;
                        temp += stiff.elem[stiff.tu].data;
                        stiff.tu++;
                    }

                    stiff.elem[stiff.tu].i = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].j = linlist[i][j] + 1;
                    stiff.elem[stiff.tu].data = -temp;
                    stiff.tu++;
                    break;
            }
        }
    }
};

int main() {
    printf("Mesh2D class defined.\n");
    
        Mesh2D mesh;
        mesh.init(10, 10, 1.0, 1.0);
        mesh.pt_type = new int*[mesh.nx];
        for (int i = 0; i < mesh.nx; i++) {
            mesh.pt_type[i] = new int[mesh.ny];
            for (int j = 0; j < mesh.ny; j++) {
                if (i == 0 || j == 0 || i == mesh.nx - 1 || j == mesh.ny - 1) {
                    mesh.pt_type[i][j] = 1; // boundary
                }  else {
                    mesh.pt_type[i][j] = 3; // inner
                }
            }
        }
        mesh.build_list(&mesh);
        mesh.build_stiff(&mesh);
        SpMat stiff = mesh.stiff;
        printf("Stiffness matrix has %d non-zero elements.\n", stiff.tu);
        return 0;
}