#include "mesh.h"
#include <cstdio>
#include <cstdlib>


// Constructor
Mesh2D::Mesh2D()
    : nx(0), ny(0), allocated_size(0), xlim(0), ylim(0),
      hx(0), hy(0), pt_type(nullptr), xlist(nullptr), ylist(nullptr),
      linlist(nullptr), stiff()
{}

// Destructor
Mesh2D::~Mesh2D() {
    delete[] xlist;
    delete[] ylist;
    if (pt_type) {
        for (int i = 0; i < nx; ++i) delete[] pt_type[i];
        delete[] pt_type;
    }
    if (linlist) {
        for (int i = 0; i < nx; ++i) delete[] linlist[i];
        delete[] linlist;
    }
}

// init
void Mesh2D::init(int nx_, int ny_, double xlim_, double ylim_) {
    nx = nx_;
    ny = ny_;
    xlim = static_cast<int>(xlim_);
    ylim = static_cast<int>(ylim_);
    if (nx > 1 && ny > 1) {
        hx = xlim / double(nx - 1);
        hy = ylim / double(ny - 1);
    } else {
        std::fprintf(stderr, "Error: nx and ny must be greater than 1.\n");
        hx = hy = 0;
    }
}

// build_list
void Mesh2D::build_list(Mesh2D* mesh) {
    // free previous
    delete[] xlist;
    delete[] ylist;
    if (linlist) {
        for (int i = 0; i < nx; ++i) delete[] linlist[i];
        delete[] linlist;
        linlist = nullptr;
    }

    pt_type = mesh->pt_type;
    nx = mesh->nx;
    ny = mesh->ny;
    xlim = mesh->xlim;
    ylim = mesh->ylim;

    allocated_size = nx * ny;
    xlist = new int[allocated_size];
    ylist = new int[allocated_size];

    linlist = new int*[nx];
    for (int i = 0; i < nx; ++i) {
        linlist[i] = new int[ny];
        for (int j = 0; j < ny; ++j) linlist[i][j] = -1;
    }

    int count = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (pt_type[i][j] != 0) {
                xlist[count] = i;
                ylist[count] = j;
                linlist[i][j] = count;
                ++count;
            }
        }
    }
    allocated_size = count;
}

// build_stiff
void Mesh2D::build_stiff(Mesh2D* mesh) {
    if (!CreateSpMat(&stiff, allocated_size, allocated_size)) {
        std::fprintf(stderr, "Error: CreateSpMat failed\n");
        return;
    }

    stiff.tu = 0;
    hx = mesh->hx;
    hy = mesh->hy;

    for (int lin = 0; lin < allocated_size && stiff.tu < MaxNum; ++lin) {
        int i = xlist[lin];
        int j = ylist[lin];

        switch (mesh->pt_type[i][j]) {
            case 0:
                break;
            case 1: { // Dirichlet boundary
                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = lin + 1;
                stiff.elem[stiff.tu].data = 1.0;
                ++stiff.tu;
                break;
            }
            case 2: {
                double diag = 0.0;
                // left
              
            
                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = linlist[i-1][j] + 1;
                stiff.elem[stiff.tu].data = -2.0 / (hx*hx);
                diag -= stiff.elem[stiff.tu].data;
                ++stiff.tu; 
             
                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = linlist[i+1][j+1] + 1;
                stiff.elem[stiff.tu].data = -1.0 / (hx*hy);
                diag -= stiff.elem[stiff.tu].data;
                ++stiff.tu;

                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = linlist[i][j-1] + 1;
                stiff.elem[stiff.tu].data = -2.0 / (hy*hy);
                diag -= stiff.elem[stiff.tu].data;
                ++stiff.tu; 

                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = linlist[i-1][j-1] + 1;
                stiff.elem[stiff.tu].data = 1.0 / (hy*hx);
                diag -= stiff.elem[stiff.tu].data;
                ++stiff.tu; 
                // diagonal
                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = lin + 1;
                stiff.elem[stiff.tu].data = diag;
                ++stiff.tu;
                break;
         
            }

            case 3: { // nearboundary / inner: 5-point stencil
                double diag = 0.0;
                // left
                if (i > 0 && linlist[i-1][j] != -1) {
                    stiff.elem[stiff.tu].i = lin + 1;
                    stiff.elem[stiff.tu].j = linlist[i-1][j] + 1;
                    stiff.elem[stiff.tu].data = -1.0 / (hx*hx);
                    diag -= stiff.elem[stiff.tu].data;
                    ++stiff.tu;
                }
                // right
                if (i < nx-1 && linlist[i+1][j] != -1) {
                    stiff.elem[stiff.tu].i = lin + 1;
                    stiff.elem[stiff.tu].j = linlist[i+1][j] + 1;
                    stiff.elem[stiff.tu].data = -1.0 / (hx*hx);
                    diag -= stiff.elem[stiff.tu].data;
                    ++stiff.tu;
                }
                // down
                if (j > 0 && linlist[i][j-1] != -1) {
                    stiff.elem[stiff.tu].i = lin + 1;
                    stiff.elem[stiff.tu].j = linlist[i][j-1] + 1;
                    stiff.elem[stiff.tu].data = -1.0 / (hy*hy);
                    diag -= stiff.elem[stiff.tu].data;
                    ++stiff.tu;
                }
                // up
                if (j < ny-1 && linlist[i][j+1] != -1) {
                    stiff.elem[stiff.tu].i = lin + 1;
                    stiff.elem[stiff.tu].j = linlist[i][j+1] + 1;
                    stiff.elem[stiff.tu].data = -1.0 / (hy*hy);
                    diag -= stiff.elem[stiff.tu].data;
                    ++stiff.tu;
                }
                // diagonal
                stiff.elem[stiff.tu].i = lin + 1;
                stiff.elem[stiff.tu].j = lin + 1;
                stiff.elem[stiff.tu].data = diag;
                ++stiff.tu;
                break;
            }
        }
    }
    std::printf("Stiffness matrix built with %d elements\n", stiff.tu);
}