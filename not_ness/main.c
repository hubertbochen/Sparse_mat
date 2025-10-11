#include "Spmat.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#include "dense.h"
#include "mesh.h"
#include "data.h"

#endif
int main() {
    SpMat M, T;
    int N = 10; // Matrix size
    if (CreateSpMat(&M, N, N)) {
        printf("Sparse matrix created successfully.\n");
        printf("Rows: %d, Columns: %d\n", M.m, M.n);
    } else {
        printf("Failed to create sparse matrix.\n");
    }
   
    //test append
   // order-1 Central Difference Matrix
   
    int rows[N*2+1], cols[N*2+1];
    double data[N*2+1];
    int idx = 0;
    for (int i = 1; i <= N; i++) {
        if (i > 1) { rows[idx] = i; cols[idx] = i-1; data[idx] = -1.0; idx++; } // left neighbor
        if (i < N) { rows[idx] = i; cols[idx] = i+1; data[idx] = 1.0; idx++; }  // right neighbor
    }
    rows[idx] = N; cols[idx] = N; data[idx] = 1.0; idx++;
    rows[idx] = N; cols[idx] = N-1; data[idx] = 1.0; idx++;
    rows[idx] = 1; cols[idx] = 1; data[idx] = -1.0; idx++;


    if (mat_append(&M, rows, cols, data, idx-1)) {
        printf("Data appended successfully.\n");
    } else {
        printf("Failed to append data.\n");
    }
    printf("Non-zero elements:\n");
    for (int k = 0; k < M.tu; k++) {
        printf("Row: %d, Column: %d, Data: %.2f\n", M.elem[k].i, M.elem[k].j, M.elem[k].data);
    }   
    //test transpose
    if (mat_transpose(&M, &T)) {
        printf("Matrix transposed successfully.\n");
        printf("Transposed matrix non-zero elements:\n");
        for (int k = 0; k < T.tu; k++) {
            printf("Row: %d, Column: %d, Data: %.2f\n", T.elem[k].i, T.elem[k].j, T.elem[k].data);
        }
    } else {
        printf("Failed to transpose matrix.\n");
    }   
    //test addition

    SpMat A;
    if (mat_copy(&M, &A)) {
        printf("Matrices copied successfully.\n");
        printf("Resultant matrix non-zero elements:\n");
        for (int k = 0; k < A.tu; k++) {
            printf("Row: %d, Column: %d, Data: %.2f\n", A.elem[k].i, A.elem[k].j, A.elem[k].data);
        }
    } else {
        printf("Failed to add matrices.\n");
    }
    //test multiplication
    SpMat D, E, F;
    CreateSpMat(&D, 2, 3);
    CreateSpMat(&E, 3, 2);
    int rowsD[] = {1, 1, 2};
    int colsD[] = {1, 3, 2};
    double dataD[] = {1.0, 2.0, 3.0};
    mat_append(&D, rowsD, colsD, dataD, 3);
    int rowsE[] = {1, 2, 3};
    int colsE[] = {2, 1, 2};
    double dataE[] = {4.0, 5.0, 6.0};
    mat_append(&E, rowsE, colsE, dataE, 3);
    if (mat_multiply(&D, &E, &F)) {
        printf("Matrices multiplied successfully.\n");
        printf("Resultant matrix non-zero elements:\n");
        for (int k = 0; k < F.tu; k++) {
            printf("Row: %d, Column: %d, Data: %.2f\n", F.elem[k].i, F.elem[k].j, F.elem[k].data);
        }
    } else {
        printf("Failed to multiply matrices.\n");
    }

    //test LU decomposition and solve
    SpMat L, U;
    return 0;
}