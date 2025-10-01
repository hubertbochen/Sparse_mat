
#ifndef SPMAT_H
#define SPMAT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#define MaxNum 1000

typedef struct triple {
    int i;
    int j;
    double data;
} triple, *MatHead;

typedef struct SpMat {
    triple elem[MaxNum];
    int m;
    int n;
    int tu;
} SpMat;

int CreateSpMat(SpMat* M, int m, int n);
int mat_append(SpMat* M, int rows[], int cols[], double data[], int count);
int mat_add(SpMat* A, SpMat* B, SpMat* C);
int mat_multiply(SpMat* A, SpMat* B, SpMat* C);
int mat_transpose(SpMat* M, SpMat* T);
int mat_LU_decompose(SpMat* A, SpMat* L, SpMat* U);
int mat_diag(SpMat* A, SpMat* D);
int mat_rmzero(SpMat* M);
int mat_print(SpMat* M);
int mat_upper_tri(SpMat* A, SpMat* U);
int mat_lower_tri(SpMat* A, SpMat* L);
int mat_get(SpMat* M, int row, int col, double* value);
int mat_set(SpMat* M, int row, int col, double value);
int mat_to_dense(SpMat* M, double** dense);
int mat_copy(SpMat* src, SpMat* dest);

#ifdef __cplusplus
}
#endif

#endif // SPMAT_H