/**
 * @file Spmat.h
 * @author Chen Bojin
 * @brief Header file for sparse matrix operations
 * @version 0.1
 * @date 2025-10-13
 */

#ifndef SPMAT_H
#define SPMAT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#define MaxNum 1000000

typedef struct triple {
    /**< triple structure for sparse matrix element */
    int i;/**< row index */
    int j;/**< column index */
    double data;/**< non-zero value */
} triple, *MatHead;

typedef struct SpMat {
    /**< sparse matrix structure */
    triple* elem; /**< dynamically allocated array of triples */
    int m;/**< number of rows */
    int n;/**< number of columns */
    int tu;/**< number of non-zero elements */
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
int mat_LU_solve(SpMat* L, SpMat* U, double* b, double* x);
int Jacobi(SpMat* A, double* b, double* x, int max_iter, double tol);
int precompute_row_indexs(SpMat* A, int** row_ptr);
int Gauss_Seidel(SpMat* A, double* b, double* x, int max_iter, double tol);
static int compare_triples(const void* a, const void* b);

#ifdef __cplusplus
}
#endif

#endif // SPMAT_H