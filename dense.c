// structure of dense matrix
#include "Spmat.h"
#include <stdio.h>
#define MaxNum 1000
typedef struct {
    int m;          // number of rows
    int n;          // number of columns
    double data[MaxNum][MaxNum]; // data array
} Dense;
int CreateDense(Dense* D, int m, int n);
int dense_add(Dense* A, Dense* B, Dense* C);
int dense_multiply(Dense* A, Dense* B, Dense* C);
int dense_print(Dense* D);
int dense_to_sparse(Dense* D, SpMat* M);
int sparse_to_dense(SpMat* M, Dense* D);

int CreateDense(Dense* D, int m, int n) {
    if (m < 1 || n < 1 || m > MaxNum || n > MaxNum) return 0;
    D->m = m;
    D->n = n;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            D->data[i][j] = 0.0; // Initialize to zero
        }
    }
    return 1;
}

int dense_add(Dense* A, Dense* B, Dense* C) {
    if (A->m != B->m || A->n != B->n) return 0; // Dimension mismatch
    if (!CreateDense(C, A->m, A->n)) return 0;
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            C->data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
    return 1;
}

int dense_multiply(Dense* A, Dense* B, Dense* C) {
    if (A->n != B->m) return 0; // Dimension mismatch
    if (!CreateDense(C, A->m, B->n)) return 0;
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < B->n; j++) {
            C->data[i][j] = 0.0;
            for (int k = 0; k < A->n; k++) {
                C->data[i][j] += A->data[i][k] * B->data[k][j];
            }
        }
    }
    return 1;
}

int dense_print(Dense* D) {
    if (!D) return 0;
    for (int i = 0; i < D->m; i++) {
        for (int j = 0; j < D->n; j++) {
            printf("%8.2f ", D->data[i][j]);
        }
        printf("\n");
    }
    return 1;
}

int dense_to_sparse(Dense* D, SpMat* M) {
    if (!D || !M) return 0;
    if (!CreateSpMat(M, D->m, D->n)) return 0;
    for (int i = 0; i < D->m; i++) {
        for (int j = 0; j < D->n; j++) {
            if (D->data[i][j] != 0.0) {
                if (M->tu >= MaxNum) return 0; // No space to add new element
                M->elem[M->tu].i = i + 1; // Convert to 1-based index
                M->elem[M->tu].j = j + 1; // Convert to 1-based index
                M->elem[M->tu].data = D->data[i][j];
                M->tu++;
            }
        }
    }
    return 1;
}

int sparse_to_dense(SpMat* M, Dense* D) {
    if (!M || !D) return 0;
    if (!CreateDense(D, M->m, M->n)) return 0;
    for (int k = 0; k < M->tu; k++) {
        int row = M->elem[k].i - 1; // Convert to 0-based index
        int col = M->elem[k].j - 1; // Convert to 0-based index
        D->data[row][col] = M->elem[k].data;
    }
    return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Test the functions
