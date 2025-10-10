#ifndef DENSE_H
#define DENSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "Spmat.h"
#include <stdio.h>

#define MaxNum 1000

typedef struct Dense {
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

#ifdef __cplusplus
}
#endif

#endif // DENSE_H
