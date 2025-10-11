/*Sparse Matrix Library

Author: ChenBojin
Date: 2025-10-03
Version: 1.2
Description: A C library for sparse matrix operations, including creation, addition,
multiplication, transposition, and conversion between sparse and dense formats.

*/

#include "..\lib\Spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int CreateSpMat(SpMat* M, int m, int n) {
    /*Initialize a sparse matrix with given dimensions
    Input:
        M: pointer to the sparse matrix structure to be initialized
        m: number of rows
        n: number of columns
    Output:
        Returns 1 on success, 0 on failure (invalid dimensions or memory allocation failure)
    */
    if (m < 1 || n < 1 || m > MaxNum || n > MaxNum) return 0;
    M->m = m;
    M->n = n;
    M->tu = 0;
    M->elem = (triple*)malloc(MaxNum * sizeof(triple));
    if (!M->elem) return 0; // 内存分配失败
    return 1;
}

int mat_append(SpMat* M, int rows[], int cols[], double data[], int count) {
    /*Append non-zero elements to the sparse matrix
    Input:
        M: pointer to the sparse matrix structure
        rows: array of row indices (1-based)
        cols: array of column indices (1-based)
        data: array of non-zero values
        count: number of elements to append
    Output:
        Returns 1 on success, 0 on failure (invalid indices or exceeding storage capacity)
    */
    if (!M || !rows || !cols || !data || count < 1) return 0;
    int success = 1;
    if (count > 1000) {
        #pragma omp parallel for default(none) shared(M, rows, cols, data, count) reduction(&&:success)
        for (int k = 0; k < count; k++) {
            if (rows[k] < 1 || rows[k] > M->m || cols[k] < 1 || cols[k] > M->n) { success = 0; continue; }
            if (data[k] != 0) {
                int idx;
                #pragma omp critical
                {
                    if (M->tu >= MaxNum) { success = 0; }
                    else { idx = M->tu++; }
                }
                if (success) {
                    M->elem[idx].i = rows[k];
                    M->elem[idx].j = cols[k];
                    M->elem[idx].data = data[k];
                }
            }
        }
        return success;
    } else {
        for (int k = 0; k < count; k++) {
            if (rows[k] < 1 || rows[k] > M->m || cols[k] < 1 || cols[k] > M->n) return 0;
            if (data[k] != 0) {
                if (M->tu >= MaxNum) return 0; // 超出存储空间
                M->elem[M->tu].i = rows[k];
                M->elem[M->tu].j = cols[k];
                M->elem[M->tu].data = data[k];
                M->tu++;
            }
        }
        return 1;
    }
}

int mat_add(SpMat* A, SpMat* B, SpMat* C) {
    /*Add two sparse matrices A and B, store the result in C
    Input:
        A: pointer to the first sparse matrix
        B: pointer to the second sparse matrix
        C: pointer to the result sparse matrix
    Output:
        Returns 1 on success, 0 on failure (dimension mismatch or exceeding storage capacity)
    */
    if (A->m != B->m || A->n != B->n) return 0;
    C->m = A->m;
    C->n = A->n;
    C->tu = 0;
    
    // Use temporary storage for parallel processing
    int temp_count = 0;
    triple temp_result[MaxNum];
    
    if (A->tu + B->tu > 500) { // Parallel threshold
        #pragma omp parallel
        {
            int local_count = 0;
            triple local_result[MaxNum/4]; // Local storage per thread
            
            #pragma omp for schedule(dynamic)
            for (int pos = 0; pos < A->m * A->n; pos++) {
                int row = pos / A->n + 1;
                int col = pos % A->n + 1;
                
                double val_A = 0, val_B = 0;
                
                // Find value in A
                for (int k = 0; k < A->tu; k++) {
                    if (A->elem[k].i == row && A->elem[k].j == col) {
                        val_A = A->elem[k].data;
                        break;
                    }
                }
                
                // Find value in B
                for (int k = 0; k < B->tu; k++) {
                    if (B->elem[k].i == row && B->elem[k].j == col) {
                        val_B = B->elem[k].data;
                        break;
                    }
                }
                
                double sum = val_A + val_B;
                if (sum != 0 && local_count < MaxNum/4) {
                    local_result[local_count].i = row;
                    local_result[local_count].j = col;
                    local_result[local_count].data = sum;
                    local_count++;
                }
            }
            
            // Merge local results to global result
            #pragma omp critical
            {
                for (int i = 0; i < local_count && temp_count < MaxNum; i++) {
                    temp_result[temp_count++] = local_result[i];
                }
            }
        }
        
        // Copy results back
        C->tu = temp_count;
        for (int i = 0; i < temp_count; i++) {
            C->elem[i] = temp_result[i];
        }
        return 1;
    } else {
        // Keep your existing serial implementation
        int i = 0, j = 0, k = 0;
        while (i < A->tu && j < B->tu) {
            if (A->elem[i].i < B->elem[j].i || (A->elem[i].i == B->elem[j].i && A->elem[i].j < B->elem[j].j)) {
                C->elem[k++] = A->elem[i++];
                C->tu++;
            } else if (A->elem[i].i > B->elem[j].i || (A->elem[i].i == B->elem[j].i && A->elem[i].j > B->elem[j].j)) {
                C->elem[k++] = B->elem[j++];
                C->tu++;
            } else {
                double sum = A->elem[i].data + B->elem[j].data;
                if (sum != 0) {
                    C->elem[k].i = A->elem[i].i;
                    C->elem[k].j = A->elem[i].j;
                    C->elem[k++].data = sum;
                    C->tu++;
                }
                i++;
                j++;
            }
        }
        while (i < A->tu) { C->elem[k++] = A->elem[i++]; C->tu++; }
        while (j < B->tu) { C->elem[k++] = B->elem[j++]; C->tu++; }
        return 1;
    }
}

int mat_multiply(SpMat* A, SpMat* B, SpMat* C) {
    /*Multiply two sparse matrices A and B, store the result in C
    Input:
        A: pointer to the first sparse matrix
        B: pointer to the second sparse matrix
        C: pointer to the result sparse matrix
    Output:
        Returns 1 on success, 0 on failure (dimension mismatch or exceeding storage capacity)
    */  
    if (A->n != B->m) return 0; // Dimension mismatch
    C->m = A->m;
    C->n = B->n;
    C->tu = 0;
    
    // Use temporary storage for results
    triple temp_result[MaxNum];
    int temp_count = 0;
    
    if (A->m * B->n > 400) { // Parallel threshold
        #pragma omp parallel
        {
            int local_count = 0;
            triple local_result[MaxNum/4];
            
            #pragma omp for schedule(dynamic, 4)
            for (int idx = 0; idx < A->m * B->n; idx++) {
                int i = idx / B->n + 1;
                int j = idx % B->n + 1;
                double sum = 0.0;
                
                // Compute dot product of row i of A with column j of B
                for (int k = 1; k <= A->n; k++) {
                    double a_ik = 0, b_kj = 0;
                    
                    // Find A[i,k]
                    for (int p = 0; p < A->tu; p++) {
                        if (A->elem[p].i == i && A->elem[p].j == k) {
                            a_ik = A->elem[p].data;
                            break;
                        }
                    }
                    
                    // Find B[k,j]
                    for (int q = 0; q < B->tu; q++) {
                        if (B->elem[q].i == k && B->elem[q].j == j) {
                            b_kj = B->elem[q].data;
                            break;
                        }
                    }
                    
                    sum += a_ik * b_kj;
                }
                
                if (sum != 0 && local_count < MaxNum/4) {
                    local_result[local_count].i = i;
                    local_result[local_count].j = j;
                    local_result[local_count].data = sum;
                    local_count++;
                }
            }
            
            // Merge local results
            #pragma omp critical
            {
                for (int i = 0; i < local_count && temp_count < MaxNum; i++) {
                    temp_result[temp_count++] = local_result[i];
                }
            }
        }
        
        // Copy results back
        C->tu = temp_count;
        for (int i = 0; i < temp_count; i++) {
            C->elem[i] = temp_result[i];
        }
        return 1;
    } else {
        // Keep your existing serial implementation
        for (int i = 1; i <= A->m; i++) {
            for (int j = 1; j <= B->n; j++) {
                double sum = 0;
                for (int k = 1; k <= A->n; k++) {
                    double a_ik = 0, b_kj = 0;
                    for (int p = 0; p < A->tu; p++) {
                        if (A->elem[p].i == i && A->elem[p].j == k) {
                            a_ik = A->elem[p].data;
                            break;
                        }
                    }
                    for (int q = 0; q < B->tu; q++) {
                        if (B->elem[q].i == k && B->elem[q].j == j) {
                            b_kj = B->elem[q].data;
                            break;
                        }
                    }
                    sum += a_ik * b_kj;
                }
                if (sum != 0) {
                    if (C->tu >= MaxNum) return 0;
                    C->elem[C->tu].i = i;
                    C->elem[C->tu].j = j;
                    C->elem[C->tu].data = sum;
                    C->tu++;
                }
            }
        }
        return 1;
    }
}

int mat_transpose(SpMat* M, SpMat* T) {
    /*Transpose a sparse matrix M, store the result in T
    Input:
        M: pointer to the sparse matrix to be transposed
        T: pointer to the result sparse matrix
    Output:
        Returns 1 on success, 0 on failure (memory allocation failure)
    */
    (*T).m = M->n;
    (*T).n = M->m;
    (*T).tu = M->tu;

    if ((*T).tu) {
        int col, p;
        int q = 0;
  
        for (col = 1; col <= M->n; col++) {
       
            for (p = 0; p < M->tu; p++) {
           
                if (M->elem[p].j == col) {
                    (*T).elem[q].i = M->elem[p].j;
                    (*T).elem[q].j = M->elem[p].i;
                    (*T).elem[q].data = M->elem[p].data;
           
                    q++;
                }
            }
        }
    }
    return 1;
}

int mat_get(SpMat* M, int row, int col, double* value) {
    /*Get the value at specified row and column in the sparse matrix
    Input:
        M: pointer to the sparse matrix
        row: row index (1-based)
        col: column index (1-based)
        value: pointer to store the retrieved value
    Output:
        Returns 1 on success, 0 on failure (invalid indices or null pointer)
    */
    if (!M || row < 1 || row > M->m || col < 1 || col > M->n || !value) return 0;
    for (int k = 0; k < M->tu; k++) {
        if (M->elem[k].i == row && M->elem[k].j == col) {
            *value = M->elem[k].data;
            return 1;
        }
    }
    *value = 0.0; // Not found, return zero
    return 1;
}

int mat_set(SpMat* M, int row, int col, double value) {
    /*Set the value at specified row and column in the sparse matrix
    Input:
        M: pointer to the sparse matrix
        row: row index (1-based)
        col: column index (1-based)
        value: value to set
    Output:
        Returns 1 on success, 0 on failure (invalid indices or exceeding storage capacity)
    */
    if (!M || row < 1 || row > M->m || col < 1 || col > M->n) return 0;
    for (int k = 0; k < M->tu; k++) {
        if (M->elem[k].i == row && M->elem[k].j == col) {
            if (value == 0) {
                // Remove the element
                for (int j = k; j < M->tu - 1; j++) {
                    M->elem[j] = M->elem[j + 1];
                }
                M->tu--;
            } else {
                M->elem[k].data = value; // Update the value
            }
            return 1;
        }
    }
    if (value != 0) {
        if (M->tu >= MaxNum) return 0; // No space to add new element
        M->elem[M->tu].i = row;
        M->elem[M->tu].j = col;
        M->elem[M->tu].data = value;
        M->tu++;
    }
    return 1;
}

int mat_to_dense(SpMat* M, double** dense) {
    /*Convert a sparse matrix to a dense matrix
    Input:
        M: pointer to the sparse matrix
        dense: 2D array to store the dense matrix (should be pre-allocated with size M->m x M->n)
    Output:
        Returns 1 on success, 0 on failure (null pointers)  
    */
    if (!M || !dense) return 0;
    for (int i = 0; i < M->m; i++) {
        for (int j = 0; j < M->n; j++) {
            dense[i][j] = 0.0; // Initialize to zero
        }
    }
    for (int k = 0; k < M->tu; k++) {
        int row = M->elem[k].i - 1; // Convert to 0-based index
        int col = M->elem[k].j - 1; // Convert to 0-based index
        dense[row][col] = M->elem[k].data;
    }
    return 1;
}

int dense_to_mat(double** dense, int m, int n, SpMat* M) {
    if (!dense || !M) return 0;
    if (!CreateSpMat(M, m, n)) return 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (dense[i][j] != 0.0) {
                if (M->tu >= MaxNum) return 0; // No space to add new element
                M->elem[M->tu].i = i + 1; // Convert to 1-based index
                M->elem[M->tu].j = j + 1; // Convert to 1-based index
                M->elem[M->tu].data = dense[i][j];
                M->tu++;
            }
        }
    }
    return 1;
}

int mat_rmzero(SpMat* M) {
    /*remove zero elements from the sparse matrix
    Input:
        M: pointer to the sparse matrix
    Output:
        Returns 1 on success, 0 on failure (null pointer)
    */
    if (!M) return 0;
    int new_tu = 0;
    for (int k = 0; k < M->tu; k++) {
        if (M->elem[k].data != 0) {
            if (new_tu != k) {
                M->elem[new_tu] = M->elem[k];
            }
            new_tu++;
        }
    }
    M->tu = new_tu;
    return 1;
}

int mat_copy(SpMat* src, SpMat* dest) {
    /*copy a sparse matrix from src to dest
    Input:
        src: pointer to the source sparse matrix
        dest: pointer to the destination sparse matrix
    Output:
        Returns 1 on success, 0 on failure (null pointers or memory allocation failure)
    */
    if (!src || !dest) return 0;
    if (!CreateSpMat(dest, src->m, src->n)) return 0;
    for (int k = 0; k < src->tu; k++) {
        if (dest->tu >= MaxNum) return 0; // No space to add new element
        dest->elem[dest->tu] = src->elem[k];
        dest->tu++;
    }
    return 1;
}

int mat_print(SpMat* M) {
    /*print the sparse matrix to standard output
    Input:
        M: pointer to the sparse matrix
    Output:
        Returns 1 on success, 0 on failure (null pointer)
    */
    if (!M) return 0;
    printf("Sparse Matrix (%d x %d) with %d non-zero elements:\n", M->m, M->n, M->tu);
    for (int k = 0; k < M->tu; k++) {
        printf("Row: %d, Column: %d, Data: %.2f\n", M->elem[k].i, M->elem[k].j, M->elem[k].data);
    }
    return 1;
}

int mat_LU_decompose(SpMat* A, SpMat* L, SpMat* U) {
    /*decompose a square sparse matrix A into lower triangular matrix L and upper triangular matrix U using LU decomposition
    Input:
        A: pointer to the square sparse matrix to be decomposed
        L: pointer to the resulting lower triangular sparse matrix
        U: pointer to the resulting upper triangular sparse matrix
    Output:
        Returns 1 on success, 0 on failure (non-square matrix or memory allocation failure)
    */
    if (!A || !L || !U || A->m != A->n) return 0; // Only square matrices can be decomposed
    int n = A->m;
    if (!CreateSpMat(L, n, n) || !CreateSpMat(U, n, n)) return 0;

    // Parallel LU decomposition for larger matrices
    if (n > 50) {
        for (int i = 1; i <= n; i++) {
            // Upper triangular part - parallelize within each row
            #pragma omp parallel for schedule(dynamic) if(n-i > 10)
            for (int j = i; j <= n; j++) {
                double sum = 0.0;
                for (int k = 1; k < i; k++) {
                    double l_ik = 0, u_kj = 0;
                    
                    // Find L[i,k]
                    for (int p = 0; p < L->tu; p++) {
                        if (L->elem[p].i == i && L->elem[p].j == k) {
                            l_ik = L->elem[p].data;
                            break;
                        }
                    }
                    
                    // Find U[k,j]
                    for (int q = 0; q < U->tu; q++) {
                        if (U->elem[q].i == k && U->elem[q].j == j) {
                            u_kj = U->elem[q].data;
                            break;
                        }
                    }
                    sum += l_ik * u_kj;
                }
                
                double a_ij = 0;
                for (int p = 0; p < A->tu; p++) {
                    if (A->elem[p].i == i && A->elem[p].j == j) {
                        a_ij = A->elem[p].data;
                        break;
                    }
                }
                
                double u_ij = a_ij - sum;
                if (u_ij != 0) {
                    #pragma omp critical
                    {
                        if (U->tu < MaxNum) {
                            U->elem[U->tu].i = i;
                            U->elem[U->tu].j = j;
                            U->elem[U->tu].data = u_ij;
                            U->tu++;
                        }
                    }
                }
            }

            // Lower triangular part
            #pragma omp parallel for schedule(dynamic) if(n-i > 10)
            for (int j = i + 1; j <= n; j++) {
                double sum = 0.0;
                for (int k = 1; k < i; k++) {
                    double l_jk = 0, u_ki = 0;
                    
                    for (int p = 0; p < L->tu; p++) {
                        if (L->elem[p].i == j && L->elem[p].j == k) {
                            l_jk = L->elem[p].data;
                            break;
                        }
                    }
                    
                    for (int q = 0; q < U->tu; q++) {
                        if (U->elem[q].i == k && U->elem[q].j == i) {
                            u_ki = U->elem[q].data;
                            break;
                        }
                    }
                    sum += l_jk * u_ki;
                }
                
                double a_ji = 0;
                for (int p = 0; p < A->tu; p++) {
                    if (A->elem[p].i == j && A->elem[p].j == i) {
                        a_ji = A->elem[p].data;
                        break;
                    }
                }
                
                // Find U(i,i)
                double u_ii = 0;
                for (int q = 0; q < U->tu; q++) {
                    if (U->elem[q].i == i && U->elem[q].j == i) {
                        u_ii = U->elem[q].data;
                        break;
                    }
                }
                
                if (u_ii != 0) {
                    double l_ji = (a_ji - sum) / u_ii;
                    if (l_ji != 0) {
                        #pragma omp critical
                        {
                            if (L->tu < MaxNum) {
                                L->elem[L->tu].i = j;
                                L->elem[L->tu].j = i;
                                L->elem[L->tu].data = l_ji;
                                L->tu++;
                            }
                        }
                    }
                }
            }

            // Set diagonal of L to 1
            if (L->tu < MaxNum) {
                L->elem[L->tu].i = i;
                L->elem[L->tu].j = i;
                L->elem[L->tu].data = 1.0;
                L->tu++;
            }
        }
    } else {
        // Keep your existing serial implementation for small matrices
        for (int i = 1; i <= n; i++) {
            // Upper Triangular
            for (int j = i; j <= n; j++) {
                double sum = 0.0;
                for (int k = 1; k < i; k++) {
                    double l_ik = 0, u_kj = 0;
                    for (int p = 0; p < L->tu; p++) {
                        if (L->elem[p].i == i && L->elem[p].j == k) {
                            l_ik = L->elem[p].data;
                            break;
                        }
                    }
                    for (int q = 0; q < U->tu; q++) {
                        if (U->elem[q].i == k && U->elem[q].j == j) {
                            u_kj = U->elem[q].data;
                            break;
                        }
                    }
                    sum += l_ik * u_kj;
                }
                double a_ij = 0;
                for (int p = 0; p < A->tu; p++) {
                    if (A->elem[p].i == i && A->elem[p].j == j) {
                        a_ij = A->elem[p].data;
                        break;
                    }
                }
                double u_ij = a_ij - sum;
                if (u_ij != 0) {
                    if (U->tu >= MaxNum) return 0;
                    U->elem[U->tu].i = i;
                    U->elem[U->tu].j = j;
                    U->elem[U->tu].data = u_ij;
                    U->tu++;
                }
            }

            // Lower Triangular
            for (int j = i + 1; j <= n; j++) {
                double sum = 0.0;
                for (int k = 1; k < i; k++) {
                    double l_jk = 0, u_ki = 0;
                    for (int p = 0; p < L->tu; p++) {
                        if (L->elem[p].i == j && L->elem[p].j == k) {
                            l_jk = L->elem[p].data;
                            break;
                        }
                    }
                    for (int q = 0; q < U->tu; q++) {
                        if (U->elem[q].i == k && U->elem[q].j == i) {
                            u_ki = U->elem[q].data;
                            break;
                        }
                    }
                    sum += l_jk * u_ki;
                }
                double a_ji = 0;
                for (int p = 0; p < A->tu; p++) {
                    if (A->elem[p].i == j && A->elem[p].j == i) {
                        a_ji = A->elem[p].data;
                        break;
                    }
                }
                double u_ii = 0;
                for (int q = 0; q < U->tu; q++) {
                    if (U->elem[q].i == i && U->elem[q].j == i) {
                        u_ii = U->elem[q].data;
                        break;
                    }
                }
                if (u_ii == 0) return 0;
                double l_ji = (a_ji - sum) / u_ii;
                if (l_ji != 0) {
                    if (L->tu >= MaxNum) return 0;
                    L->elem[L->tu].i = j;
                    L->elem[L->tu].j = i;
                    L->elem[L->tu].data = l_ji;
                    L->tu++;
                }
            }
            // Set diagonal of L to 1
            if (L->tu >= MaxNum) return 0;
            L->elem[L->tu].i = i;
            L->elem[L->tu].j = i;
            L->elem[L->tu].data = 1.0;
            L->tu++;
        }
    }
    
    mat_rmzero(L);
    mat_rmzero(U);
    return 1; // Decomposition successful
}

int mat_upper_tri(SpMat* A, SpMat* U){
    if (!A || !U || A->m != A->n) return 0; // Only square matrices can be processed
    int n = A->m;
    if (!CreateSpMat(U, n, n)) return 0;
    for (int k = 0; k < A->tu; k++) {
        if (A->elem[k].i < A->elem[k].j) {
            U->elem[U->tu++] = A->elem[k];
        }
    }
    return 1;
}

int mat_lower_tri(SpMat* A, SpMat* L){
    if (!A || !L || A->m != A->n) return 0; // Only square matrices can be processed
    int n = A->m;
    if (!CreateSpMat(L, n, n)) return 0;
    for (int k = 0; k < A->tu; k++) {
        if (A->elem[k].i > A->elem[k].j) {
            L->elem[L->tu++] = A->elem[k];
        }
    }
    return 1;
}

int mat_diag(SpMat* A, SpMat* D) {
    if (!A || !D || A->m != A->n) return 0; // Only square matrices can have a diagonal matrix
    int n = A->m;
    if (!CreateSpMat(D, n, n)) return 0;
    for (int k = 0; k < A->tu; k++) {
        if (A->elem[k].i == A->elem[k].j) {
            D->elem[D->tu++] = A->elem[k];
        }
    }
    return 1;
}

// Precompute row offsets and index list (CSR) and values for fast per-row access.
// Allocates *row_start (size n+1), *row_idx (size nnz), *row_val (size nnz).
// Caller must free() row_start, row_idx, row_val.
static int precompute_row_index(SpMat* A, int** row_start_out, int** row_idx_out, double** row_val_out) {
    if (!A || !row_start_out || !row_idx_out || !row_val_out) return 0;
    int n = A->m;
    int nnz = A->tu;
    int *counts = (int*)calloc(n, sizeof(int));
    if (!counts) return 0;

    // count entries per row
    for (int k = 0; k < nnz; ++k) {
        int r = A->elem[k].i - 1;
        if (r >= 0 && r < n) ++counts[r];
    }

    int *row_start = (int*)malloc((n + 1) * sizeof(int));
    if (!row_start) { free(counts); return 0; }
    row_start[0] = 0;
    for (int i = 0; i < n; i++) row_start[i + 1] = row_start[i] + counts[i];

    int *row_idx = (int*)malloc(nnz * sizeof(int));
    double *row_val = (double*)malloc(nnz * sizeof(double));
    if (!row_idx || !row_val) {
        free(counts);
        free(row_start);
        if (row_idx) free(row_idx);
        if (row_val) free(row_val);
        *row_start_out = NULL;
        *row_idx_out = NULL;
        *row_val_out = NULL;
        return 0;
    }

    // fill indices and values by scanning elements and placing them into CSR arrays
    for (int i = 0; i < n; ++i) counts[i] = 0;
    for (int k = 0; k < nnz; ++k) {
        int r = A->elem[k].i - 1;
        if (r >= 0 && r < n) {
            int pos = row_start[r] + counts[r];
            row_idx[pos] = A->elem[k].j - 1;   // store 0-based column index
            row_val[pos] = A->elem[k].data;
            ++counts[r];
        }
    }

    free(counts);
    *row_start_out = row_start;
    *row_idx_out = row_idx;
    *row_val_out = row_val;
    return 1;
}

// Comparator for qsort: sort by row then column (1-based stored in elem)
static int cmp_triple_rowcol(const void* a, const void* b) {
    const triple* A = (const triple*)a;
    const triple* B = (const triple*)b;
    if (A->i < B->i) return -1;
    if (A->i > B->i) return 1;
    if (A->j < B->j) return -1;
    if (A->j > B->j) return 1;
    return 0;
}

/* Gauss-Seidel method for solving Ax = b
   Implementation: sort elements by row once, then for each sweep loop A->elem
   sequentially (single pass over all nonzeros produces all row updates).
   This avoids per-row index arrays while keeping in-place GS semantics. */
int Gauss_Seidel(SpMat* A, double* b, double* x, int max_iter, double tol) {
    if (!A || !b || !x || A->m != A->n) return 0;
    if (A->tu <= 0) return 0;
    int n = A->m;

    // sort nonzeros by row then column so a single linear scan visits rows in order
    qsort(A->elem, (size_t)A->tu, sizeof(triple), cmp_triple_rowcol);

    const double diag_eps = 1e-12;
    int warned_zero_diag = 0;

    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        double diff_norm = 0.0;

        int k = 0;
        // single pass over all nonzeros (rows are contiguous after qsort)
        while (k < A->tu) {
            int row = A->elem[k].i - 1;
            double diag = 0.0;
            double sum = 0.0;

            // process all elements in this row
            while (k < A->tu && (A->elem[k].i - 1) == row) {
                int col = A->elem[k].j - 1;
                double v = A->elem[k].data;
                if (col == row) {
                    diag = v;
                } else {
                    // use current x[col] (in-place GS semantics)
                    sum += v * x[col];
                }
                ++k;
            }

            if (diag == 0.0) {
                diag = diag_eps;
                warned_zero_diag = 1;
            }

            double x_old = x[row];
            x[row] = (b[row] - sum) / diag;
            double d = x[row] - x_old;
            diff_norm += d * d;
        }

        diff_norm = sqrt(diff_norm);
        if (diff_norm < tol) break;
    }
    
    if (warned_zero_diag) {
        fprintf(stderr, "Warning: Gauss_Seidel replaced missing diagonal entries with eps=%g\n", diag_eps);
    }
    // print number of iterations taken
    printf("Gauss-Seidel converged in %d iterations\n", iter);
    return 1;
}

// Jacobi method for solving Ax = b
int Jacobi(SpMat* A, double* b, double* x, int max_iter, double tol) {
    if (!A || !b || !x || A->m != A->n) return 0;
    int n = A->m;

    double* x_new = (double*)malloc(n * sizeof(double));
    if (!x_new) return 0;

    /* pre-store diagonal elements (one-time) */
    double* diag = (double*)malloc(n * sizeof(double));
    if (!diag) { free(x_new); return 0; }
    for (int i = 0; i < n; ++i) diag[i] = 0.0;
    for (int k = 0; k < A->tu; ++k) {
        if (A->elem[k].i == A->elem[k].j) {
            int r = A->elem[k].i - 1;
            if (r >= 0 && r < n) diag[r] = A->elem[k].data;
        }
    }
    const double diag_eps = 1e-12;
    int zero_diag_count = 0;
    for (int i = 0; i < n; ++i) {
        if (diag[i] == 0.0) { diag[i] = diag_eps; ++zero_diag_count; }
    }
    if (zero_diag_count) {
        fprintf(stderr, "Warning: Jacobi found %d zero diagonal entries; using eps=%g\n",
                zero_diag_count, diag_eps);
    }

    /* workspace: row-wise accumulated off-diagonal contributions */
    double* row_sum = (double*)malloc(n * sizeof(double));
    if (!row_sum) { free(x_new); free(diag); return 0; }

    for (int iter = 0; iter < max_iter; ++iter) {
        /* zero row_sum fast */
        memset(row_sum, 0, n * sizeof(double));

        /* single pass over nonzeros to accumulate row sums:
           row_sum[i] += A(i,j) * x_old[j] for j != i.
           Use OpenMP atomics to allow parallel element-wise accumulation. */
        #pragma omp parallel for schedule(static)
        for (int kk = 0; kk < A->tu; ++kk) {
            int i = A->elem[kk].i - 1;
            int j = A->elem[kk].j - 1;
            if (i < 0 || i >= n || j < 0 || j >= n) continue;
            if (i == j) continue;
            double v = A->elem[kk].data;
            double term = v * x[j];
            #pragma omp atomic
            row_sum[i] += term;
        }

        /* form new iterate and compute difference norm (parallel) */
        double diff_norm = 0.0;
        #pragma omp parallel for reduction(+:diff_norm) schedule(static)
        for (int i = 0; i < n; ++i) {
            x_new[i] = (b[i] - row_sum[i]) / diag[i];
            double d = x_new[i] - x[i];
            diff_norm += d * d;
        }

        /* copy new values back to x (parallel) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) x[i] = x_new[i];

        diff_norm = sqrt(diff_norm);
        if (diff_norm < tol) break;
    }

    free(row_sum);
    free(diag);
    free(x_new);
    return 1;
}

