#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

typedef struct ElNode{
    int row,col;
    double data;
    struct ElNode * rnext,*cnext;
}ElNode, *ElLink;

typedef struct
{
    ElLink* rhead, * chead; 
    int m, n, total;  
}MatList;

int CreateSparse(MatList* M,int m,int n);
int addrnode2(ElLink* rhead, ElLink* chead, int row, int col, double data);
int sparse_append(MatList* ML, int* rows, int* cols, double* data, int count);
int sparse_add(MatList* A, MatList* B, MatList* C);
int sparse_multiply(MatList* A, MatList* B, MatList* C);
int sparse_transpose(MatList* A, MatList* AT);
int sparse_free(MatList* M);
int sparse_print(MatList* M);
int sparse_get(MatList* M, int row, int col, double* value);
int sparse_set(MatList* M, int row, int col, double value);
int sparse_to_dense(MatList* M, double** dense);
int dense_to_sparse(double** dense, int m, int n, MatList* M);
int sparse_copy(MatList* src, MatList* dest);
int sparse_LU_decompose(MatList* A, MatList* L, MatList* U);
int sparse_rmzero(MatList* M);



int main() {
    MatList ML;
    if (CreateSparse(&ML, 5, 5)) {
        printf("sparse list created successfully.\n");
        printf("Rows: %d, Columns: %d\n", ML.m, ML.n);
    } else {
        printf("Failed to create sparse list.\n");
    }

    addrnode2(ML.rhead, ML.chead, 1, 2, 5.0);
    addrnode2(ML.rhead, ML.chead, 1, 3, 8.0);
    addrnode2(ML.rhead, ML.chead, 1, 2, -5.0); // This should remove the node at column 2

    printf("Row 1 after adding nodes:\n");
    ElLink current = ML.rhead[1]; 

    while (current != NULL) {
        printf("Column: %d, Data: %.2f\n", current->col, current->data);
        current = current->rnext;
    }

    sparse_append(&ML, (int[]){1,2,2,3}, (int[]){1,2,3,3}, (double[]){1.0,2.0,5.0,3.0}, 4);

    printf("Total non-zero elements after appending: %d\n", ML.total);
    printf("Row 2 after appending nodes:\n");
    current = ML.rhead[2];
    while (current != NULL) {
        printf("Column: %d, Data: %.2f\n", current->col, current->data);
        current = current->rnext;
    }
  

  
   

    sparse_free(&ML);
    if (CreateSparse(&ML, 10, 10)) {
        ML.total = 0; // Ensure total is initialized to 0
        // Create 1D first-order central difference matrix (off-diagonals: +1/-1, diagonal: 0)
        int rows[24], cols[24];
        double data[24];
        int idx = 0;
        for (int i = 1; i <= 10; ++i) {
            if (i > 1) { rows[idx] = i; cols[idx] = i-1; data[idx] = -1.0; idx++; } // left neighbor
            if (i < 10) { rows[idx] = i; cols[idx] = i+1; data[idx] = 1.0; idx++; }  // right neighbor
        }
        rows[idx] = 10; cols[idx] = 10; data[idx] = 1.0; idx++; 
        rows[idx] = 10; cols[idx] = 9; data[idx] = 1.0; idx++;
        rows[idx] = 1; cols[idx] = 1; data[idx] = -1.0; idx++;  
        

        if (sparse_append(&ML, rows, cols, data, idx)) {
            printf("Order one central difference matrix created successfully.\n");
            sparse_print(&ML);
        } else {
            printf("Failed to append data to the sparse matrix.\n");
        }
    } else {
        printf("Failed to create sparse matrix.\n");
    }

    MatList ML2;
    if (sparse_copy(&ML, &ML2)) {
        printf("Sparse matrix copied successfully.\n");
        sparse_print(&ML2);
    } else {
        printf("Failed to copy sparse matrix.\n");
    }
    // test addition and multiplication

    MatList ML_sum, ML_prod;
    if (sparse_add(&ML, &ML2, &ML_sum)) {
        printf("Sparse matrix addition successful.\n");
        sparse_print(&ML_sum);
    } else {
        printf("Failed to add sparse matrices.\n");
    }

    if (sparse_multiply(&ML, &ML2, &ML_prod)) {
        printf("Sparse matrix multiplication successful.\n");
        sparse_print(&ML_prod);
    } else {
        printf("Failed to multiply sparse matrices.\n");
    }
    // test transpose
    MatList ML_transpose;
    if (sparse_transpose(&ML, &ML_transpose)) {
        printf("Sparse matrix transpose successful.\n");
        sparse_print(&ML_transpose);
    } else {
        printf("Failed to transpose sparse matrix.\n");
    }
    // test LU decomposition
    MatList L, U;
    if (sparse_LU_decompose(&ML, &L, &U)) {
        printf("Sparse matrix LU decomposition successful.\n");
        printf("L matrix:\n");
        sparse_print(&L);
        printf("U matrix:\n");
        sparse_print(&U);
        sparse_free(&L);
        sparse_free(&U);
    } else {
        printf("Failed to perform LU decomposition on sparse matrix.\n");
    }
    sparse_free(&ML2);
    sparse_free(&ML_sum);
    sparse_free(&ML_prod);
    sparse_free(&ML_transpose);
    sparse_free(&ML);
    return 0;
}



int CreateSparse(MatList* ML, int m, int n) {
    int i;
    ML->m = m;
    ML->n = n;
    ML->total = 0;

    ML->rhead = (ElLink*)malloc((m + 1) * sizeof(ElLink));
    if (!ML->rhead) return 0;
    for (i = 1; i <= m; i++) ML->rhead[i] = NULL;

    ML->chead = (ElLink*)malloc((n + 1) * sizeof(ElLink));
    if (!ML->chead) {
        free(ML->rhead);
        return 0;
    }
    for (i = 1; i <= n; i++) ML->chead[i] = NULL;

    return 1;
}


// Insert a node into both row and column lists
int addrnode2(ElLink* rhead, ElLink* chead, int row, int col, double data) {
    // Insert into row list (rnext)
    ElLink *rowp = &rhead[row];
    ElLink *colp = &chead[col];
    ElLink prow = *rowp, pcol = *colp;
    ElLink prev_row = NULL, prev_col = NULL;
    while (prow && prow->col < col) { prev_row = prow; prow = prow->rnext; }
    while (pcol && pcol->row < row) { prev_col = pcol; pcol = pcol->cnext; }
    // If already exists, update value
    if (prow && prow->col == col) {
        prow->data += data;
        if (prow->data == 0) {
            // Remove from row list
            if (prev_row) prev_row->rnext = prow->rnext; else *rowp = prow->rnext;
            // Remove from col list
            if (prev_col) prev_col->cnext = prow->cnext; else *colp = prow->cnext;
            free(prow);
        }
        return 1;
    }
    // Create new node
    ElLink newNode = (ElLink)malloc(sizeof(ElNode));
    if (!newNode) return 0;
    newNode->row = row;
    newNode->col = col;
    newNode->data = data;
    // Insert into row list
    newNode->rnext = prow;
    if (prev_row) prev_row->rnext = newNode; else *rowp = newNode;
    // Insert into col list
    newNode->cnext = pcol;
    if (prev_col) prev_col->cnext = newNode; else *colp = newNode;
    return 1;
}

int sparse_append(MatList* ML, int* rows, int* cols, double* data, int count) {
    int success = 1;
    if (count > 1000) {
        // Use OpenMP parallel for large counts
        #pragma omp parallel for default(none) shared(ML, rows, cols, data, count) reduction(&&:success)
        for (int i = 0; i < count; i++) {
            int local_success = 1;
            if (rows[i] == 0 || cols[i] == 0) continue;
            if (rows[i] < 1 || rows[i] > ML->m || cols[i] < 1 || cols[i] > ML->n) {
                local_success = 0;
            } else if (data[i] != 0) {
                int ok;
                #pragma omp critical(rchead)
                {
                    ok = addrnode2(ML->rhead, ML->chead, rows[i], cols[i], data[i]);
                }
                if (!ok) local_success = 0;
                #pragma omp atomic
                ML->total++;
            }
            if (!local_success) success = 0;
        }
    } else {
        // Use normal for loop for small counts
        for (int i = 0; i < count; i++) {
            int local_success = 1;
            if (rows[i] == 0 || cols[i] == 0) continue;
            if (rows[i] < 1 || rows[i] > ML->m || cols[i] < 1 || cols[i] > ML->n) {
                local_success = 0;
            } else if (data[i] != 0) {
                int ok = addrnode2(ML->rhead, ML->chead, rows[i], cols[i], data[i]);
                if (!ok) local_success = 0;
                ML->total++;
            }
            if (!local_success) success = 0;
        }
    }
    return success;
}

int sparse_add(MatList* A, MatList* B, MatList* C) {
    if (A->m != B->m || A->n != B->n) return 0; // Dimension mismatch
    if (!CreateSparse(C, A->m, A->n)) return 0;

    for (int i = 1; i <= A->m; i++) {
        ElLink pa = A->rhead[i];
        ElLink pb = B->rhead[i];
        while (pa != NULL && pb != NULL) {
            if (pa->col < pb->col) {
                addrnode2(C->rhead, C->chead, i, pa->col, pa->data);
                pa = pa->rnext;
            } else if (pa->col > pb->col) {
                addrnode2(C->rhead, C->chead, i, pb->col, pb->data);
                pb = pb->rnext;
            } else {
                double sum = pa->data + pb->data;
                if (sum != 0) {
                    addrnode2(C->rhead, C->chead, i, pa->col, sum);
                }
                pa = pa->rnext;
                pb = pb->rnext;
            }
        }
        while (pa != NULL) {
            addrnode2(C->rhead, C->chead, i, pa->col, pa->data);
            pa = pa->rnext;
        }
        while (pb != NULL) {
            addrnode2(C->rhead, C->chead, i, pb->col, pb->data);
            pb = pb->rnext;
        }
        
    }

    // Update total count in C
    for (int i = 1; i <= C->m; i++) {
        ElLink current = C->rhead[i];
        while (current != NULL) {
            C->total++;
            current = current->rnext;
        }
    }

    return 1;
}

int sparse_multiply(MatList* A, MatList* B, MatList* C) {
    if (A->n != B->m) return 0; // Dimension mismatch
    if (!CreateSparse(C, A->m, B->n)) return 0;

    int m = A->m, n = B->n;
    // Parallelize only if the problem size is large
    if (m * n > 1000) {
        #pragma omp parallel for default(none) shared(A, B, C, m, n)
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                double sum = 0.0;
                ElLink pa = A->rhead[i];
                ElLink pb = B->chead[j];
                while (pa != NULL && pb != NULL) {
                    if (pa->col < pb->row) {
                        pa = pa->rnext;
                    } else if (pa->col > pb->row) {
                        pb = pb->cnext;
                    } else {
                        sum += pa->data * pb->data;
                        pa = pa->rnext;
                        pb = pb->cnext;
                    }
                }
                if (sum != 0.0) {
                    #pragma omp critical(sparse_mult_addr)
                    {
                        addrnode2(C->rhead, C->chead, i, j, sum);
                        C->total++;
                    }
                }
            }
        }
    } else {
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                double sum = 0.0;
                ElLink pa = A->rhead[i];
                ElLink pb = B->chead[j];
                //printf("%d %d\n", pa != NULL, pb != NULL);
                while (pa != NULL && pb != NULL) {
                    if (pa->col < pb->row) {
                        pa = pa->rnext;
                    } else if (pa->col > pb->row) {
                        pb = pb->cnext;
                    } else {
                        sum += pa->data * pb->data;
                        pa = pa->rnext;
                        pb = pb->cnext;
                    }
                }
                if (sum != 0.0) {
                    addrnode2(C->rhead, C->chead, i, j, sum);
                    C->total++;
                }
            }
        }
    }

    return 1;
}

int sparse_transpose(MatList* A, MatList* AT) {
    if (!CreateSparse(AT, A->n, A->m)) return 0;

    for (int i = 1; i <= A->m; i++) {
        ElLink current = A->rhead[i];
        while (current != NULL) {
            addrnode2(AT->rhead, AT->chead, current->col, i, current->data);
            AT->total++;
            current = current->rnext;
        }
    }

    return 1;
}

int sparse_free(MatList* M) {
    if (!M) return 0;
    for (int i = 1; i <= M->m; i++) {
        ElLink current = M->rhead[i];
        while (current != NULL) {
            ElLink temp = current;
            current = current->rnext;
            free(temp);
        }
    }
    free(M->rhead);
    free(M->chead);
    M->rhead = NULL;
    M->chead = NULL;
    M->m = 0;
    M->n = 0;
    M->total = 0;
    return 1;
}


int sparse_print(MatList* M) {
    if (!M) return 0;
    for (int i = 1; i <= M->m; i++) {
        ElLink current = M->rhead[i];
        printf("Row %d: ", i);
        while (current != NULL) {
            printf("(Col: %d, Data: %.2f) ", current->col, current->data);
            current = current->rnext;
        }
        printf("\n");
    }
    return 1;
}

int sparse_get(MatList* M, int row, int col, double* value) {
    if (row < 1 || row > M->m || col < 1 || col > M->n) return 0;
    ElLink current = M->rhead[row];
    while (current != NULL) {
        if (current->col == col) {
            *value = current->data;
            return 1;
        } else if (current->col > col) {
            break;
        }
        current = current->rnext;
    }
    *value = 0.0; // Element not found, return 0
    return 1;
}

int sparse_set(MatList* M, int row, int col, double value) {
    if (row < 1 || row > M->m || col < 1 || col > M->n) return 0;
    ElLink* head = &M->rhead[row];
    ElLink current = *head;
    ElLink prev = NULL;

    while (current != NULL && current->col < col) {
        prev = current;
        current = current->rnext;
    }

    if (current != NULL && current->col == col) {
        if (value == 0) {
            // Remove the node
            if (prev == NULL) {
                *head = current->rnext;
            } else {
                prev->rnext = current->rnext;
            }
            free(current);
            M->total--;
        } else {
            current->data = value; // Update the value
        }
    } else {
        if (value != 0) {
            // Insert new node
            ElLink newNode = (ElLink)malloc(sizeof(ElNode));
            if (!newNode) return 0;
            newNode->col = col;
            newNode->data = value;
            newNode->rnext = current;

            if (prev == NULL) {
                *head = newNode;
            } else {
                prev->rnext = newNode;
            }
            M->total++;
        }
    }
    return 1;
}

int sparse_to_dense(MatList* M, double** dense) {
    if (!M || !dense) return 0;
    for (int i = 0; i < M->m; i++) {
        for (int j = 0; j < M->n; j++) {
            dense[i][j] = 0.0;
        }
    }
    for (int i = 1; i <= M->m; i++) {
        ElLink current = M->rhead[i];
        while (current != NULL) {
            dense[i - 1][current->col - 1] = current->data;
            current = current->rnext;
        }
    }
    return 1;
}

int dense_to_sparse(double** dense, int m, int n, MatList* M) {
    if (!dense || !M) return 0;
    if (!CreateSparse(M, m, n)) return 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (dense[i][j] != 0.0) {
                addrnode2(M->rhead, M->chead, i + 1, j + 1, dense[i][j]);
                M->total++;
            }
        }
    }
    return 1;
}

int sparse_copy(MatList* src, MatList* dest) {
    if (!src || !dest) return 0;
    if (!CreateSparse(dest, src->m, src->n)) return 0;
    for (int i = 1; i <= src->m; i++) {
        ElLink current = src->rhead[i];
        while (current != NULL) {
            addrnode2(dest->rhead, dest->chead, i, current->col, current->data);
            dest->total++;
            current = current->rnext;
        }
    }
    return 1;
}

int sparse_LU_decompose(MatList* A, MatList* L, MatList* U) {
    if (A->m != A->n) return 0; // LU decomposition requires square sparse
    if (!CreateSparse(L, A->m, A->n)) return 0;
    if (!CreateSparse(U, A->m, A->n)) {
        sparse_free(L);
        return 0;
    }
    // Doolittle LU: U upper, L lower with 1 on diagonal
    for (int i = 1; i <= A->m; i++) {
        // Compute U(i, j) for j >= i
        for (int j = i; j <= A->n; j++) {
            double sum = 0.0;
            for (int k = 1; k < i; k++) {
                double l_ik, u_kj;
                sparse_get(L, i, k, &l_ik);
                sparse_get(U, k, j, &u_kj);
                sum += l_ik * u_kj;
            }
            double a_ij;
            sparse_get(A, i, j, &a_ij);
            addrnode2(U->rhead, U->chead, i, j, a_ij - sum);
        }
        // Compute L(j, i) for j > i
        for (int j = i + 1; j <= A->m; j++) {
            double sum = 0.0;
            for (int k = 1; k < i; k++) {
                double l_jk, u_ki;
                sparse_get(L, j, k, &l_jk);
                sparse_get(U, k, i, &u_ki);
                sum += l_jk * u_ki;
            }
            double a_ji, u_ii;
            sparse_get(A, j, i, &a_ji);
            sparse_get(U, i, i, &u_ii);
            if (u_ii == 0) {
                sparse_free(L);
                sparse_free(U);
                return 0; // Singular sparse
            }
            addrnode2(L->rhead, L->chead, j, i, (a_ji - sum) / u_ii);
        }
        // Set diagonal of L to 1
        addrnode2(L->rhead, L->chead, i, i, 1.0);
    }
    // Remove explicit zeros from L and U
    sparse_rmzero(L);
    sparse_rmzero(U);
    return 1;
}

int sparse_rmzero(MatList* M) {
    if (!M) return 0;
    for (int i = 1; i <= M->m; i++) {
        ElLink current = M->rhead[i];
        ElLink prev = NULL;
        while (current != NULL) {
            if (current->data == 0.0) {
                if (prev == NULL) {
                    M->rhead[i] = current->rnext;
                } else {
                    prev->rnext = current->rnext;
                }
                ElLink temp = current;
                current = current->rnext;
                free(temp);
                M->total--;
            } else {
                prev = current;
                current = current->rnext;
            }
        }
    }
    return 1;
}

// take the upper triangular part of A as U

int upper_triangular(MatList* A, MatList* U) {
    if (A->m != A->n) return 0; // LU decomposition requires square sparse
    if (!CreateSparse(U, A->m, A->n)) return 0;
    for (int i = 1; i <= A->m; i++) {
        for (int j = 1; j <= A->n; j++) {
            double a_ij;
            sparse_get(A, i, j, &a_ij);
            if (i <= j) {
                addrnode2(U->rhead, U->chead, i, j, a_ij);
            }
        }
    }
    return 1;
}