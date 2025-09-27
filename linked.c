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
int addrnode(ElLink* head, int col, double data);
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



int main() {
    MatList ML;
    if (CreateSparse(&ML, 3, 4)) {
        printf("sparse list created successfully.\n");
        printf("Rows: %d, Columns: %d\n", ML.m, ML.n);
    } else {
        printf("Failed to create sparse list.\n");
    }

    addrnode(&ML.rhead[1], 2, 5.0);
    addrnode(&ML.rhead[1], 3, 8.0);
    addrnode(&ML.rhead[1], 2, -5.0); // This should remove the node at column 2

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
    // Free allocated memory (not shown here for brevity)
    free(ML.rhead);
    free(ML.chead);

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

int addrnode(ElLink* head, int col, double data) {
    ElLink newNode = (ElLink)malloc(sizeof(ElNode));
    if (!newNode) return 0;
    newNode->col = col;
    newNode->data = data;
    newNode->rnext = NULL;

    if (*head == NULL || (*head)->col > col) {
        newNode->rnext = *head;
        *head = newNode;
    }else if ((*head)->col == col)
    {
        (*head)->data += data;
        if ((*head)->data == 0) {
            ElLink temp = *head;
            *head = (*head)->rnext;
            free(temp);
        }
    }
     else {
        ElLink current = *head;
        while (current->rnext != NULL && current->rnext->col < col) {
            current = current->rnext;
        }// col rise
        if (current->rnext != NULL && current->rnext->col == col) {
            current->rnext->data += data;// if same col, add data
            // If the data becomes zero, remove the node
            if (current->rnext->data == 0) {
                ElLink temp = current->rnext;
                current->rnext = current->rnext->rnext;
                free(temp);
            }
            free(newNode);
        } else {
            newNode->rnext = current->rnext;
            current->rnext = newNode;
        }
    }
    return 1;
}

int sparse_append(MatList* ML, int* rows, int* cols, double* data, int count) {
    int success = 1;
    if (count > 1000) {
        // Use OpenMP parallel for large counts
        #pragma omp parallel for default(none) shared(ML, rows, cols, data, count) reduction(&&:success)
        for (int i = 0; i < count; i++) {
            int local_success = 1;
            if (rows[i] < 1 || rows[i] > ML->m || cols[i] < 1 || cols[i] > ML->n) {
                local_success = 0;
            } else if (data[i] != 0) {
                int ok1, ok2;
                #pragma omp critical(rhead)
                {
                    ok1 = addrnode(&ML->rhead[rows[i]], cols[i], data[i]);
                }
                #pragma omp critical(chead)
                {
                    ok2 = addrnode(&ML->chead[cols[i]], rows[i], data[i]);
                }
                if (!ok1 || !ok2) local_success = 0;
                #pragma omp atomic
                ML->total++;
            }
            if (!local_success) success = 0;
        }
    } else {
        // Use normal for loop for small counts
        for (int i = 0; i < count; i++) {
            int local_success = 1;
            if (rows[i] < 1 || rows[i] > ML->m || cols[i] < 1 || cols[i] > ML->n) {
                local_success = 0;
            } else if (data[i] != 0) {
                int ok1 = addrnode(&ML->rhead[rows[i]], cols[i], data[i]);
                int ok2 = addrnode(&ML->chead[cols[i]], rows[i], data[i]);
                if (!ok1 || !ok2) local_success = 0;
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
                addrnode(&C->rhead[i], pa->col, pa->data);
                pa = pa->rnext;
            } else if (pa->col > pb->col) {
                addrnode(&C->rhead[i], pb->col, pb->data);
                pb = pb->rnext;
            } else {
                double sum = pa->data + pb->data;
                if (sum != 0) {
                    addrnode(&C->rhead[i], pa->col, sum);
                }
                pa = pa->rnext;
                pb = pb->rnext;
            }
        }
        while (pa != NULL) {
            addrnode(&C->rhead[i], pa->col, pa->data);
            pa = pa->rnext;
        }
        while (pb != NULL) {
            addrnode(&C->rhead[i], pb->col, pb->data);
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

    for (int i = 1; i <= A->m; i++) {
        for (int j = 1; j <= B->n; j++) {
            double sum = 0.0;
            ElLink pa = A->rhead[i];
            ElLink pb = B->chead[j];
            while (pa != NULL && pb != NULL) {
                if (pa->col < pb->col) {
                    pa = pa->rnext;
                } else if (pa->col > pb->col) {
                    pb = pb->cnext;
                } else {
                    sum += pa->data * pb->data;
                    pa = pa->rnext;
                    pb = pb->cnext;
                }
            }
            if (sum != 0.0) {
                addrnode(&C->rhead[i], j, sum);
                C->total++;
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
            addrnode(&AT->rhead[current->col], i, current->data);
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
                addrnode(&M->rhead[i + 1], j + 1, dense[i][j]);
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
            addrnode(&dest->rhead[i], current->col, current->data);
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
    for (int i = 1; i <= A->m; i++) {
        for (int j = 1; j <= A->n; j++) {
            double sum = 0.0;
            ElLink pa = A->rhead[i];
            ElLink pb = U->chead[j];
            while (pa != NULL && pb != NULL) {
                if (pa->col < pb->col) {
                    pa = pa->rnext;
                } else if (pa->col > pb->col) {
                    pb = pb->cnext;
                } else {
                    sum += pa->data * pb->data;
                    pa = pa->rnext;
                    pb = pb->cnext;
                }
            }
            double a_ij;
            sparse_get(A, i, j, &a_ij);
            if (i == j) {
                addrnode(&U->rhead[i], j, a_ij - sum);
                addrnode(&L->rhead[i], j, 1.0);
            } else if (i < j) {
                addrnode(&U->rhead[i], j, a_ij - sum);
            } else {
                double u_jj;
                sparse_get(U, j, j, &u_jj);
                if (u_jj == 0) {
                    sparse_free(L);
                    sparse_free(U);
                    return 0; // Singular sparse
                }
                addrnode(&L->rhead[i], j, (a_ij - sum) / u_jj);
            }
        }
    }

    return 1;
}