#include<stdio.h>
#define MaxNum 1000
//三元组
typedef struct {
    int i;
    int j;
    double data;
}triple, *MatHead;
//三元组顺序表
typedef struct {
    triple elem[MaxNum];
    int m;
    int n;
    int tu;
}SpMat;

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

int main() {
    SpMat M, T;
    if (CreateSpMat(&M, 3, 3)) {
        printf("Sparse matrix created successfully.\n");
        printf("Rows: %d, Columns: %d\n", M.m, M.n);
    } else {
        printf("Failed to create sparse matrix.\n");
    }
   
    //test append
    int rows[] = {1, 1, 2, 3, 3};
    int cols[] = {1, 3, 2, 1, 3};
    double data[] = {5.0, 8.0, 3.0, 6.0, 9.0};
    if (mat_append(&M, rows, cols, data, 5)) {
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

    SpMat A, B, C;
    CreateSpMat(&A, 3, 3);
    CreateSpMat(&B, 3, 3);
    int rowsA[] = {1, 2, 3};
    int colsA[] = {1, 2, 3};
    double dataA[] = {1.0, 2.0, 3.0};
    mat_append(&A, rowsA, colsA, dataA, 3);
    int rowsB[] = {1, 2, 3};
    int colsB[] = {2, 3, 1};
    double dataB[] = {4.0, 5.0, 6.0};
    mat_append(&B, rowsB, colsB, dataB, 3);
    if (mat_add(&A, &B, &C)) {
        printf("Matrices added successfully.\n");
        printf("Resultant matrix non-zero elements:\n");
        for (int k = 0; k < C.tu; k++) {
            printf("Row: %d, Column: %d, Data: %.2f\n", C.elem[k].i, C.elem[k].j, C.elem[k].data);
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
    return 0;
}



int CreateSpMat(SpMat* M, int m, int n) {
    if (m < 1 || n < 1 || m > MaxNum || n > MaxNum) return 0;
    M->m = m;
    M->n = n;
    M->tu = 0;
    return 1;
}

int mat_append(SpMat* M, int rows[], int cols[], double data[], int count) {
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
    if (A->m != B->m || A->n != B->n) return 0; //维数不同，无法相加
    C->m = A->m;
    C->n = A->n;
    C->tu = 0;
    int i = 0, j = 0, k = 0;
    if (A->tu + B->tu > 1000) {
        // Parallelize the merge for large input (not fully optimal, but illustrative)
        #pragma omp parallel for default(none) shared(A, B, C) firstprivate(i, j, k)
        for (int idx = 0; idx < A->tu + B->tu; idx++) {
            // Not a true parallel merge, but for illustration only
        }
        // Fallback to serial for correctness
    }
    while (i < A->tu && j < B->tu) {
        if (A->elem[i].i < B->elem[j].i || (A->elem[i].i == B->elem[j].i && A->elem[i].j < B->elem[j].j)) {
            C->elem[k++] = A->elem[i++];
            C->tu++;
        } else if (A->elem[i].i > B->elem[j].i || (A->elem[i].i == B->elem[j].i && A->elem[i].j > B->elem[j].j)) {
            C->elem[k++] = B->elem[j++];
            C->tu++;
        } else {
            int sum = A->elem[i].data + B->elem[j].data;
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
    while (i < A->tu) {
        C->elem[k++] = A->elem[i++];
        C->tu++;
    }
    while (j < B->tu) {
        C->elem[k++] = B->elem[j++];
        C->tu++;
    }
    return 1;
}

int mat_multiply(SpMat* A, SpMat* B, SpMat* C) {
    if (A->n != B->m) return 0; //维数不匹配，无法相乘
    C->m = A->m;
    C->n = B->n;
    C->tu = 0;
    int m = A->m, n = B->n;
    if (m * n > 1000) {
        #pragma omp parallel for default(none) shared(A, B, C, m, n)
        for (int idx = 0; idx < m * n; idx++) {
            int i = idx / n + 1;
            int j = idx % n + 1;
            int sum = 0;
            for (int k = 1; k <= A->n; k++) {
                int a_ik = 0, b_kj = 0;
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
                int idxc;
                #pragma omp critical
                {
                    if (C->tu >= MaxNum) return 0;
                    idxc = C->tu++;
                }
                C->elem[idxc].i = i;
                C->elem[idxc].j = j;
                C->elem[idxc].data = sum;
            }
        }
        return 1;
    } else {
        for (int i = 1; i <= A->m; i++) {
            for (int j = 1; j <= B->n; j++) {
                int sum = 0;
                for (int k = 1; k <= A->n; k++) {
                    int a_ik = 0, b_kj = 0;
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
                    if (C->tu >= MaxNum) return 0; // 超出存储空间
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

//稀疏矩阵的转置
int mat_transpose(SpMat* M, SpMat* T) {
    //1.稀疏矩阵的行数和列数互换
    (*T).m = M->n;
    (*T).n = M->m;
    (*T).tu = M->tu;

    if ((*T).tu) {
        int col, p;
        int q = 0;
        //2.遍历原表中的各个三元组
        for (col = 1; col <= M->n; col++) {
            //重复遍历原表 M.m 次，将所有三元组都存放到新表中
            for (p = 0; p < M->tu; p++) {
                //3.每次遍历，找到 j 列值最小的三元组，将它的行、列互换后存储到新表中
                if (M->elem[p].j == col) {
                    (*T).elem[q].i = M->elem[p].j;
                    (*T).elem[q].j = M->elem[p].i;
                    (*T).elem[q].data = M->elem[p].data;
                    //为存放下一个三元组做准备
                    q++;
                }
            }
        }
    }
    return 1;
}


int mat_get(SpMat* M, int row, int col, double* value) {
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
    if (!M) return 0;
    printf("Sparse Matrix (%d x %d) with %d non-zero elements:\n", M->m, M->n, M->tu);
    for (int k = 0; k < M->tu; k++) {
        printf("Row: %d, Column: %d, Data: %.2f\n", M->elem[k].i, M->elem[k].j, M->elem[k].data);
    }
    return 1;
}

int mat_LU_decompose(SpMat* A, SpMat* L, SpMat* U) {
    if (!A || !L || !U || A->m != A->n) return 0; // Only square matrices can be decomposed
    int n = A->m;
    if (!CreateSpMat(L, n, n) || !CreateSpMat(U, n, n)) return 0;

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
                if (U->tu >= MaxNum) return 0; // No space to add new element
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
                if (U->elem[i - 1].data == 0) return 0; // Division by zero check
                double l_ji = (a_ji - sum) / U->elem[i - 1].data;
                if (l_ji != 0) {
                    if (L->tu >= MaxNum) return 0; // No space to add new element
                    L->elem[L->tu].i = j;
                    L->elem[L->tu].j = i;
                    L->elem[L->tu].data = l_ji;
                    L->tu++;
                }
            }
        }
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
};
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