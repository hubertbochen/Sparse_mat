#include<stdio.h>
#define NUM 10
//三元组
typedef struct {
    int i;
    int j;
    int data;
}triple;
//三元组顺序表
typedef struct {
    triple data[NUM];
    int mu;
    int nu;
    int tu;
}TSMatrix;
//稀疏矩阵的转置
void transposeMatrix(TSMatrix M, TSMatrix* T) {
    //1.稀疏矩阵的行数和列数互换
    (*T).mu = M.nu;
    (*T).nu = M.mu;
    (*T).tu = M.tu;
    
    if ((*T).tu) {
        int col, p;
        int q = 0;
        //2.遍历原表中的各个三元组
        for (col = 1; col <= M.nu; col++) {
            //重复遍历原表 M.m 次，将所有三元组都存放到新表中
            for (p = 0; p < M.tu; p++) {
                //3.每次遍历，找到 j 列值最小的三元组，将它的行、列互换后存储到新表中
                if (M.data[p].j == col) {
                    (*T).data[q].i = M.data[p].j;
                    (*T).data[q].j = M.data[p].i;
                    (*T).data[q].data = M.data[p].data;
                    //为存放下一个三元组做准备
                    q++;
                }
            }
        }
    }
}

int main() {
    int i, k;
    TSMatrix M, T;
    M.mu = 3;
    M.nu = 2;
    M.tu = 4;
    
    M.data[0].i = 1;
    M.data[0].j = 2;
    M.data[0].data = 1;
    
    M.data[1].i = 2;
    M.data[1].j = 2;
    M.data[1].data = 3;
    
    M.data[2].i = 3;
    M.data[2].j = 1;
    M.data[2].data = 6;
    
    M.data[3].i = 3;
    M.data[3].j = 2;
    M.data[3].data = 5;
    
    for (k = 0; k < NUM; k++) {
        T.data[k].i = 0;
        T.data[k].j = 0;
        T.data[k].data = 0;
    }
    transposeMatrix(M, &T);
    for (i = 0; i < T.tu; i++) {
        printf("(%d,%d,%d)\n", T.data[i].i, T.data[i].j, T.data[i].data);
    }
    return 0;
}