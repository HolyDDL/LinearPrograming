#include<stdio.h>

#include"utils.h"

int main(){
    double** ori_matrix;
    printf("input rows number(m): ");
    int m;
    scanf("%d", &m);
    printf("input cols number(n): ");
    int n;
    scanf("%d", &n);
    ori_matrix = Get_Matrix(m, n);
    double* obj;
    obj = Get_objfunc(n);
    double** UnitMatrix;
    UnitMatrix = Judge_Matrix(ori_matrix, m, n);
    double z0 = 0;
    SimplexAlgorithm(ori_matrix, obj, m, n, &z0, obj);
    printf("min z = %lf\n", z0);
}