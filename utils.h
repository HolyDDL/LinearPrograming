#include<stdio.h>
#include<stdlib.h>

double* Get_objfunc(int n){
    printf("input object func: \n");
    double* obj = (double*)malloc(n*sizeof(double));
    for(int i =0; i<n ;i++){
        double temp;
        scanf("%lf", &temp);
        obj[i] = temp;
    }
    printf("end input of obj func. \n");
    return obj;
}

double** Get_Matrix(int m, int n){
    printf("------INPUT MATRIX-----\n");
    
    double** rowhead = (double**)malloc(m*sizeof(double*));
    for(int i = 0; i<m; i++){
        printf("input the %d row of matrix(include b vector): \n", i);
        rowhead[i] = (double*)malloc(n*sizeof(double));
        for(int j = 0; j<n; j++){
            double temp;
            scanf("%lf", &temp);
            rowhead[i][j] = temp;
        }
        printf("finish the input of the %d row.\n", i);
    }
    printf("------END INPUT------\n");
    return rowhead;
}

// use it to judge whether if the matrix need support var , return UnitVarMatrix
double** Judge_Matrix(double** matrix, int m, int n){
    double** UnitMatrix = (double**)malloc(n*sizeof(double*));
    for(int b = 0; b<n ;b++)
        UnitMatrix[b] = NULL;
    int unitnum = 0;
    for(int j = 0; j<n; j++){
        //every cols
        int isUnitVector = 1;
        double* col = (double*)malloc(m*sizeof(double));
        for(int i = 0; i<m; i++){
            col[i] = matrix[i][j];
        }
        int OneNum = 0;
        for(int k = 0; k<m; k++){
            //judge col
            if(!(col[k]>1.0||col[k]<1.0))
                OneNum++;
            if(((col[k]>0.0 || col[k]<0.0) && (col[k]>1.0 || col[k]<1.0) )|| OneNum > 1){
                isUnitVector = 0;
                break;
            }
        }
        if(isUnitVector){
            UnitMatrix[j] = col;
            unitnum++;
        }
    }
    printf("this Matrix has %d unit vectors. \n", unitnum);
    return UnitMatrix;
}

int SimplexAlgorithm(double** Matrix, double* obj, int m, int n, double* minz, double* vector){
    typedef struct {
        int isUnit;
        int isZzero; // judge whether if z is zero
        int unitROW;
    }Unit;
    double** LPMatrix = (double**)malloc((m+1)*sizeof(double*));
    double* tempobj = (double*)malloc(n*sizeof(double));
    tempobj = obj;
    LPMatrix[0] = tempobj;
    for(int i = 1; i<m+1; i++){
        double* temp = (double*)malloc(n*sizeof(double));
        LPMatrix[i] = temp;
    }
    Unit* Unitcols = (Unit*)calloc((n-1), sizeof(Unit));
    for(int col = 0; col<n-1; col++){
        // this loop is used to judge Unit Vectors
        int isUnit = 1;
        int unitROW = 0;
        for(int row = 1; row<m+1; row++){
            int OneNum = 0;
            if(!(LPMatrix[row][col]>1.0||LPMatrix[row][col]<1.0)){
                unitROW = row;
                OneNum++;
            }
            if(((LPMatrix[row][col]>0.0 || LPMatrix[row][col]<0.0) && (LPMatrix[row][col]>1.0 || LPMatrix[row][col]<1.0) )|| OneNum > 1){
                isUnit = 0;
                unitROW = 0;
                break;
            }
        }
        if(isUnit){
            Unitcols[col].isUnit = 1;
            Unitcols[col].unitROW = unitROW;
            if(!(LPMatrix[0][col] > 0.0 || LPMatrix[0][col] < 0.0))
                Unitcols[col].isZzero = 1;
        }
    }//end unit loop
    int* haveonenum = (int*)calloc(m, sizeof(int));
    for(int col = 0; col<n-1; col++){
        //delete extra cols
        if(Unitcols[col].isUnit == 1){
            haveonenum[Unitcols[col].unitROW]++;
            if(haveonenum[Unitcols[col].unitROW] > 1){
                Unitcols[col].isUnit = 0;
            }
        }
    }// end delete extra cols loop
    for(int col = 0; col<n-1; col++){
        //Init Simplex Table
        if(Unitcols[col].isUnit == 1 && Unitcols[col].isZzero != 1){
            double* temprow = (double*)malloc(n*sizeof(double));
            temprow = LPMatrix[Unitcols[col].unitROW];
            for(int i=0; i<n; i++){
                temprow[i] *= -LPMatrix[0][col];
            }
            for(int i=0; i<n; i++){
                LPMatrix[0][i] += temprow[i];
            }
        }
    } // end init
    int finish = 1;
    for(int times = 0; times<n*m; times++){
        //Simplex
        int neigative_z = 1;
        int max_z = LPMatrix[0][0];
        int max_z_col = 0;
        for(int col = 0; col<n-1; col++){
            if(LPMatrix[0][col] > 0.0)
                neigative_z = 0;
            if(LPMatrix[0][col] > max_z){
                max_z = LPMatrix[0][col];
                max_z_col = col;
            }
        }
        if(neigative_z)
            break;
        int all_neigaative_a = 1;
        int min_b_a = -1;
        int min_row = 1;
        for(int row = 1; row<m+1; row++){
            if(LPMatrix[row][max_z_col] > 0.0)
                all_neigaative_a = 0;
            if((LPMatrix[row][n-1] / LPMatrix[row][max_z_col] < min_b_a) || min_b_a == -1){
                min_b_a = LPMatrix[row][n-1] / LPMatrix[row][max_z_col];
                min_row = row;
            }
        }
        if(all_neigaative_a){
            printf("NO ANSWER!\n");
            return 0;
        }
        for(int col = 0; col<n; col++){
            LPMatrix[min_row][col] /= LPMatrix[min_row][max_z_col];
        }
        for(int row = 0; row<m+1; row++){
            if(LPMatrix[row][max_z_col] > 0.0){
                double* temprow = (double*)malloc(n*sizeof(double));
                for(int i = 0; i<n; i++){
                    temprow[i] = -LPMatrix[row][max_z_col]*LPMatrix[min_row][i];
                }
                for(int col = 0; col<n; col++){
                    LPMatrix[row][col] += temprow[col];
                }
            }
        }
        if(times == n*m-1)
            finish = 0;
    }
    if(!finish)
        printf("######## NOT FINISH! ########\n");
    *minz = LPMatrix[0][n-1];
}