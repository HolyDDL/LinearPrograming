#include<stdio.h>
#include<stdlib.h>

typedef struct un{
    int col;
    int one_row;
    int isZzero;
    struct un* next;
}Unit;

int FindUnitVector(double** Matrix, Unit* unit, double* obj, int m, int n){
    Unit* head = unit;
    int* onerows = (int*)calloc(m, sizeof(int));
    for(int col = 0; col<n-1; col++){
        int isUnit = 1;
        int unitROW = 0;
        int OneNum = 0;
        for(int row = 0; row<m; row++){

            if(!(Matrix[row][col]>1.0||Matrix[row][col]<1.0)){
                unitROW = row;
                OneNum++;
            }
            if(((Matrix[row][col]>0.0 || Matrix[row][col]<0.0) && (Matrix[row][col]>1.0 || Matrix[row][col]<1.0) )|| OneNum > 1){
                isUnit = 0;
                unitROW = 0;
                break;
            }
        }
        if(isUnit){
            if(onerows[unitROW] == 0){
                onerows[unitROW]++;
                Unit* p = (Unit*)malloc(sizeof(Unit));
                p->col = col;
                p->one_row = unitROW;
                p->next = head->next;
                head->next = p;
                head = p;
                if((obj[col]>0.0)||(obj[col]<0.0))
                    p->isZzero = 0;
                else
                    p->isZzero = 1;
            }
        }
    }
   return 0;
}

int Simplex(double** Matrix, double* obj, int m, int n, Unit* unit){
    Unit* head = unit->next;
    while(head){
        if(head->isZzero == 0){
            double sc = obj[head->col];
            for(int i = 0; i<n; i++){
                double temprow = -sc * Matrix[head->one_row][i];
                obj[i] += temprow;
            }
        }
        head = head->next;
    }
    for(int times = 0; times < m*n+1; times++){
        int allnegative = 1;
        double max_z = 0.0;
        int max_z_col = -1;
        for(int col = 0; col<n-1; col++){
            if(obj[col] > max_z ){
                allnegative = 0;
                int negative = 1;
                for(int row=0; row<m; row++){
                    if(Matrix[row][col] > 0.0){
                        negative = 0;
                        break;
                    }
                }
                if(!negative){
                    max_z = obj[col];
                    max_z_col = col;
                }
            }
        }
        if(max_z_col == -1 && !allnegative){
            printf("NO SULOTION!\n");
            return 0;
        }
        if(allnegative){
            printf("****** SOLUTION ******\n");
            printf("The min z is %.2f.\n", obj[n-1]);
            Unit* lasthead = (Unit*)malloc(sizeof(Unit));
            lasthead->next = NULL;
            FindUnitVector(Matrix, lasthead, obj, m, n);
            Unit* p = lasthead->next;
            printf("the solution vector is: \n");
            double* vector = (double*)calloc(n, sizeof(double));
            while(p){
                vector[p->col] = Matrix[p->one_row][n-1];
                p = p->next;
            }
            for(int i = 0; i<n; i++){
                if(i == 0)
                    printf("[%.2f",vector[i]);
                else if(i == n-1)
                    printf(", %.2f]\n", vector[i]);
                else
                    printf(", %.2f", vector[i]);
            }
            return 0;
        }

        double min_ba;
        int firsttimein = 1;
        int min_row = -1;
        for(int row = 0; row<m; row++){
            if(Matrix[row][max_z_col]>0.0){
                if(firsttimein || (Matrix[row][n-1]/Matrix[row][max_z_col]<min_ba)){
                    firsttimein = 0;
                    min_ba = Matrix[row][n-1]/Matrix[row][max_z_col];
                    min_row = row;
                }
            }
        }
        double factor = Matrix[min_row][max_z_col];
        for(int i = 0; i<n; i++){
            Matrix[min_row][i] /= factor;  
        }
        double scale = obj[max_z_col];
        for(int i = 0; i<n; i++){
            double temp = -scale*Matrix[min_row][i];
            obj[i] += temp;
        }
        for(int i = 0; i<m; i++){
            if(i != min_row){
                double _scale = Matrix[i][max_z_col];
                for(int j = 0; j<n; j++){
                    double temp = -_scale * Matrix[min_row][j];
                    Matrix[i][j] += temp;
                }
            }
        }
    }
    return 0;
}

double** StageONE(double** Matrix, double* obj, int m, int n, Unit* unit, int supportnum){
    Unit* head = unit->next;
    while(head){
        if(head->isZzero == 0){
            double scale = obj[head->col];
            for(int i = 0; i<n; i++){
                double temprow = -scale * Matrix[head->one_row][i];
                obj[i] += temprow;
            }
        }
        head = head->next;
    }
    for(int times = 0; times < m*n+1; times++){
        int allnegative = 1;
        double max_z = 0.0;
        int max_z_col = -1;
        for(int col = 0; col<n-1; col++){
            if(obj[col] > max_z ){
                allnegative = 0;
                int negative = 1;
                for(int row=0; row<m; row++){
                    if(Matrix[row][col] > 0.0){
                        negative = 0;
                        break;
                    }
                }
                if(!negative){
                    max_z = obj[col];
                    max_z_col = col;
                }
            }
        }
        if(allnegative){
           if(!((obj[n-1]-0.0)<1e-7)){
                printf("NO SULOTION!\n");
                return NULL;
           }
           else{
                double** AfterMatrix = (double**)malloc(m*sizeof(double*));
                for(int i =0; i<m; i++){
                    double* line = (double*)malloc((n - supportnum)*sizeof(double));
                    for(int j = 0; j<n-supportnum-1; j++){
                        line[j] = Matrix[i][j];
                    }
                    line[n-supportnum-1] = Matrix[i][n-1];
                    AfterMatrix[i] = line;
                }
                return AfterMatrix;
           }
        }
        double min_ba;
        int firsttimein = 1;
        int min_row = -1;
        for(int row = 0; row<m; row++){
            if(Matrix[row][max_z_col]>0.0){
                if(firsttimein || (Matrix[row][n-1]/Matrix[row][max_z_col]<min_ba)){
                    firsttimein = 0;
                    min_ba = Matrix[row][n-1]/Matrix[row][max_z_col];
                    min_row = row;
                }
            }
        }
        double factor = Matrix[min_row][max_z_col];
        for(int i = 0; i<n; i++){
            Matrix[min_row][i] /= factor;  
        }
        double scale = obj[max_z_col];
        for(int i = 0; i<n; i++){
            double temp = -scale*Matrix[min_row][i];
            obj[i] += temp;
        }
        for(int i = 0; i<m; i++){
            if(i != min_row){
                double _scale = Matrix[i][max_z_col];
                for(int j = 0; j<n; j++){
                    double temp = -_scale * Matrix[min_row][j];
                    Matrix[i][j] += temp;
                }
            }
        }
    }
    return NULL;
}

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
    printf("------ INPUT MATRIX -----\n");
    
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
    printf("------ END INPUT MATRIX ------\n");
    return rowhead;
}

double TwoStageMethod(double** Matrix, int m, int n, double* obj){
    Unit* head = (Unit*)malloc(sizeof(Unit));
    head->next = NULL;
    FindUnitVector(Matrix, head, obj, m, n);
    int* UnitROW = (int*)calloc(m, sizeof(int));
    int unitrow_num = 0;
    Unit* p = head->next;
    while(p){
        UnitROW[p->one_row] = 1;
        p = p->next;
        unitrow_num++;
    }
    // int* UnitMASK = (int*)malloc(m*sizeof(int));
    // for(int i = 0; i<m; i++)
    //     UnitMASK[i] = 1 - UnitROW[i];
    int needVARSnum = m - unitrow_num;
    double** BuildMatrix = (double**)malloc(m*sizeof(double*));
    for(int i = 0; i<m; i++){
        double* line = (double*)calloc(needVARSnum, sizeof(double));
        BuildMatrix[i] = line;
    }
    for(int j = 0; j<needVARSnum; j++){
        for(int i = 0; i<m; i++){
            if(UnitROW[i] == 0){
                UnitROW[i] = 1;
                BuildMatrix[i][j] = 1.0;
                break;
            }
        }
    }
    double** LPMatrix = (double**)malloc(m*sizeof(double*));
    for(int i = 0; i<m; i++){
        double* line = (double*)malloc((n+needVARSnum)*sizeof(double));
        for(int j = 0; j<n-1; j++){
            line[j] = Matrix[i][j];
        }
        for(int j = 0; j<needVARSnum; j++){
            line[n+j-1] = BuildMatrix[i][j];
        }
        line[n+needVARSnum-1] = Matrix[i][n-1];
        LPMatrix[i] = line;
    }
    double* supportobj = (double*)calloc((n+needVARSnum), sizeof(double));
    for(int j = 0; j<needVARSnum; j++){
        supportobj[n+j-1] = -1;
    }
    int supportm = m;
    int supportn = n+needVARSnum;
    Unit* supporthead = (Unit*)malloc(sizeof(Unit));
    supporthead->next = NULL;
    FindUnitVector(LPMatrix, supporthead, supportobj, supportm, supportn);
    double** AfterMatrix = StageONE(LPMatrix, supportobj, supportm, supportn, supporthead, needVARSnum);
    if(!AfterMatrix)
        return -1;
    Unit* Afterhead = (Unit*)malloc(sizeof(Unit));
    Afterhead->next = NULL;
    FindUnitVector(AfterMatrix, Afterhead, obj, m ,n);
    Simplex(AfterMatrix, obj, m, n, Afterhead);
}

int main(){
    // int m = 2;
    // int n = 6;
    // double obj[6] = {-5, 0, -21, 0, 0, 0};
    // double line1[6] = {1, -1, 6, -1, 0, 2};
    // double line2[6] = {1, 1, 2, 0, -1, 1};
    // double** Matrix = (double**)malloc(m*sizeof(double*));
    // Matrix[0] = line1;
    // Matrix[1] = line2;
    printf("input rows number(m): ");
    int m;
    scanf("%d", &m);
    printf("input cols number(n): ");
    int n;
    scanf("%d", &n);
    double** Matrix = Get_Matrix(m, n);
    double* obj = Get_objfunc(n);
    // Unit* header = (Unit*)malloc(sizeof(Unit));
    // header->next =  NULL;
    // FindUnitVector(Matrix, header, obj, m ,n);
    // Simplex(Matrix, obj, m, n, header);
    TwoStageMethod(Matrix, m, n, obj);
}
