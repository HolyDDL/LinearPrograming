#include<stdio.h>
#include<stdlib.h>

int main(){
    int m[5] = {0,1,2,3,4};
    int* n = (int*)calloc(5, sizeof(int));
    n = m;
    for(int i = 0; i<5; i++){
        printf("%d ", n[i]);
    }
}