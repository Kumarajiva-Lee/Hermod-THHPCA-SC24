#include<stdlib.h>
#include<string.h>

#include "Memory.h"

//按C的连续维度排列，x为最外层，z为最内层
int ***allocate_3d_array_I(int x, int y, int z) {
    int ***arr = (int ***) malloc(sizeof(int **) * x);
    int *data = (int *) malloc(sizeof(int) * x * y * z);
    memset(data,0,sizeof(int) * x * y * z);
    for (int i = 0; i < x; i++) {
        arr[i] = (int **) malloc(sizeof(int *) * y);
        for (int j = 0; j < y; j++) {
            arr[i][j] = &data[(i * y * z) + (j * z)];
        }
    }
    return arr;
}

double ***allocate_3d_array_D(int x, int y, int z) {
    double ***arr = (double ***) malloc(sizeof(double **) * x);
    double *data = (double *) malloc(sizeof(double) * x * y * z);
    memset(data,0,sizeof(double) * x * y * z);
    for (int i = 0; i < x*y*z; i++) {
        data[i] = 0;
    }
    for (int i = 0; i < x; i++) {
        arr[i] = (double **) malloc(sizeof(double *) * y);
        for (int j = 0; j < y; j++) {
            arr[i][j] = &data[(i * y * z) + (j * z)];
        }
    }
    return arr;
}

void free_3d_array_I(int ***arr, int x, int y, int z) {
    free(arr[0][0]);
    for (int i = 0; i < x; i++) {
        free(arr[i]);
    }
    free(arr);
}

void free_3d_array_D(double ***arr, int x, int y, int z) {
    free(arr[0][0]);
    for (int i = 0; i < x; i++) {
        free(arr[i]);
    }
    free(arr);
}

double **allocate_2d_array_D(int x, int y) {
    double **arr = (double **) malloc(sizeof(double *) * x);
    double *data = (double *) malloc(sizeof(double) * x * y);
    memset(data,0,sizeof(double) * x * y);
    for (int i = 0; i < x; i++) {
        arr[i] = &data[i * y];
    }
    return arr;
}

void free_2d_array_D(double **arr, int x, int y) {
    free(arr[0]);
    free(arr);
}

double *allocate_1d_array_D(int x) {
    double *arr = (double*) malloc(sizeof(double) * x);
    memset(arr,0,sizeof(double) * x);
    return arr;
}

int *allocate_1d_array_I(int x) {
    int *arr = (int*) malloc(sizeof(int) * x);
    memset(arr,0,sizeof(int) * x);
    return arr;
}

void free_1d_array_D(double *arr, int x) {
    free(arr);
}

void free_1d_array_I(int *arr, int x) {
    free(arr);
}