#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED 1

int ***allocate_3d_array_I(int x, int y, int z);
double ***allocate_3d_array_D(int x, int y, int z);
double **allocate_2d_array_D(int x, int y);
double *allocate_1d_array_D(int x);
int *allocate_1d_array_I(int x);

void free_3d_array_I(int ***arr, int x, int y, int z);
void free_3d_array_D(double ***arr, int x, int y, int z);
void free_2d_array_D(double **arr, int x, int y);
void free_1d_array_D(double *arr, int x);
void free_1d_array_I(int *arr, int x);

#endif