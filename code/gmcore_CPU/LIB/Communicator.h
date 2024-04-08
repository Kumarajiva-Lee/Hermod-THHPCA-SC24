#ifndef COMMUNICATORS_H_INCLUDED
#define COMMUNICATORS_H_INCLUDED

#include<stdio.h>
#include<stdbool.h>
#include<string.h>

#include"mpi.h"
#include"ParaParam.h"
#include"Process.h"
#include"Async.h"


void UpdateHalo_2d_I(ProcType Proc, int *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo);
void UpdateHalo_3d_I(ProcType Proc, int *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo);
void UpdateHalo_2d_S(ProcType Proc, float *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo);
void UpdateHalo_3d_S(ProcType Proc, float *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo);
void UpdateHalo_2d_D(ProcType Proc, double *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo, bool small);
void UpdateHalo_3d_D(ProcType Proc, double *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo, bool small);

void UpdateHaloCS_3d_D(ProcType Proc, double *array, double *arrayX, double *arrayY, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo);

void HaloWait(ProcType Proc, AsyncType *async);

//void Zonal_Sum_1d_D(ProcType Proc, double* sum, double *res,int size);
void Zonal_Sum_1d_D(ProcType Proc, double** sum, int size);

#endif