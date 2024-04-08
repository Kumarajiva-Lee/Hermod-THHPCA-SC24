#ifndef NC_H_INCLUDED
#define NC_H_INCLUDED 1

#include<stdio.h>
#include<stdbool.h>
#include<string.h>
#include<netcdf.h>

int ncOpenFile(char *file_name, int is_write);
void ncDefDim(int nid, int num_lev, int num_lat, int num_lon);
void ncDefVar(int nid, char *var_name, int has_lev, int full_lat, int full_lon);
void ncGetVar(int nid, char *var_name, int lev_beg, int lev_cnt, int lev_hw, int lat_beg, int lat_cnt, int lat_hw, int lon_beg, int lon_cnt, int lon_hw, double *data);
void ncPutVar(int nid, char *var_name, int lev_beg, int lev_cnt, int lev_hw, int lat_beg, int lat_cnt, int lat_hw, int lon_beg, int lon_cnt, int lon_hw, double *data);
void ncCloseFile(int ncid);

#endif