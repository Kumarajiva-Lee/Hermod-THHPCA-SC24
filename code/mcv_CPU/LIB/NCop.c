#include<stdlib.h>
#include<string.h>
#include<netcdf.h>

#include "NCop.h"


int ncOpenFile(char *file_name, int is_write)
{
    int ncid, res;
    if (is_write == 0) {
        res = nc_open(file_name, NC_NOWRITE || NC_MPIIO, &ncid);
    } else {
        res = nc_create(file_name, NC_NETCDF4 || NC_MPIIO, &ncid);
        nc_enddef(ncid);
    }
    return ncid;
}

void ncGetVar(int nid, char *var_name, int lev_beg, int lev_cnt, int lev_hw, int lat_beg, int lat_cnt, int lat_hw, int lon_beg, int lon_cnt, int lon_hw, double *data)
{
    int varid, res;
    long strd[4];
    unsigned long start[4], cnt[4];
    int i,j,k,tot,tid;
    tot = 0;
    strd[0] = 1;
    strd[1] = 1;
    strd[2] = 1;
    strd[3] = 1;
    start[0] = 0;
    start[1] = lev_beg;
    start[2] = lat_beg;
    start[3] = lon_beg;
    cnt[0] = 1;
    cnt[1] = lev_cnt;
    cnt[2] = lat_cnt;
    cnt[3] = lon_cnt;
    res = nc_inq_varid(nid, var_name, &varid);
    double *tmp;
    tmp = (double*)malloc(sizeof(double)*cnt[1]*cnt[2]*cnt[3]);
    long nstrd[3];
    unsigned long nstart[3], ncnt[3];
    if (cnt[1] == 1) {
        nstrd[0] = 1;
        nstrd[1] = 1;
        nstrd[2] = 1;
        ncnt[0] = 1;
        ncnt[1] = lat_cnt;
        ncnt[2] = lon_cnt;
        nstart[0] = 0;
        nstart[1] = lat_beg;
        nstart[2] = lon_beg;
        res = nc_get_vars(nid, varid, nstart, ncnt, nstrd, tmp);
    }
    else {
        res = nc_get_vars(nid, varid, start, cnt, strd, tmp);
    }
    
    for (i = 0; i < lev_cnt; i++)
        for (j = 0; j < lat_cnt; j++)
            for (k = 0; k < lon_cnt; k++) {
                tid = (i+lev_hw)*(lat_cnt+2*lat_hw)*(lon_cnt+2*lon_hw) + (j+lat_hw)*(lon_cnt+2*lon_hw) + k+lon_hw;
                data[tid] = tmp[tot];
                tot++;
            }
    free(tmp);
}

void ncDefDim(int nid, int num_lev, int num_lat, int num_lon)
{
    int temp_id;
    nc_redef(nid);
    nc_def_dim(nid, "time", NC_UNLIMITED, &temp_id);
    nc_def_dim(nid, "lev", num_lev, &temp_id);
    nc_def_dim(nid, "lon", num_lon, &temp_id);
    nc_def_dim(nid, "lat", num_lat, &temp_id);
    nc_def_dim(nid, "ilon", num_lon, &temp_id);
    nc_def_dim(nid, "ilat", num_lat - 1, &temp_id);
    nc_enddef(nid);
}

void ncDefVar(int nid, char *var_name, int has_lev, int full_lat, int full_lon)
{
    nc_redef(nid);
    int dim_id, var_id;
    if (has_lev == 1) {
        int temp_ids[4];
        nc_inq_dimid(nid, "time", &dim_id);
        temp_ids[0] = dim_id;
        nc_inq_dimid(nid, "lev", &dim_id);
        temp_ids[1] = dim_id;
        if (full_lat == 1) {
            nc_inq_dimid(nid, "lat", &dim_id);
        } else {
            nc_inq_dimid(nid, "ilat", &dim_id);
        }
        temp_ids[2] = dim_id;
        if (full_lon == 1) {
            nc_inq_dimid(nid, "lon", &dim_id);
        } else {
            nc_inq_dimid(nid, "ilon", &dim_id);
        }
        temp_ids[3] = dim_id;
        nc_def_var(nid, var_name, NC_DOUBLE, 4, temp_ids, &var_id);
    }
    else {
        int temp_ids[3];
        nc_inq_dimid(nid, "time", &dim_id);
        temp_ids[0] = dim_id;
        if (full_lat == 1) {
            nc_inq_dimid(nid, "lat", &dim_id);
        } else {
            nc_inq_dimid(nid, "ilat", &dim_id);
        }
        temp_ids[1] = dim_id;
        if (full_lon == 1) {
            nc_inq_dimid(nid, "lon", &dim_id);
        } else {
            nc_inq_dimid(nid, "ilon", &dim_id);
        }
        temp_ids[2] = dim_id;
        nc_def_var(nid, var_name, NC_DOUBLE, 3, temp_ids, &var_id);
    }
    nc_enddef(nid);
}

void ncPutVar(int nid, char *var_name, int lev_beg, int lev_cnt, int lev_hw, int lat_beg, int lat_cnt, int lat_hw, int lon_beg, int lon_cnt, int lon_hw, double *data)
{
    int varid, res;
    long strd[4];
    unsigned long start[4], cnt[4];
    int i,j,k,tot,tid;
    tot = 0;
    strd[0] = 1;
    strd[1] = 1;
    strd[2] = 1;
    strd[3] = 1;
    start[0] = 0;
    start[1] = lev_beg;
    start[2] = lat_beg;
    start[3] = lon_beg;
    cnt[0] = 1;
    cnt[1] = lev_cnt;
    cnt[2] = lat_cnt;
    cnt[3] = lon_cnt;
    res = nc_inq_varid(nid, var_name, &varid);
    double *tmp;
    tmp = (double*)malloc(sizeof(double)*cnt[1]*cnt[2]*cnt[3]);
    for (i = 0; i < lev_cnt; i++)
        for (j = 0; j < lat_cnt; j++)
            for (k = 0; k < lon_cnt; k++) {
                tid = (i+lev_hw)*(lat_cnt+2*lat_hw)*(lon_cnt+2*lon_hw) + (j+lat_hw)*(lon_cnt+2*lon_hw) + k+lon_hw;
                tmp[tot] = data[tid];
                tot++;
            }
    long nstrd[3];
    unsigned long nstart[3], ncnt[3];
    if (cnt[1] == 1) {
        nstrd[0] = 1;
        nstrd[1] = 1;
        nstrd[2] = 1;
        ncnt[0] = 1;
        ncnt[1] = lat_cnt;
        ncnt[2] = lon_cnt;
        nstart[0] = 0;
        nstart[1] = lat_beg;
        nstart[2] = lon_beg;
        res = nc_put_vars(nid, varid, nstart, ncnt, nstrd, tmp);
    }
    else {
        res = nc_put_vars(nid, varid, start, cnt, strd, tmp);
        // printf("%s %d %d %d %d %d %d\n",var_name,start[0],start[1],start[2],cnt[0],cnt[1],cnt[2]);
    }
    //printf("res = %d \n",res);
    free(tmp);
}


void ncCloseFile(int ncid)
{
    nc_close(ncid);
}
