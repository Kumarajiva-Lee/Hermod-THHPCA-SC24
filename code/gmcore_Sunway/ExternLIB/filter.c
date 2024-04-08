#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>
#include<crts.h>

#include"../LIB/Process.h"
#include"../LIB/Memory.h"
#include"../src/physical_variable.h"
#include"filter.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

typedef struct{
    int kds;
    int kde;
    int jds;
    int jde;
    int ids;
    int ide;
    int lon;
    int lat;
    int * ngrid;
    double ** wgt;
    double * input_x;
    double * x;
    double * x_bak;
    int pid;
}s_filter;

extern SLAVE_FUN(filter_kernel)(void*);
extern SLAVE_FUN(filter_kernel_1)(void*);
extern SLAVE_FUN(filter_kernel_ex)(void*);
extern SLAVE_FUN(filter_kernel_ex_1)(void*);

double max_wave_speed       = 300;
double max_cfl              = 0.5;
double filter_coef_a        = 1.0;
double filter_coef_b        = 0.4;
double filter_coef_c        = 0.2;
  
double filter_coef1_phs     = 0.375;
double filter_coef2_phs     = 0;
double filter_lat0_phs      = 60.0;
double filter_coef1_pt      = 0.375;
double filter_coef2_pt      = 0;
double filter_lat0_pt       = 60.0;
double filter_coef1_uv      = 0.375;
double filter_coef2_uv      = 0;
double filter_lat0_uv       = 60.0;
double pi2                  = 3.1415926535897 * 2;

double exp_two_values(double val1, double val2, double x0, double x1, double x){
    double w,res;

    w = exp(pow((x - x0),2) * log(1e-3) / pow((x1-x0),2));
    res = w * val1 + (1 - w) * val2;
    return res;
}

void gaussian_weight(double width, int ngrid, double* w){
    double s,x,sum;
    double vexp;
    int i,p;

    s = filter_coef_b * width / 2.0;
    p = 0;
    sum = 0;

    for (i = 0; i < ngrid; i++){
        x = i + 1 - (ngrid + 1) / 2;
        vexp = -x*x / (2 * s * s);
        if (vexp < -500) *(w + i) = 0;
        else *(w + i) = exp(vexp) / (s * sqrt(pi2));
        // *(w + i) = exp(-x*x / (2 * s * s)) / (s * sqrt(pi2));
        sum += *(w+i);
    }

    for (i = 0; i < ngrid; i++){
        *(w + i) = *(w + i) / sum;
    }
    
}

void Filter_Init(struct HybridMeshField* mesh, double dt){
    int j,p;
    int full_jds_no_pole,full_jde_no_pole,half_jds,half_jde;
    double dx,dy,cfl,w;
    int n;
    int maxlon,maxlat;
    double lat0;

    full_jds_no_pole = Has_South_Pole()?1:0;
    full_jde_no_pole = Has_North_Pole()?Proc.full_nlat-1:Proc.full_nlat;

    half_jds = 0;
    half_jde = Proc.half_nlat;

    for (p = 0 ; p < 4 ; p++){
        Filter[p].width_lon = allocate_1d_array_D(Proc.full_nlat + 2*Proc.lat_hw);
        Filter[p].ngrid_lon = allocate_1d_array_I(Proc.full_nlat + 2*Proc.lat_hw);
        Filter[p].width_lat = allocate_1d_array_D(Proc.half_nlat + 2*Proc.lat_hw);
        Filter[p].ngrid_lat = allocate_1d_array_I(Proc.half_nlat + 2*Proc.lat_hw);
    }

    maxlon = 0;
    maxlat = 0;

    if (max_wave_speed > 0){
        for (p = 0 ;  p < 4 ; p++){
            for (j = full_jds_no_pole + Proc.lat_hw; j < full_jde_no_pole + Proc.lat_hw; j++){
                dx = mesh->de_lon[0][j][0];
                dy = mesh->le_lon[0][j][0];
                if (dx > 0){
                    cfl = max_wave_speed * dt / dx;
                    w = filter_coef_a * cfl / max_cfl * (filter_coef_c * (tanh(90-fabs(mesh->full_lat_deg[0][j][0])) - 1) + 1);
                    n = ceil(w) + 2;
                    if (n % 2 == 0){
                        n += 1;
                    }
                }
                Filter[p].width_lon[j] = w;
                Filter[p].ngrid_lon[j] = n;
                maxlon = MAX(maxlon,n);
            }

            for (j = half_jds+Proc.lat_hw; j < half_jde+Proc.lat_hw; j++){
                dx = mesh->le_lat[0][j][0];
                dy = mesh->de_lat[0][j][0];
                if (dx > 0){
                    cfl = max_wave_speed * dt / dx;
                    w = filter_coef_a * cfl / max_cfl * (filter_coef_c * (tanh(90 - fabs(mesh->half_lat_deg[0][j][0])) - 1) + 1);
                    n = ceil(w) + 2;
                    if (n % 2 == 0){
                        n += 1;
                    }
                }
                Filter[p].width_lat[j] = w;
                Filter[p].ngrid_lat[j] = n;
                maxlat = MAX(maxlat,n);
            }
        }
    }

    for (p = 0 ; p < 4 ; p++){
        Filter[p].wgt_lon = allocate_2d_array_D(Proc.full_nlat + 2*Proc.lat_hw,maxlon);
        Filter[p].wgt_lat = allocate_2d_array_D(Proc.half_nlat + 2*Proc.lat_hw,maxlat);
    }  

    //Big
    for (j = full_jds_no_pole + Proc.lat_hw; j < full_jde_no_pole + Proc.lat_hw; j++){
        if (Filter[0].ngrid_lon[j] > 1){
            gaussian_weight(Filter[0].width_lon[j], Filter[0].ngrid_lon[j], Filter[0].wgt_lon[j]);
        }
    }
    for (j = half_jds + Proc.lat_hw; j <= half_jde + Proc.lat_hw; j++){
        if (Filter[0].ngrid_lat[j] > 1){
            gaussian_weight(Filter[0].width_lat[j], Filter[0].ngrid_lat[j], Filter[0].wgt_lat[j]);
        }
    }
    //phs
    lat0 = - global_mesh[0].full_lat_deg[0][Proc.lat_hw+1][0];
    for (j = full_jds_no_pole + Proc.lat_hw; j < full_jde_no_pole + Proc.lat_hw; j++){
        if (Filter[1].ngrid_lon[j] > 1){
            w = exp_two_values(filter_coef1_phs, filter_coef2_phs, lat0, filter_lat0_phs, fabs(mesh->full_lat_deg[0][j][0]));
            n = ceil(w * Filter[1].ngrid_lon[j]) + 1;
            if (n % 2 == 0) n++;
            Filter[1].ngrid_lon[j] = n;
            Filter[1].width_lon[j] = w * Filter[1].width_lon[j];
            gaussian_weight(Filter[1].width_lon[j], Filter[1].ngrid_lon[j], Filter[1].wgt_lon[j]);
        }
    }
    //pt
    lat0 = - global_mesh[0].full_lat_deg[0][Proc.lat_hw+1][0];
    for (j = full_jds_no_pole + Proc.lat_hw; j < full_jde_no_pole + Proc.lat_hw; j++){
        if (Filter[2].ngrid_lon[j] > 1){
            w = exp_two_values(filter_coef1_pt, filter_coef2_pt, lat0, filter_lat0_pt, fabs(mesh->full_lat_deg[0][j][0]));
            n = ceil(w * Filter[2].ngrid_lon[j]) + 1;
            if (n % 2 == 0) n++;
            Filter[2].ngrid_lon[j] = n;
            Filter[2].width_lon[j] = w * Filter[2].width_lon[j];
            gaussian_weight(Filter[2].width_lon[j], Filter[2].ngrid_lon[j], Filter[2].wgt_lon[j]);
        }
    }
    //uv
    lat0 = - global_mesh[0].full_lat_deg[0][Proc.lat_hw+1][0];
    for (j = full_jds_no_pole + Proc.lat_hw; j < full_jde_no_pole + Proc.lat_hw; j++){
        if (Filter[3].ngrid_lon[j] > 1){
            w = exp_two_values(filter_coef1_uv, filter_coef2_uv, lat0, filter_lat0_uv, fabs(mesh->full_lat_deg[0][j][0]));
            n = ceil(w * Filter[3].ngrid_lon[j]) + 1;
            if (n % 2 == 0) n++;
            Filter[3].ngrid_lon[j] = n;
            Filter[3].width_lon[j] = w * Filter[3].width_lon[j];
            gaussian_weight(Filter[3].width_lon[j], Filter[3].ngrid_lon[j], Filter[3].wgt_lon[j]);
        }
    }
    lat0 = - global_mesh[0].half_lat_deg[0][Proc.lat_hw][0];
    for (j = half_jds + Proc.lat_hw; j < half_jde + Proc.lat_hw; j++){
        if (Filter[3].ngrid_lat[j] > 1){
            w = exp_two_values(filter_coef1_uv, filter_coef2_uv, lat0, filter_lat0_uv, fabs(mesh->half_lat_deg[0][j][0]));
            n = ceil(w * Filter[3].ngrid_lat[j]) + 1;
            if (n % 2 == 0) n++;
            Filter[3].ngrid_lat[j] = n;
            Filter[3].width_lat[j] = w * Filter[3].width_lat[j];
            gaussian_weight(Filter[3].width_lat[j], Filter[3].ngrid_lat[j], Filter[3].wgt_lat[j]);
        }
    }
    if (Proc.id == 0){
        printf("%d %d %d\n",full_jde_no_pole,half_jde,Proc.full_nlat);
    }
    for (p = 0 ; p < 4 ; p++)
        if (Proc.id == 0){
            printf("%d %d\n",Filter[p].ngrid_lon[full_jds_no_pole + Proc.lat_hw], Filter[p].ngrid_lat[half_jds + Proc.lat_hw]);
            printf("%.8f %.8f\n",Filter[p].wgt_lon[full_jds_no_pole + Proc.lat_hw][0],Filter[p].wgt_lon[full_jds_no_pole + Proc.lat_hw][Filter[p].ngrid_lon[full_jds_no_pole + Proc.lat_hw ]-1]);
            printf("%.8f %.8f\n",Filter[p].wgt_lat[half_jds + Proc.lat_hw][0],Filter[p].wgt_lat[half_jds + Proc.lat_hw][Filter[p].ngrid_lat[half_jds + Proc.lat_hw]-1]);
        }
}

void Filter_On_Cell(int id, bool is3d, double *x){
    int lon,lat,lev;
    int full_ids, full_ide, full_jds, full_jde, full_kds, full_kde;
    int i,j,k,l;
    int n,hn,p;
    

    lon = Proc.full_nlon + 2 * Proc.lon_hw;
    full_ids = 0 + Proc.lon_hw;
    full_ide = Proc.full_nlon + Proc.lon_hw;

    lat = Proc.full_nlat + 2 * Proc.lat_hw;
    full_jds = 0 + Proc.lat_hw;
    full_jde = Proc.full_nlat + Proc.lat_hw;

    double tmp[lon];


    if (is3d){
        full_kds = 0 + Proc.lev_hw;
        full_kde = Proc.full_nlev + Proc.lev_hw;
        s_filter f;
        f.kds = full_kds;
        f.kde = full_kde;
        f.jds = full_jds;
        f.jde = full_jde;
        f.ids = full_ids;
        f.ide = full_ide;
        f.lon = lon;
        f.lat = lat;
        f.ngrid = Filter[id].ngrid_lon;
        f.wgt = Filter[id].wgt_lon;
        f.x = x;
        f.pid = Proc.id;
#ifdef slave_acc 
        if(Proc.full_nlev==32)
            athread_spawn(filter_kernel_ex_1, &f);
        else
			athread_spawn(filter_kernel_ex, &f);
        athread_join();
#else
        for (k = full_kds; k < full_kde; k++)
          for (j = full_jds; j< full_jde; j++)
            if (Filter[id].ngrid_lon[j] > 1){
              memset(tmp,0,sizeof(tmp));
              n = Filter[id].ngrid_lon[j];
              hn = (n - 1) / 2;
              for (i = full_ids; i < full_ide; i++){
                p = 0;
                for (l = i - hn; l <= i + hn; l++){
                  tmp[i] += Filter[id].wgt_lon[j][p++] * (*(x + k*lon*lat + j*lon + l));
                  //if (Proc.id == 19 && is3d == 1 && i == full_ids && j == full_jde - 1 && k == full_kds) printf("%.8f %.8f\n",Filter[id].wgt_lon[j][p-1],*(x + k*lon*lat + j*lon + l));
                }
              }
              
              for (i = full_ids; i < full_ide; i++)
                *(x + k*lon*lat + j*lon + i) = tmp[i];
            }
#endif
    }
    else{
        full_kds = 0;
        full_kde = 1;

        for (k = full_kds; k < full_kde; k++)
          for (j = full_jds; j< full_jde; j++)
            if (Filter[id].ngrid_lon[j] > 1){
              memset(tmp,0,sizeof(tmp));
              n = Filter[id].ngrid_lon[j];
              hn = (n - 1) / 2;
              for (i = full_ids; i < full_ide; i++){
                p = 0;
                for (l = i - hn; l <= i + hn; l++){
                  tmp[i] += Filter[id].wgt_lon[j][p++] * (*(x + k*lon*lat + j*lon + l));
                  //if (Proc.id == 19 && is3d == 1 && i == full_ids && j == full_jde - 1 && k == full_kds) printf("%.8f %.8f\n",Filter[id].wgt_lon[j][p-1],*(x + k*lon*lat + j*lon + l));
                }
              }
              
              for (i = full_ids; i < full_ide; i++)
                *(x + k*lon*lat + j*lon + i) = tmp[i];
            }
    }

}

void Filter_on_lon_edge(int id, double *x){
    int lon,lat,lev;
    int half_ids, half_ide, full_jds, full_jde, full_kds, full_kde;
    int i,j,k,l;
    int n,hn,p;
    

    lon = Proc.half_nlon + 2 * Proc.lon_hw;
    half_ids = 0 + Proc.lon_hw;
    half_ide = Proc.half_nlon + Proc.lon_hw;

    lat = Proc.full_nlat + 2 * Proc.lat_hw;
    full_jds = 0 + Proc.lat_hw;
    full_jde = Proc.full_nlat + Proc.lat_hw;

    full_kds = 0 + Proc.lev_hw;
    full_kde = Proc.full_nlev + Proc.lev_hw;


#ifdef slave_acc 
    s_filter f;
    f.kds = full_kds;
    f.kde = full_kde;
    f.jds = full_jds;
    f.jde = full_jde;
    f.ids = half_ids;
    f.ide = half_ide;
    f.lon = lon;
    f.lat = lat;
    f.ngrid = Filter[id].ngrid_lon;
    f.wgt = Filter[id].wgt_lon;
    f.x = x;
    f.pid = Proc.id;
    if(Proc.full_nlev==32)
        athread_spawn(filter_kernel_ex_1, &f);
    else
		athread_spawn(filter_kernel_ex, &f);
    athread_join();

#else
    double tmp[lon];
    for (k = full_kds; k < full_kde; k++)
      for (j = full_jds; j< full_jde; j++)
        if (Filter[id].ngrid_lon[j] > 1){
          memset(tmp,0,sizeof(tmp));
          n = Filter[id].ngrid_lon[j];
          hn = (n - 1) / 2;
          for (i = half_ids; i < half_ide; i++){
            p = 0;
            for (l = i - hn; l <= i + hn; l++)
              tmp[i] += Filter[id].wgt_lon[j][p++] * (*(x + k*lon*lat + j*lon + l));
          }
          for (i = half_ids; i < half_ide; i++)
            *(x + k*lon*lat + j*lon + i) = tmp[i];
        }
#endif
}

void Filter_on_lat_edge(int id, double *x){
    int lon,lat,lev;
    int full_ids, full_ide, half_jds, half_jde, full_kds, full_kde;
    int i,j,k,l;
    int n,hn,p;
    

    lon = Proc.full_nlon + 2 * Proc.lon_hw;
    full_ids = 0 + Proc.lon_hw;
    full_ide = Proc.full_nlon + Proc.lon_hw;

    lat = Proc.half_nlat + 2 * Proc.lat_hw;
    half_jds = 0 + Proc.lat_hw;
    half_jde = Proc.half_nlat + Proc.lat_hw;

    full_kds = 0 + Proc.lev_hw;
    full_kde = Proc.full_nlev + Proc.lev_hw;


#ifdef slave_acc 
    s_filter f;
    f.kds = full_kds;
    f.kde = full_kde;
    f.jds = half_jds;
    f.jde = half_jde;
    f.ids = full_ids;
    f.ide = full_ide;
    f.lon = lon;
    f.lat = lat;
    f.ngrid = Filter[id].ngrid_lat;
    f.wgt = Filter[id].wgt_lat;
    f.x = x;
    f.pid = Proc.id;
    if(Proc.full_nlev==32)
		athread_spawn(filter_kernel_ex_1, &f);
    else
		athread_spawn(filter_kernel_ex, &f);
    athread_join();

#else

    //memset(tmp,0,sizeof(tmp));
    double tmp[lon];
    memset(tmp,0,sizeof(tmp));
    for (k = full_kds; k < full_kde; k++)
      for (j = half_jds; j< half_jde; j++)//lat
        if (Filter[id].ngrid_lat[j] > 1){
          memset(tmp,0,sizeof(tmp));
          n = Filter[id].ngrid_lat[j];
          hn = (n - 1) / 2;
          for (i = full_ids; i < full_ide; i++){//lon
            p = 0;
            for (l = i - hn; l <= i + hn; l++){
              tmp[i] += Filter[id].wgt_lat[j][p++] * (*(x + k*lon*lat + j*lon + l));
            }
          }
          for (i = full_ids; i < full_ide; i++)
            *(x + k*lon*lat + j*lon + i) = tmp[i];
        }
#endif


}

