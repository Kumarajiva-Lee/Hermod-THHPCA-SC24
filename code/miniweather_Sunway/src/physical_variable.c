#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include"../LIB/Process.h"
#include"../LIB/Memory.h"

#include"physical_variable.h"

double TimeBeg;
double TimeEnd;
double CommTime;
double CompTime;
double MeshInitTime;
double StateInitialrhoTime;
double StateInitialTime;
double ComputeTendenciesXfrhoTime;
double ComputeTendenciesXdrhoTime;
double SetBoundaryZrhoTime;
double SetBoundaryZuTime;
double SetBoundaryZwTime;
double SetBoundaryZptTime;
double ComputeTendenciesZfrhoTime;
double ComputeTendenciesZdrhoTime;
double UpdateStaterhoTime;
double OutputPreparedensTime;
double ncDefDimTime;
double ncCloseFileTime;

void PhysicalVariableInit(){
  int lev,lat,lon;
  int plat,plon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].xp = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    global_mesh[p].zpf = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    global_mesh[p].zph = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].xp = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    mesh[p].zpf = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    mesh[p].zph = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    state[p].u = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    state[p].w = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    state[p].rho = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    state[p].pt = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    global_staticv[p].hy_dens_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    global_staticv[p].hy_dens_theta_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    global_staticv[p].hy_dens_int = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    global_staticv[p].hy_dens_theta_int = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    global_staticv[p].hy_pressure_int = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    staticv[p].hy_dens_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    staticv[p].hy_dens_theta_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    staticv[p].hy_dens_int = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    staticv[p].hy_dens_theta_int = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    staticv[p].hy_pressure_int = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    DirPara.xd = false;
    DirPara.zd = false;
    DirPara.stepnum = 0;
  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    stateout[p].dens = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    stateout[p].uwnd = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    stateout[p].wwnd = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    stateout[p].theta = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    tend[p].du = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    tend[p].dw = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    tend[p].drho = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    tend[p].dpt = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    flux[p].fu = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    flux[p].fw = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    flux[p].frho = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    flux[p].fpt = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    OutPara.step = 0;
    OutPara.interval = 2000;
  }

}

void PhysicalVariableFinish(){
  int lev,lat,lon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].xp, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].zpf, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].zph, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].xp, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].zpf, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].zph, lev, lat, lon);

  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(state[p].u, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(state[p].w, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(state[p].rho, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(state[p].pt, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(global_staticv[p].hy_dens_cell, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(global_staticv[p].hy_dens_theta_cell, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(global_staticv[p].hy_dens_int, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(global_staticv[p].hy_dens_theta_int, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(global_staticv[p].hy_pressure_int, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(staticv[p].hy_dens_cell, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 200 + 2*Proc.lev_hw;
    free_3d_array_D(staticv[p].hy_dens_theta_cell, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(staticv[p].hy_dens_int, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(staticv[p].hy_dens_theta_int, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 201 + 2*Proc.lev_hw;
    free_3d_array_D(staticv[p].hy_pressure_int, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(stateout[p].dens, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(stateout[p].uwnd, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(stateout[p].wwnd, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(stateout[p].theta, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].du, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dw, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].drho, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.full_nlev + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dpt, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    free_3d_array_D(flux[p].fu, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    free_3d_array_D(flux[p].fw, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    free_3d_array_D(flux[p].frho, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = Proc.half_nlev + 2*Proc.lev_hw;
    free_3d_array_D(flux[p].fpt, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
  }

}

void MeshInit_cp(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh){
  int i,j,k,p;

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->xp[k][j][i] = global_mesh->xp[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->zpf[k][j][i] = global_mesh->zpf[k][j][i];

  for (k = 0 ; k < Proc.half_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->zph[k][j][i] = global_mesh->zph[k][j][i];

}

void ProfilingOutPut(int id){
  char filename[30] = "Profiling_";
  char ts[10] = ".txt";
  char sid[10];

  sprintf(sid,"%d",id);
  strcat(filename,sid);
  strcat(filename,ts);
  freopen(filename, "w", stdout);
  printf("Total %.6f\n", CommTime+CompTime);
  printf("Comp %.6f\n", CompTime);
  printf("Comm %.6f\n", CommTime);
  printf("MeshInitTime %.6f\n", MeshInitTime);
  printf("StateInitialrhoTime %.6f\n", StateInitialrhoTime);
  printf("StateInitialTime %.6f\n", StateInitialTime);
  printf("ComputeTendenciesXfrhoTime %.6f\n", ComputeTendenciesXfrhoTime);
  printf("ComputeTendenciesXdrhoTime %.6f\n", ComputeTendenciesXdrhoTime);
  printf("SetBoundaryZrhoTime %.6f\n", SetBoundaryZrhoTime);
  printf("SetBoundaryZuTime %.6f\n", SetBoundaryZuTime);
  printf("SetBoundaryZwTime %.6f\n", SetBoundaryZwTime);
  printf("SetBoundaryZptTime %.6f\n", SetBoundaryZptTime);
  printf("ComputeTendenciesZfrhoTime %.6f\n", ComputeTendenciesZfrhoTime);
  printf("ComputeTendenciesZdrhoTime %.6f\n", ComputeTendenciesZdrhoTime);
  printf("UpdateStaterhoTime %.6f\n", UpdateStaterhoTime);
  printf("OutputPreparedensTime %.6f\n", OutputPreparedensTime);
  printf("ncDefDimTime %.6f\n", ncDefDimTime);
  printf("ncCloseFileTime %.6f\n", ncCloseFileTime);
  fclose(stdout);
}

