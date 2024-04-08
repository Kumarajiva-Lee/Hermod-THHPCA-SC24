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
double StateInitialh_11Time;
double StateInitialh_21Time;
double StateInitialh_31Time;
double StateInitialh_12Time;
double StateInitialh_22Time;
double StateInitialh_32Time;
double StateInitialh_13Time;
double StateInitialh_23Time;
double StateInitialh_33Time;
double updateStateh_11Time;
double updateStateh_12Time;
double updateStateh_13Time;
double updateStateh_21Time;
double updateStateh_22Time;
double updateStateh_23Time;
double updateStateh_31Time;
double updateStateh_32Time;
double updateStateh_33Time;
double updateXfh_1Time;
double updateXfh_2Time;
double updateXfu_2Time;
double updateXfv_2Time;
double updateXdxh_11Time;
double updateXdxh_12Time;
double updateXdxh_13Time;
double updateYfh_1Time;
double updateYfh_2Time;
double updateYfu_2Time;
double updateYfv_2Time;
double updateYdyh_11Time;
double update_rk1h_11Time;
double update_rk1h_12Time;
double update_rk1h_13Time;
double update_rk1h_21Time;
double update_rk1h_22Time;
double update_rk1h_23Time;
double update_rk1h_31Time;
double update_rk1h_32Time;
double update_rk1h_33Time;
double update_rk2h_11Time;
double update_rk2h_12Time;
double update_rk2h_13Time;
double update_rk2h_21Time;
double update_rk2h_22Time;
double update_rk2h_23Time;
double update_rk2h_31Time;
double update_rk2h_32Time;
double update_rk2h_33Time;
double update_rk3h_11Time;
double update_rk3h_12Time;
double update_rk3h_13Time;
double update_rk3h_21Time;
double update_rk3h_22Time;
double update_rk3h_23Time;
double update_rk3h_31Time;
double update_rk3h_32Time;
double update_rk3h_33Time;

void PhysicalVariableInit(){
  int lev,lat,lon;
  int plat,plon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].x_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].x_2 = allocate_3d_array_D(lev, lat, lon);

    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].x_3 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].y_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].y_2 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].y_3 = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].x_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].x_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].x_3 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].y_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].y_2 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].y_3 = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].h_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].u_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].v_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab5_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab6_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].jab7_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].uc_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    state[p].vc_33 = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxh_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxu_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dxv_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_11 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_21 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_31 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_12 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_22 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_32 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_13 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_23 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyh_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyu_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].dyv_33 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fh_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fh_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qh_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qh_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fu_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fu_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qu_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qu_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fv_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].fv_2 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qv_1 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    plon = Proc.full_nlon + 2*Proc.p_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    plon = Proc.full_nlat + 2*Proc.p_hw;
    lev = 1;
    tend[p].qv_2 = allocate_3d_array_D(lev, lat, lon);

  }

}

void PhysicalVariableFinish(){
  int lev,lat,lon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].x_1, lev, lat, lon);

    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].x_2, lev, lat, lon);

    lon = 600 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].x_3, lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].y_1, lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].y_2, lev, lat, lon);

    lon = 1;
    lat = 600 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].y_3, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].x_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].x_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].x_3, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].y_1, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].y_2, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].y_3, lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].h_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].u_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].v_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab5_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab6_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].jab7_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].uc_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].vc_33, lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxh_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxu_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dxv_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_11, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_21, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_31, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_12, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_22, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_32, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_13, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_23, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyh_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyu_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dyv_33, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fh_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fh_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qh_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qh_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fu_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fu_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qu_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qu_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fv_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].fv_2, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qv_1, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].qv_2, lev, lat, lon);

  }

}

void MeshInit_cp(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh){
  int i,j,k,p;

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->x_1[k][j][i] = global_mesh->x_1[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->x_2[k][j][i] = global_mesh->x_2[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->x_3[k][j][i] = global_mesh->x_3[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->y_1[k][j][i] = global_mesh->y_1[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->y_2[k][j][i] = global_mesh->y_2[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->y_3[k][j][i] = global_mesh->y_3[k][j + Proc.lat_beg][i];

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
  printf("StateInitialh_11Time %.6f\n", StateInitialh_11Time);
  printf("StateInitialh_21Time %.6f\n", StateInitialh_21Time);
  printf("StateInitialh_31Time %.6f\n", StateInitialh_31Time);
  printf("StateInitialh_12Time %.6f\n", StateInitialh_12Time);
  printf("StateInitialh_22Time %.6f\n", StateInitialh_22Time);
  printf("StateInitialh_32Time %.6f\n", StateInitialh_32Time);
  printf("StateInitialh_13Time %.6f\n", StateInitialh_13Time);
  printf("StateInitialh_23Time %.6f\n", StateInitialh_23Time);
  printf("StateInitialh_33Time %.6f\n", StateInitialh_33Time);
  printf("updateStateh_11Time %.6f\n", updateStateh_11Time);
  printf("updateStateh_12Time %.6f\n", updateStateh_12Time);
  printf("updateStateh_13Time %.6f\n", updateStateh_13Time);
  printf("updateStateh_21Time %.6f\n", updateStateh_21Time);
  printf("updateStateh_22Time %.6f\n", updateStateh_22Time);
  printf("updateStateh_23Time %.6f\n", updateStateh_23Time);
  printf("updateStateh_31Time %.6f\n", updateStateh_31Time);
  printf("updateStateh_32Time %.6f\n", updateStateh_32Time);
  printf("updateStateh_33Time %.6f\n", updateStateh_33Time);
  printf("updateXfh_1Time %.6f\n", updateXfh_1Time);
  printf("updateXfh_2Time %.6f\n", updateXfh_2Time);
  printf("updateXfu_2Time %.6f\n", updateXfu_2Time);
  printf("updateXfv_2Time %.6f\n", updateXfv_2Time);
  printf("updateXdxh_11Time %.6f\n", updateXdxh_11Time);
  printf("updateXdxh_12Time %.6f\n", updateXdxh_12Time);
  printf("updateXdxh_13Time %.6f\n", updateXdxh_13Time);
  printf("updateYfh_1Time %.6f\n", updateYfh_1Time);
  printf("updateYfh_2Time %.6f\n", updateYfh_2Time);
  printf("updateYfu_2Time %.6f\n", updateYfu_2Time);
  printf("updateYfv_2Time %.6f\n", updateYfv_2Time);
  printf("updateYdyh_11Time %.6f\n", updateYdyh_11Time);
  printf("update_rk1h_11Time %.6f\n", update_rk1h_11Time);
  printf("update_rk1h_12Time %.6f\n", update_rk1h_12Time);
  printf("update_rk1h_13Time %.6f\n", update_rk1h_13Time);
  printf("update_rk1h_21Time %.6f\n", update_rk1h_21Time);
  printf("update_rk1h_22Time %.6f\n", update_rk1h_22Time);
  printf("update_rk1h_23Time %.6f\n", update_rk1h_23Time);
  printf("update_rk1h_31Time %.6f\n", update_rk1h_31Time);
  printf("update_rk1h_32Time %.6f\n", update_rk1h_32Time);
  printf("update_rk1h_33Time %.6f\n", update_rk1h_33Time);
  printf("update_rk2h_11Time %.6f\n", update_rk2h_11Time);
  printf("update_rk2h_12Time %.6f\n", update_rk2h_12Time);
  printf("update_rk2h_13Time %.6f\n", update_rk2h_13Time);
  printf("update_rk2h_21Time %.6f\n", update_rk2h_21Time);
  printf("update_rk2h_22Time %.6f\n", update_rk2h_22Time);
  printf("update_rk2h_23Time %.6f\n", update_rk2h_23Time);
  printf("update_rk2h_31Time %.6f\n", update_rk2h_31Time);
  printf("update_rk2h_32Time %.6f\n", update_rk2h_32Time);
  printf("update_rk2h_33Time %.6f\n", update_rk2h_33Time);
  printf("update_rk3h_11Time %.6f\n", update_rk3h_11Time);
  printf("update_rk3h_12Time %.6f\n", update_rk3h_12Time);
  printf("update_rk3h_13Time %.6f\n", update_rk3h_13Time);
  printf("update_rk3h_21Time %.6f\n", update_rk3h_21Time);
  printf("update_rk3h_22Time %.6f\n", update_rk3h_22Time);
  printf("update_rk3h_23Time %.6f\n", update_rk3h_23Time);
  printf("update_rk3h_31Time %.6f\n", update_rk3h_31Time);
  printf("update_rk3h_32Time %.6f\n", update_rk3h_32Time);
  printf("update_rk3h_33Time %.6f\n", update_rk3h_33Time);
  fclose(stdout);
}

