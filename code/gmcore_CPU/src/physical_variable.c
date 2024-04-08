#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include"../LIB/Process.h"
#include"../LIB/Memory.h"

#include"physical_variable.h"

double TimeBeg ;
double TimeEnd ;
double CommTime ;
double CompTime ;
double LatLonMeshInitTime ;
double ncCloseFileTime ;
double prepareStaticdzsdlonTime ;
double prepareStaticdzsdlatTime ;
double preparegzlevgz_levTime ;
double c2auTime ;
double calc_phph_levTime ;
double calc_phph_exn_levTime ;
double calc_phphTime ;
double calc_mmTime ;
double calc_mm_levTime ;
double averageCellToLonEdgem_lonTime ;
double averageCellToLatEdgem_latTime ;
double interpCellToVtxm_vtxTime ;
double calc_ttTime ;
double accum_uv_celluuTime ;
double accum_uv_cellvvTime ;
double accum_uv_cellu0Time ;
double accum_uv_cellv0Time ;
double accum_uv_cellcflxTime ;
double accum_uv_cellcflyTime ;
double accum_uv_celldivxTime ;
double accum_uv_celltmpsumTime ;
double accum_uv_celldivyTime ;
double calc_mfmfx_lonTime ;
double calc_mfmfy_latTime ;
double accum_mf_cellmfxTime ;
double accum_mf_cellmfyTime ;
double calc_mfmfx_latTime ;
double calc_mfmfy_lonTime ;
double calc_kekeTime ;
double calc_ketmpsumTime ;
double calc_vorvorTime ;
double calc_vortmpsumTime ;
double calc_pvpvTime ;
double interp_pv_upwindpv_lonTime ;
double interp_pv_upwindpv_latTime ;
double calc_divdivTime ;
double calc_divtmpsumTime ;
double calc_gz_levgz_levTime ;
double copy_old_mold_mTime ;
double Filter_InitTime ;
double calc_grad_mfdmfdlonTime ;
double calc_grad_mfdmfdlatTime ;
double calc_grad_mftmpsumTime ;
double calc_dphsdtdphsTime ;
double calc_we_levwe_levTime ;
double accum_we_levweTime ;
double accum_we_levcflzTime ;
double interp_lev_edge_to_lev_lon_edgewe_lev_lonTime ;
double interp_lev_edge_to_lev_lat_edgewe_lev_latTime ;
double calc_wedudlev_wedvdlevwedudlevTime ;
double calc_wedudlev_wedvdlevwedvdlevTime ;
double hflx_ppm_innerqlxTime ;
double hflx_ppm_innerptf_lonTime ;
double hflx_ppm_innerptf_latTime ;
double ffsl_calc_tracer_hflxqxTime ;
double ffsl_calc_tracer_hflxtmpsumTime ;
double hflx_ppm_outerqlxTime ;
double hflx_ppm_outerptf_lonTime ;
double hflx_ppm_outerptf_latTime ;
double calc_grad_ptfdptfdlonTime ;
double calc_grad_ptfdptfdlatTime ;
double calc_grad_ptftmpsumTime ;
double calc_grad_ptfptTime ;
double vflx_ppmqlxTime ;
double vflx_ppmptf_levTime ;
double calc_grad_ptfdptfdlevTime ;
double calc_coriolisqhuTime ;
double calc_coriolisqhvTime ;
double calc_grad_kedkedlonTime ;
double calc_grad_kedkedlatTime ;
double calc_tend_forwardduTime ;
double calc_tend_forwarddvTime ;
double calc_tend_forwarddptTime ;
double Filter_On_CellTime ;
double update_statephsTime ;
double update_stateptTime ;
double update_stategzTime ;
double Filter_on_lon_edgeTime ;
double Filter_on_lat_edgeTime ;
double update_stateu_lonTime ;
double update_statev_latTime ;
double pgf_lin97pgf_lonTime ;
double pgf_lin97pgf_latTime ;
double calc_tend_backwardduTime ;
double calc_tend_backwarddvTime ;
double div_damp_runu_lonTime ;
double div_damp_runv_latTime ;
double smag_damp_runsmag_tTime ;
double smag_damp_runsmag_sTime ;
double smag_damp_runkmh_lonTime ;
double smag_damp_runkmh_latTime ;
double smag_damp_runsmag_dudtTime ;
double smag_damp_runu_lonTime ;
double smag_damp_runsmag_dvdtTime ;
double smag_damp_runv_latTime ;
double trickyptptTime ;
double pole_damp_runptTime ;
double Time_AdvanceTime ;
double hflx_ppm_innerqmf_lonTime ;
double hflx_ppm_innerqmf_latTime ;
double hflx_ppm_outerqmf_lonTime ;
double hflx_ppm_outerqmf_latTime ;
double adv_runqvTime ;
double vflx_ppmqmf_levTime ;
double DiagnoseTime ;
double adv_runtmpsumTime;

void PhysicalVariableInit(){
  int lev,lat,lon;
  int plat,plon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].dlat = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].full_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].half_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    global_mesh[p].full_lev = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    global_mesh[p].half_lev = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].full_cos_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].half_cos_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].full_sin_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].half_sin_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_cos_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_cos_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_sin_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_sin_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].full_lon_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    global_mesh[p].half_lon_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_lat_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_lat_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lon_west = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lon_east = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lon_north = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lon_south = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lat_west = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lat_east = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lat_north = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_lat_south = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_vtx = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_subcell_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].area_subcell_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].de_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].de_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].le_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].le_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_f = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_f = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_tangent_wgt_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].full_tangent_wgt_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_tangent_wgt_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    global_mesh[p].half_tangent_wgt_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    global_mesh[p].hyai = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    global_mesh[p].hybi = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    global_mesh[p].hyam = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    global_mesh[p].hybm = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    global_mesh[p].c_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    global_mesh[p].c_lat = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].dlat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].full_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].half_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    mesh[p].full_lev = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    mesh[p].half_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].full_cos_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].half_cos_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].full_sin_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].half_sin_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_cos_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_cos_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_sin_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_sin_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].full_lon_deg = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    mesh[p].half_lon_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_lat_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_lat_deg = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_cell = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lon_west = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lon_east = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lon_north = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lon_south = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lat_west = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lat_east = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lat_north = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_lat_south = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_vtx = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_subcell_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].area_subcell_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].de_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].de_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].le_lat = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].le_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_f = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_f = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_tangent_wgt_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].full_tangent_wgt_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_tangent_wgt_0 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    mesh[p].half_tangent_wgt_1 = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    mesh[p].hyai = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    mesh[p].hybi = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    mesh[p].hyam = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    mesh[p].hybm = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    mesh[p].c_lon = allocate_3d_array_D(lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    mesh[p].c_lat = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].u = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].u_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].u_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].v = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].v_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].v_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].we_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].we_lev_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].we_lev_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].gz = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].gz_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].m = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].m_vtx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].m_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].m_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].m_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].mfx_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].mfy_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].mfx_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].mfy_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].pv = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].pv_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].pv_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].ke = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].pt = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].ptf_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].ptf_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].ptf_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].t = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].ph = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].ph_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    state[p].ph_exn_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    state[p].phs = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].div = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].vor = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].qv = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].qm = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].smag_t = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].smag_s = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].kmh = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].kmh_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].kmh_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].q = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    state[p].tmpsum = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    staticv[p].gzs = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    staticv[p].dzsdlon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    staticv[p].dzsdlat = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qmf_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qmf_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].qmf_lev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].old_m = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].mfx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].mfy = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].mm = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].m0 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].uu = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].u0 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].vv = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].v0 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].we = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].we0 = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].cflx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].cfly = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    adv[p].cflz = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].divx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].divy = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qlx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qly = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].dqx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].dqy = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].q6x = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].q6y = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qx = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    adv[p].qy = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    advptPara.dynamic = true;
    advptPara.nstep = 1;
    advptPara.uv_step = 0;
    advptPara.we_step = 0;
    advptPara.mf_step = 0;
  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].du = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dv = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dgz = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dpt = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    tend[p].dphs = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].qhv = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].qhu = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dkedlon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dkedlat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dmfdlon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dmfdlat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dptfdlon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dptfdlat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].dptfdlev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].pgf_lon = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].pgf_lat = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].wedudlev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].wedvdlev = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].smag_dptdt = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].smag_dudt = allocate_3d_array_D(lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    tend[p].smag_dvdt = allocate_3d_array_D(lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    tendPara.phs = false;
    tendPara.pt = false;
    tendPara.gz = false;
    tendPara.u = false;
    tendPara.v = false;
  }

  for (p = 0 ; p < 1 ; p++){
    advmPara.dynamic = false;
    advmPara.nstep = 1;
    advmPara.uv_step = 0;
    advmPara.we_step = 0;
    advmPara.mf_step = 0;
  }

}

void PhysicalVariableFinish(){
  int lev,lat,lon;
  int i,j,k,p;

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].dlat, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_lon, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_lon, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_lat, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_lat, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].full_lev, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].half_lev, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_cos_lon, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_cos_lon, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_sin_lon, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_sin_lon, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_cos_lat, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_cos_lat, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_sin_lat, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_sin_lat, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_lon_deg, lev, lat, lon);

    lon = 2400 + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_lon_deg, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_lat_deg, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_lat_deg, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_cell, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lon, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lon_west, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lon_east, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lon_north, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lon_south, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lat, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lat_west, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lat_east, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lat_north, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_lat_south, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_vtx, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_subcell_0, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].area_subcell_1, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].de_lon, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].de_lat, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].le_lat, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].le_lon, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_f, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_f, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_tangent_wgt_0, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].full_tangent_wgt_1, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_tangent_wgt_0, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(global_mesh[p].half_tangent_wgt_1, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].hyai, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].hybi, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].hyam, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].hybm, lev, lat, lon);

    lon = 1;
    lat = 1200 + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].c_lon, lev, lat, lon);

    lon = 1;
    lat = 1199 + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(global_mesh[p].c_lat, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].dlat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].full_lon, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].half_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_lat, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].full_lev, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].half_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].full_cos_lon, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].half_cos_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].full_sin_lon, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].half_sin_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_cos_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_cos_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_sin_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_sin_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].full_lon_deg, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = 1;
    lev = 1;
    free_3d_array_D(mesh[p].half_lon_deg, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_lat_deg, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_lat_deg, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_cell, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lon_west, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lon_east, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lon_north, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lon_south, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lat_west, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lat_east, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lat_north, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_lat_south, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_vtx, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_subcell_0, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].area_subcell_1, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].de_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].de_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].le_lat, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].le_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_f, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_f, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_tangent_wgt_0, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].full_tangent_wgt_1, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_tangent_wgt_0, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(mesh[p].half_tangent_wgt_1, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].hyai, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].hybi, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].hyam, lev, lat, lon);

    lon = 1;
    lat = 1;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].hybm, lev, lat, lon);

    lon = 1;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].c_lon, lev, lat, lon);

    lon = 1;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(mesh[p].c_lat, lev, lat, lon);

  }

  for (p = 0 ; p < 3 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].u, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].u_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].u_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].v, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].v_lat, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].v_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].we_lev, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].we_lev_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].we_lev_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].gz, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].gz_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].m, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].m_vtx, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].m_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].m_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].m_lev, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].mfx_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].mfy_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].mfx_lat, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].mfy_lon, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].pv, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].pv_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].pv_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ke, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].pt, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ptf_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ptf_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ptf_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].t, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ph, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ph_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].ph_exn_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(state[p].phs, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].div, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].vor, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].qv, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].qm, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].smag_t, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].smag_s, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].kmh, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].kmh_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].kmh_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].q, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(state[p].tmpsum, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(staticv[p].gzs, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(staticv[p].dzsdlon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(staticv[p].dzsdlat, lev, lat, lon);

  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qmf_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qmf_lat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qmf_lev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].old_m, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].mfx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].mfy, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].mm, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].m0, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].uu, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].u0, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].vv, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].v0, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].we, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].we0, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].cflx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].cfly, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 33 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].cflz, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].divx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].divy, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qlx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qly, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].dqx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].dqy, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].q6x, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].q6y, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qx, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(adv[p].qy, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
  }

  for (p = 0 ; p < 2 ; p++){
    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].du, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dv, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dgz, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dpt, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 1;
    free_3d_array_D(tend[p].dphs, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].qhv, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].qhu, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dkedlon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dkedlat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dmfdlon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dmfdlat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dptfdlon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dptfdlat, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].dptfdlev, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].pgf_lon, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].pgf_lat, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].wedudlev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].wedvdlev, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].smag_dptdt, lev, lat, lon);

    lon = Proc.half_nlon + 2*Proc.lon_hw;
    lat = Proc.full_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].smag_dudt, lev, lat, lon);

    lon = Proc.full_nlon + 2*Proc.lon_hw;
    lat = Proc.half_nlat + 2*Proc.lat_hw;
    lev = 32 + 2*Proc.lev_hw;
    free_3d_array_D(tend[p].smag_dvdt, lev, lat, lon);

  }

  for (p = 0 ; p < 1 ; p++){
  }

  for (p = 0 ; p < 1 ; p++){
  }

}

void LatLonMeshInit_cp(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh){
  int i,j,k,p;

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->dlat[k][j][i] = global_mesh->dlat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->full_lon[k][j][i] = global_mesh->full_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.half_nlon + 2*Proc.lon_hw ; i ++)
        mesh->half_lon[k][j][i] = global_mesh->half_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_lat[k][j][i] = global_mesh->full_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_lat[k][j][i] = global_mesh->half_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_lev[k][j][i] = global_mesh->full_lev[k][j][i];

  for (k = 0 ; k < Proc.half_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_lev[k][j][i] = global_mesh->half_lev[k][j][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->full_cos_lon[k][j][i] = global_mesh->full_cos_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.half_nlon + 2*Proc.lon_hw ; i ++)
        mesh->half_cos_lon[k][j][i] = global_mesh->half_cos_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->full_sin_lon[k][j][i] = global_mesh->full_sin_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.half_nlon + 2*Proc.lon_hw ; i ++)
        mesh->half_sin_lon[k][j][i] = global_mesh->half_sin_lon[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_cos_lat[k][j][i] = global_mesh->full_cos_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_cos_lat[k][j][i] = global_mesh->half_cos_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_sin_lat[k][j][i] = global_mesh->full_sin_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_sin_lat[k][j][i] = global_mesh->half_sin_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.full_nlon + 2*Proc.lon_hw ; i ++)
        mesh->full_lon_deg[k][j][i] = global_mesh->full_lon_deg[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < Proc.half_nlon + 2*Proc.lon_hw ; i ++)
        mesh->half_lon_deg[k][j][i] = global_mesh->half_lon_deg[k][j][i + Proc.lon_beg];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_lat_deg[k][j][i] = global_mesh->full_lat_deg[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_lat_deg[k][j][i] = global_mesh->half_lat_deg[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_cell[k][j][i] = global_mesh->area_cell[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lon[k][j][i] = global_mesh->area_lon[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lon_west[k][j][i] = global_mesh->area_lon_west[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lon_east[k][j][i] = global_mesh->area_lon_east[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lon_north[k][j][i] = global_mesh->area_lon_north[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lon_south[k][j][i] = global_mesh->area_lon_south[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lat[k][j][i] = global_mesh->area_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lat_west[k][j][i] = global_mesh->area_lat_west[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lat_east[k][j][i] = global_mesh->area_lat_east[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lat_north[k][j][i] = global_mesh->area_lat_north[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_lat_south[k][j][i] = global_mesh->area_lat_south[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_vtx[k][j][i] = global_mesh->area_vtx[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_subcell_0[k][j][i] = global_mesh->area_subcell_0[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->area_subcell_1[k][j][i] = global_mesh->area_subcell_1[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->de_lon[k][j][i] = global_mesh->de_lon[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->de_lat[k][j][i] = global_mesh->de_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->le_lat[k][j][i] = global_mesh->le_lat[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->le_lon[k][j][i] = global_mesh->le_lon[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_f[k][j][i] = global_mesh->full_f[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_f[k][j][i] = global_mesh->half_f[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_tangent_wgt_0[k][j][i] = global_mesh->full_tangent_wgt_0[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->full_tangent_wgt_1[k][j][i] = global_mesh->full_tangent_wgt_1[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_tangent_wgt_0[k][j][i] = global_mesh->half_tangent_wgt_0[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < 1 ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->half_tangent_wgt_1[k][j][i] = global_mesh->half_tangent_wgt_1[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < Proc.half_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->hyai[k][j][i] = global_mesh->hyai[k][j][i];

  for (k = 0 ; k < Proc.half_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->hybi[k][j][i] = global_mesh->hybi[k][j][i];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->hyam[k][j][i] = global_mesh->hyam[k][j][i];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < 1 ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->hybm[k][j][i] = global_mesh->hybm[k][j][i];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < Proc.full_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->c_lon[k][j][i] = global_mesh->c_lon[k][j + Proc.lat_beg][i];

  for (k = 0 ; k < Proc.full_nlev + 2*Proc.lev_hw ; k ++)
    for (j = 0 ; j < Proc.half_nlat + 2*Proc.lat_hw ; j ++)
      for (i = 0 ; i < 1 ; i ++)
        mesh->c_lat[k][j][i] = global_mesh->c_lat[k][j + Proc.lat_beg][i];

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
  printf("LatLonMeshInitTime %.6f\n", LatLonMeshInitTime);
  printf("ncCloseFileTime %.6f\n", ncCloseFileTime);
  printf("prepareStaticdzsdlonTime %.6f\n", prepareStaticdzsdlonTime);
  printf("prepareStaticdzsdlatTime %.6f\n", prepareStaticdzsdlatTime);
  printf("preparegzlevgz_levTime %.6f\n", preparegzlevgz_levTime);
  printf("c2auTime %.6f\n", c2auTime);
  printf("calc_phph_levTime %.6f\n", calc_phph_levTime);
  printf("calc_phph_exn_levTime %.6f\n", calc_phph_exn_levTime);
  printf("calc_phphTime %.6f\n", calc_phphTime);
  printf("calc_mmTime %.6f\n", calc_mmTime);
  printf("calc_mm_levTime %.6f\n", calc_mm_levTime);
  printf("averageCellToLonEdgem_lonTime %.6f\n", averageCellToLonEdgem_lonTime);
  printf("averageCellToLatEdgem_latTime %.6f\n", averageCellToLatEdgem_latTime);
  printf("interpCellToVtxm_vtxTime %.6f\n", interpCellToVtxm_vtxTime);
  printf("calc_ttTime %.6f\n", calc_ttTime);
  printf("accum_uv_celluuTime %.6f\n", accum_uv_celluuTime);
  printf("accum_uv_cellvvTime %.6f\n", accum_uv_cellvvTime);
  printf("accum_uv_cellu0Time %.6f\n", accum_uv_cellu0Time);
  printf("accum_uv_cellv0Time %.6f\n", accum_uv_cellv0Time);
  printf("accum_uv_cellcflxTime %.6f\n", accum_uv_cellcflxTime);
  printf("accum_uv_cellcflyTime %.6f\n", accum_uv_cellcflyTime);
  printf("accum_uv_celldivxTime %.6f\n", accum_uv_celldivxTime);
  printf("accum_uv_celltmpsumTime %.6f\n", accum_uv_celltmpsumTime);
  printf("accum_uv_celldivyTime %.6f\n", accum_uv_celldivyTime);
  printf("calc_mfmfx_lonTime %.6f\n", calc_mfmfx_lonTime);
  printf("calc_mfmfy_latTime %.6f\n", calc_mfmfy_latTime);
  printf("accum_mf_cellmfxTime %.6f\n", accum_mf_cellmfxTime);
  printf("accum_mf_cellmfyTime %.6f\n", accum_mf_cellmfyTime);
  printf("calc_mfmfx_latTime %.6f\n", calc_mfmfx_latTime);
  printf("calc_mfmfy_lonTime %.6f\n", calc_mfmfy_lonTime);
  printf("calc_kekeTime %.6f\n", calc_kekeTime);
  printf("calc_ketmpsumTime %.6f\n", calc_ketmpsumTime);
  printf("calc_vorvorTime %.6f\n", calc_vorvorTime);
  printf("calc_vortmpsumTime %.6f\n", calc_vortmpsumTime);
  printf("calc_pvpvTime %.6f\n", calc_pvpvTime);
  printf("interp_pv_upwindpv_lonTime %.6f\n", interp_pv_upwindpv_lonTime);
  printf("interp_pv_upwindpv_latTime %.6f\n", interp_pv_upwindpv_latTime);
  printf("calc_divdivTime %.6f\n", calc_divdivTime);
  printf("calc_divtmpsumTime %.6f\n", calc_divtmpsumTime);
  printf("calc_gz_levgz_levTime %.6f\n", calc_gz_levgz_levTime);
  printf("copy_old_mold_mTime %.6f\n", copy_old_mold_mTime);
  printf("Filter_InitTime %.6f\n", Filter_InitTime);
  printf("calc_grad_mfdmfdlonTime %.6f\n", calc_grad_mfdmfdlonTime);
  printf("calc_grad_mfdmfdlatTime %.6f\n", calc_grad_mfdmfdlatTime);
  printf("calc_grad_mftmpsumTime %.6f\n", calc_grad_mftmpsumTime);
  printf("calc_dphsdtdphsTime %.6f\n", calc_dphsdtdphsTime);
  printf("calc_we_levwe_levTime %.6f\n", calc_we_levwe_levTime);
  printf("accum_we_levweTime %.6f\n", accum_we_levweTime);
  printf("accum_we_levcflzTime %.6f\n", accum_we_levcflzTime);
  printf("interp_lev_edge_to_lev_lon_edgewe_lev_lonTime %.6f\n", interp_lev_edge_to_lev_lon_edgewe_lev_lonTime);
  printf("interp_lev_edge_to_lev_lat_edgewe_lev_latTime %.6f\n", interp_lev_edge_to_lev_lat_edgewe_lev_latTime);
  printf("calc_wedudlev_wedvdlevwedudlevTime %.6f\n", calc_wedudlev_wedvdlevwedudlevTime);
  printf("calc_wedudlev_wedvdlevwedvdlevTime %.6f\n", calc_wedudlev_wedvdlevwedvdlevTime);
  printf("hflx_ppm_innerqlxTime %.6f\n", hflx_ppm_innerqlxTime);
  printf("hflx_ppm_innerptf_lonTime %.6f\n", hflx_ppm_innerptf_lonTime);
  printf("hflx_ppm_innerptf_latTime %.6f\n", hflx_ppm_innerptf_latTime);
  printf("ffsl_calc_tracer_hflxqxTime %.6f\n", ffsl_calc_tracer_hflxqxTime);
  printf("ffsl_calc_tracer_hflxtmpsumTime %.6f\n", ffsl_calc_tracer_hflxtmpsumTime);
  printf("hflx_ppm_outerqlxTime %.6f\n", hflx_ppm_outerqlxTime);
  printf("hflx_ppm_outerptf_lonTime %.6f\n", hflx_ppm_outerptf_lonTime);
  printf("hflx_ppm_outerptf_latTime %.6f\n", hflx_ppm_outerptf_latTime);
  printf("calc_grad_ptfdptfdlonTime %.6f\n", calc_grad_ptfdptfdlonTime);
  printf("calc_grad_ptfdptfdlatTime %.6f\n", calc_grad_ptfdptfdlatTime);
  printf("calc_grad_ptftmpsumTime %.6f\n", calc_grad_ptftmpsumTime);
  printf("calc_grad_ptfptTime %.6f\n", calc_grad_ptfptTime);
  printf("vflx_ppmqlxTime %.6f\n", vflx_ppmqlxTime);
  printf("vflx_ppmptf_levTime %.6f\n", vflx_ppmptf_levTime);
  printf("calc_grad_ptfdptfdlevTime %.6f\n", calc_grad_ptfdptfdlevTime);
  printf("calc_coriolisqhuTime %.6f\n", calc_coriolisqhuTime);
  printf("calc_coriolisqhvTime %.6f\n", calc_coriolisqhvTime);
  printf("calc_grad_kedkedlonTime %.6f\n", calc_grad_kedkedlonTime);
  printf("calc_grad_kedkedlatTime %.6f\n", calc_grad_kedkedlatTime);
  printf("calc_tend_forwardduTime %.6f\n", calc_tend_forwardduTime);
  printf("calc_tend_forwarddvTime %.6f\n", calc_tend_forwarddvTime);
  printf("calc_tend_forwarddptTime %.6f\n", calc_tend_forwarddptTime);
  printf("Filter_On_CellTime %.6f\n", Filter_On_CellTime);
  printf("update_statephsTime %.6f\n", update_statephsTime);
  printf("update_stateptTime %.6f\n", update_stateptTime);
  printf("update_stategzTime %.6f\n", update_stategzTime);
  printf("Filter_on_lon_edgeTime %.6f\n", Filter_on_lon_edgeTime);
  printf("Filter_on_lat_edgeTime %.6f\n", Filter_on_lat_edgeTime);
  printf("update_stateu_lonTime %.6f\n", update_stateu_lonTime);
  printf("update_statev_latTime %.6f\n", update_statev_latTime);
  printf("pgf_lin97pgf_lonTime %.6f\n", pgf_lin97pgf_lonTime);
  printf("pgf_lin97pgf_latTime %.6f\n", pgf_lin97pgf_latTime);
  printf("calc_tend_backwardduTime %.6f\n", calc_tend_backwardduTime);
  printf("calc_tend_backwarddvTime %.6f\n", calc_tend_backwarddvTime);
  printf("div_damp_runu_lonTime %.6f\n", div_damp_runu_lonTime);
  printf("div_damp_runv_latTime %.6f\n", div_damp_runv_latTime);
  printf("smag_damp_runsmag_tTime %.6f\n", smag_damp_runsmag_tTime);
  printf("smag_damp_runsmag_sTime %.6f\n", smag_damp_runsmag_sTime);
  printf("smag_damp_runkmh_lonTime %.6f\n", smag_damp_runkmh_lonTime);
  printf("smag_damp_runkmh_latTime %.6f\n", smag_damp_runkmh_latTime);
  printf("smag_damp_runsmag_dudtTime %.6f\n", smag_damp_runsmag_dudtTime);
  printf("smag_damp_runu_lonTime %.6f\n", smag_damp_runu_lonTime);
  printf("smag_damp_runsmag_dvdtTime %.6f\n", smag_damp_runsmag_dvdtTime);
  printf("smag_damp_runv_latTime %.6f\n", smag_damp_runv_latTime);
  printf("trickyptptTime %.6f\n", trickyptptTime);
  printf("pole_damp_runptTime %.6f\n", pole_damp_runptTime);
  printf("Time_AdvanceTime %.6f\n", Time_AdvanceTime);
  printf("hflx_ppm_innerqmf_lonTime %.6f\n", hflx_ppm_innerqmf_lonTime);
  printf("hflx_ppm_innerqmf_latTime %.6f\n", hflx_ppm_innerqmf_latTime);
  printf("hflx_ppm_outerqmf_lonTime %.6f\n", hflx_ppm_outerqmf_lonTime);
  printf("hflx_ppm_outerqmf_latTime %.6f\n", hflx_ppm_outerqmf_latTime);
  printf("adv_runqvTime %.6f\n", adv_runqvTime);
  printf("adv_runtmpsumTime %.6f\n", adv_runtmpsumTime);
  printf("vflx_ppmqmf_levTime %.6f\n", vflx_ppmqmf_levTime);
  printf("DiagnoseTime %.6f\n", DiagnoseTime);
  fclose(stdout);
}

