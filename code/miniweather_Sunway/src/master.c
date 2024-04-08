#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <crts.h>


#include"../LIB/Process.h"
#include"../LIB/Communicator.h"
#include"../LIB/ParaParam.h"
#include"../LIB/Async.h"
#include"../LIB/NCop.h"
#include"../LIB/Time.h"
#include"../LIB/Memory.h"
#include"../LIB/Diagnose.h"

#include"namelist.h"
#include"physical_variable.h"

#include "slave_struct.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

extern SLAVE_FUN(global_prepare)(void*);
extern SLAVE_FUN(StateInitial_0_0)(void*);
extern SLAVE_FUN(StateInitial_0_1)(void*);
extern SLAVE_FUN(ComputeTendenciesX_0_20)(void*);
extern SLAVE_FUN(ComputeTendenciesX_0_21)(void*);
extern SLAVE_FUN(ComputeTendenciesZ_0_20)(void*);
extern SLAVE_FUN(ComputeTendenciesZ_0_21)(void*);
extern SLAVE_FUN(ComputeTendenciesZ_0_22)(void*);
extern SLAVE_FUN(ComputeTendenciesZ_0_23)(void*);
extern SLAVE_FUN(ComputeTendenciesX_0_22)(void*);
extern SLAVE_FUN(ComputeTendenciesX_0_23)(void*);
extern SLAVE_FUN(UpdateState_0_5)(void*);
extern SLAVE_FUN(OutputPrepare_0_0)(void*);

struct MeshVector3 disect(int inx, struct MeshVector6* loopinfo)
{
  int k, tmp, j, i;
  k = (inx % (loopinfo->lz));
  tmp = (inx / (loopinfo->lz));
  j = (tmp % (loopinfo->ly));
  tmp = (tmp / (loopinfo->ly));
  i = (tmp % (loopinfo->lx));
  return (struct MeshVector3){(i + (loopinfo->sx)),(j + (loopinfo->sy)),(k + (loopinfo->sz))};
}

double sample_ellipse_cosine(double x, double z, double amp, double x0, double z0, double xrad, double zrad)
{
  double dist;
  dist = ((sqrt(((((x - x0) / xrad) * ((x - x0) / xrad)) + (((z - z0) / zrad) * ((z - z0) / zrad)))) * 3.141592653589793) / 2.0);
  if ((dist <= (3.141592653589793 / 2.0)))
  {
    return (amp * pow(cos(dist),2));
  }
  else
  {
    return 0.0;
  }
}

struct Vector2 hydro_const_theta(double z)
{
  double theta0, exner0, t, exner, p, rt, r;
  theta0 = 300.0;
  exner0 = 1.0;
  t = theta0;
  exner = (exner0 - ((9.8 * z) / (1004.0 * theta0)));
  p = (100000.0 * pow(exner,(1004.0 / 287.0)));
  rt = pow((p / 27.562941092972594),(1.0 / 1.400278940027894));
  r = (rt / t);
  return (struct Vector2){r,t};
}

struct InitVector collision(double x, double z)
{
  struct Vector2 buf;
  struct InitVector ret;
  buf = (struct Vector2){0.0,0.0};
  buf = hydro_const_theta(z);
  ret = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
  ret.r = 0.0;
  ret.t = 0.0;
  ret.u = 0.0;
  ret.w = 0.0;
  ret.t = (ret.t + sample_ellipse_cosine(x,z,20.0,(20000.0 / 2),2000.0,2000.0,2000.0));
  ret.t = (ret.t + sample_ellipse_cosine(x,z,-20.0,(20000.0 / 2),8000.0,2000.0,2000.0));
  ret.hr = buf.x;
  ret.ht = buf.y;
  return ret;
}

struct InitVector thermal(double x, double z)
{
  struct Vector2 buf;
  struct InitVector ret;
  buf = (struct Vector2){0.0,0.0};
  buf = hydro_const_theta(z);
  ret = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
  ret.r = 0.0;
  ret.u = 0.0;
  ret.w = 0.0;
  ret.t = (0.0 + sample_ellipse_cosine(x,z,3.0,(20000.0 / 2),2000.0,2000.0,2000.0));
  ret.hr = buf.x;
  ret.ht = buf.y;
  return ret;
}

struct Vector2 hydro_const_bvfreq(double z, double bv_freq0)
{
  double theta0, exner0, t, exner, p, rt, r;
  theta0 = 300.0;
  exner0 = 1.0;
  t = (theta0 * exp((((bv_freq0 * bv_freq0) / 9.8) * z)));
  exner = (exner0 - ((((9.8 * 9.8) / ((1004.0 * bv_freq0) * bv_freq0)) * (t - theta0)) / (t * theta0)));
  p = (100000.0 * pow(exner,(1004.0 / 287.0)));
  rt = pow((p / 27.562941092972594),(1.0 / 1.400278940027894));
  r = (rt / t);
  return (struct Vector2){r,rt};
}

struct InitVector gravity_waves(double x, double z)
{
  struct Vector2 buf;
  double bv_freq0;
  struct InitVector ret;
  buf = (struct Vector2){0.0,0.0};
  bv_freq0 = 0.02;
  buf = hydro_const_bvfreq(z,bv_freq0);
  ret = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
  ret.r = 0.0;
  ret.t = 0.0;
  ret.u = 15.0;
  ret.w = 0.0;
  ret.hr = buf.x;
  ret.ht = buf.y;
  return ret;
}

void MeshInit(struct HybridMeshField* mesh)
{
  int t, i, k;
  struct MeshVector6 loopinfo;
  struct MeshVector3 v_inx;
  t = 0;
  loopinfo = (struct MeshVector6){1,1,400,0,0,-2};
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=0; i<1+2*Proc.lon_hw; i+=1){
    v_inx = disect(t,&loopinfo);
    t = (t + 1);
    mesh->xp[0][0][i] = ((v_inx.z + 0.5) * 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=1+Proc.lon_hw; i<400+Proc.lon_hw; i+=1){
    mesh->xp[0][0][i] = (mesh->xp[0][0][(i - 1)] + 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  t = 0;
  loopinfo = (struct MeshVector6){200,1,1,-2,0,0};
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    v_inx = disect(t,&loopinfo);
    t = (t + 1);
    mesh->zpf[k][0][0] = ((v_inx.x + 0.5) * 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
    mesh->zpf[k][0][0] = (mesh->zpf[(k - 1)][0][0] + 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  t = 0;
  loopinfo = (struct MeshVector6){200,1,1,-2,0,0};
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    v_inx = disect(t,&loopinfo);
    t = (t + 1);
    mesh->zph[k][0][0] = (v_inx.x * 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
    mesh->zph[k][0][0] = (mesh->zph[(k - 1)][0][0] + 50.0);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void spaceOperatorInit_0(struct HybridStateField* state, struct HybridStateField* state_tmp, struct HybridStaticField* staticv, struct HybridMeshField* mesh, struct HybridDirPara* DirPara)
{
  int kk, ii, i, k, j;
  double qpi, qwi, qpk, qwk, x, z;
  struct InitVector InitV;
  DirPara = DirPara;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_rho]);
  HaloWait(Proc,&Proc.FieldReq[async_u]);
  HaloWait(Proc,&Proc.FieldReq[async_w]);
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  StateInitial_0_0_info StateInitial_0_0_para;
  StateInitial_0_0_para.lz = 200 - 0;
  StateInitial_0_0_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
  StateInitial_0_0_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
  StateInitial_0_0_para.oz = 0;
  StateInitial_0_0_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  StateInitial_0_0_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  StateInitial_0_0_para.hx = Proc.lon_hw;
  StateInitial_0_0_para.hy = Proc.lat_hw;
  StateInitial_0_0_para.hz = 0;
  StateInitial_0_0_para.bx = 256;
  StateInitial_0_0_para.by = 8;
  StateInitial_0_0_para.bz = 1;
  StateInitial_0_0_para.mx = 1;
  StateInitial_0_0_para.my = 2;
  StateInitial_0_0_para.mz = 32;
  StateInitial_0_0_para.staterho = &state->rho[Proc.lev_hw - 0+StateInitial_0_0_para.oz][Proc.lat_hw -StateInitial_0_0_para.hy +StateInitial_0_0_para.oy][Proc.lon_hw -StateInitial_0_0_para.hx +StateInitial_0_0_para.ox];
  StateInitial_0_0_para.stateu = &state->u[Proc.lev_hw - 0+StateInitial_0_0_para.oz][Proc.lat_hw -StateInitial_0_0_para.hy +StateInitial_0_0_para.oy][Proc.lon_hw -StateInitial_0_0_para.hx +StateInitial_0_0_para.ox];
  StateInitial_0_0_para.statew = &state->w[Proc.lev_hw - 0+StateInitial_0_0_para.oz][Proc.lat_hw -StateInitial_0_0_para.hy +StateInitial_0_0_para.oy][Proc.lon_hw -StateInitial_0_0_para.hx +StateInitial_0_0_para.ox];
  StateInitial_0_0_para.statept = &state->pt[Proc.lev_hw - 0+StateInitial_0_0_para.oz][Proc.lat_hw -StateInitial_0_0_para.hy +StateInitial_0_0_para.oy][Proc.lon_hw -StateInitial_0_0_para.hx +StateInitial_0_0_para.ox];
  StateInitial_0_0_para.xp = &mesh->xp[0][0][0];
  StateInitial_0_0_para.zpf = &mesh->zpf[0][0][0];
  
  athread_spawn(StateInitial_0_0, &StateInitial_0_0_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialrhoTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_rho]);
  HaloWait(Proc,&Proc.FieldReq[async_u]);
  HaloWait(Proc,&Proc.FieldReq[async_w]);
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  StateInitial_0_1_info StateInitial_0_1_para;
  StateInitial_0_1_para.lz = 200 - 0;
  StateInitial_0_1_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
  StateInitial_0_1_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
  StateInitial_0_1_para.oz = 0;
  StateInitial_0_1_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  StateInitial_0_1_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  StateInitial_0_1_para.hx = 2;
  StateInitial_0_1_para.hy = 0;
  StateInitial_0_1_para.hz = 0;
  StateInitial_0_1_para.bx = 256;
  StateInitial_0_1_para.by = 8;
  StateInitial_0_1_para.bz = 1;
  StateInitial_0_1_para.mx = 1;
  StateInitial_0_1_para.my = 2;
  StateInitial_0_1_para.mz = 32;
  StateInitial_0_1_para.staterho = &state->rho[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.stateu = &state->u[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.statew = &state->w[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.statept = &state->pt[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.state_tmprho = &state_tmp->rho[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.state_tmpu = &state_tmp->u[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.state_tmpw = &state_tmp->w[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  StateInitial_0_1_para.state_tmppt = &state_tmp->pt[Proc.lev_hw - 0+StateInitial_0_1_para.oz][Proc.lat_hw -0+StateInitial_0_1_para.oy][Proc.lon_hw -2+StateInitial_0_1_para.ox];
  
  athread_spawn(StateInitial_0_1, &StateInitial_0_1_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialrhoTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_tmp->rho[0][0][0], &Proc.FieldReq[async_rho], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_tmp->u[0][0][0], &Proc.FieldReq[async_u], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_tmp->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_tmp->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
    for(kk=0; kk<3; kk+=1){
      if ((kk == 0))
      {
        qpk = 0.11270166537925831;
        qwk = 0.2777777777777778;
      }
      else
      {
        if ((kk == 1))
        {
          qpk = 0.5;
          qwk = 0.4444444444444444;
        }
        else
        {
          if ((kk == 2))
          {
            qpk = 0.8872983346207417;
            qwk = 0.2777777777777778;
          }
        }
      }
      z = mesh->zpf[k][0][0];
      InitV = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
      InitV = thermal(0.0,z);
      staticv->hy_dens_cell[k][0][0] = (staticv->hy_dens_cell[k][0][0] + (InitV.hr * qwk));
      staticv->hy_dens_theta_cell[k][0][0] = (staticv->hy_dens_theta_cell[k][0][0] + ((InitV.hr * InitV.ht) * qwk));
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
    z = mesh->zph[k][0][0];
    InitV = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
    InitV = thermal(0.0,z);
    staticv->hy_dens_int[k][0][0] = InitV.hr;
    staticv->hy_dens_theta_int[k][0][0] = (InitV.hr * InitV.ht);
    staticv->hy_pressure_int[k][0][0] = (27.562941092972594 * pow((InitV.hr * InitV.ht),1.400278940027894));
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (DirPara->xd) = 1;
  (DirPara->zd) = 0;
  (DirPara->stepnum) = 0;
}

void SemiDiscreteStep_0(struct HybridStateField* state_init, struct HybridStateField* state_forcing, struct HybridStateField* state_out, struct HybridStaticField* staticv, struct HybridFluxField* flux, struct HybridTendField* tend, struct HybridDirPara* DirPara, double dt)
{
  double hv_coef, r, u, w, t, p;
  struct Vector4 vals, d3_vals;
  int k, j, i;
  if ((DirPara->xd))
  {
    if (((DirPara->stepnum) < 3))
    {
      dt = dt;
      hv_coef = ((-0.05 * 50.0) / (16 * dt));
      vals = (struct Vector4){0.0,0.0,0.0,0.0};
      d3_vals = (struct Vector4){0.0,0.0,0.0,0.0};
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesX_0_20_info ComputeTendenciesX_0_20_para;
      ComputeTendenciesX_0_20_para.lz = 200 - 0;
      ComputeTendenciesX_0_20_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesX_0_20_para.lx = MIN(Proc.lon_end+1, 401) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesX_0_20_para.oz = 0;
      ComputeTendenciesX_0_20_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesX_0_20_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesX_0_20_para.hx = 2;
      ComputeTendenciesX_0_20_para.hy = 0;
      ComputeTendenciesX_0_20_para.hz = 0;
      ComputeTendenciesX_0_20_para.bx = 256;
      ComputeTendenciesX_0_20_para.by = 8;
      ComputeTendenciesX_0_20_para.bz = 1;
      ComputeTendenciesX_0_20_para.mx = 1;
      ComputeTendenciesX_0_20_para.my = 2;
      ComputeTendenciesX_0_20_para.mz = 32;
      ComputeTendenciesX_0_20_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.state_forcingu = &state_forcing->u[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.state_forcingw = &state_forcing->w[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.state_forcingpt = &state_forcing->pt[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.fluxfrho = &flux->frho[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.fluxfu = &flux->fu[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.fluxfw = &flux->fw[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.fluxfpt = &flux->fpt[Proc.lev_hw - 0+ComputeTendenciesX_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_20_para.ox];
      ComputeTendenciesX_0_20_para.hy_dens_cell = &staticv->hy_dens_cell[0][0][0];
      ComputeTendenciesX_0_20_para.hy_dens_theta_cell = &staticv->hy_dens_theta_cell[0][0][0];
      ComputeTendenciesX_0_20_para.hv_coef = hv_coef;
      
      athread_spawn(ComputeTendenciesX_0_20, &ComputeTendenciesX_0_20_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesXfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_frho]);
      HaloWait(Proc,&Proc.FieldReq[async_fu]);
      HaloWait(Proc,&Proc.FieldReq[async_fw]);
      HaloWait(Proc,&Proc.FieldReq[async_fpt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesX_0_21_info ComputeTendenciesX_0_21_para;
      ComputeTendenciesX_0_21_para.lz = 200 - 0;
      ComputeTendenciesX_0_21_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesX_0_21_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesX_0_21_para.oz = 0;
      ComputeTendenciesX_0_21_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesX_0_21_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesX_0_21_para.hx = 1;
      ComputeTendenciesX_0_21_para.hy = 0;
      ComputeTendenciesX_0_21_para.hz = 0;
      ComputeTendenciesX_0_21_para.bx = 256;
      ComputeTendenciesX_0_21_para.by = 8;
      ComputeTendenciesX_0_21_para.bz = 1;
      ComputeTendenciesX_0_21_para.mx = 1;
      ComputeTendenciesX_0_21_para.my = 2;
      ComputeTendenciesX_0_21_para.mz = 32;
      ComputeTendenciesX_0_21_para.fluxfrho = &flux->frho[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.fluxfu = &flux->fu[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.fluxfw = &flux->fw[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.fluxfpt = &flux->fpt[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.tenddrho = &tend->drho[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.tenddu = &tend->du[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.tenddw = &tend->dw[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      ComputeTendenciesX_0_21_para.tenddpt = &tend->dpt[Proc.lev_hw - 0+ComputeTendenciesX_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_21_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_21_para.ox];
      
      athread_spawn(ComputeTendenciesX_0_21, &ComputeTendenciesX_0_21_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesXdrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k + 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k + 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k - 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k - 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k + 2)][j][i] / staticv->hy_dens_cell[(k + 2)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k + 1)][j][i] / staticv->hy_dens_cell[(k + 1)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k - 1)][j][i] / staticv->hy_dens_cell[(k - 1)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k - 2)][j][i] / staticv->hy_dens_cell[(k - 2)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->w[k][j][i] = 0.0;
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZwTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &state_forcing->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->w[k][j][i] = 0.0;
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZwTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &state_forcing->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k + 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k + 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k - 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k - 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      dt = dt;
      hv_coef = ((-0.05 * 50.0) / (16 * dt));
      vals = (struct Vector4){0.0,0.0,0.0,0.0};
      d3_vals = (struct Vector4){0.0,0.0,0.0,0.0};
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            vals.r = ((((-state_forcing->rho[(k - 2)][j][i] / 12) + ((7 * state_forcing->rho[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->rho[k][j][i]) / 12)) - (state_forcing->rho[(k + 1)][j][i] / 12));
            d3_vals.r = (((-state_forcing->rho[(k - 2)][j][i] + (3 * state_forcing->rho[(k - 1)][j][i])) - (3 * state_forcing->rho[k][j][i])) + state_forcing->rho[(k + 1)][j][i]);
            vals.u = ((((-state_forcing->u[(k - 2)][j][i] / 12) + ((7 * state_forcing->u[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->u[k][j][i]) / 12)) - (state_forcing->u[(k + 1)][j][i] / 12));
            d3_vals.u = (((-state_forcing->u[(k - 2)][j][i] + (3 * state_forcing->u[(k - 1)][j][i])) - (3 * state_forcing->u[k][j][i])) + state_forcing->u[(k + 1)][j][i]);
            vals.w = ((((-state_forcing->w[(k - 2)][j][i] / 12) + ((7 * state_forcing->w[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->w[k][j][i]) / 12)) - (state_forcing->w[(k + 1)][j][i] / 12));
            d3_vals.w = (((-state_forcing->w[(k - 2)][j][i] + (3 * state_forcing->w[(k - 1)][j][i])) - (3 * state_forcing->w[k][j][i])) + state_forcing->w[(k + 1)][j][i]);
            vals.t = ((((-state_forcing->pt[(k - 2)][j][i] / 12) + ((7 * state_forcing->pt[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->pt[k][j][i]) / 12)) - (state_forcing->pt[(k + 1)][j][i] / 12));
            d3_vals.t = (((-state_forcing->pt[(k - 2)][j][i] + (3 * state_forcing->pt[(k - 1)][j][i])) - (3 * state_forcing->pt[k][j][i])) + state_forcing->pt[(k + 1)][j][i]);
            r = (vals.r + staticv->hy_dens_int[k][0][0]);
            u = (vals.u / r);
            w = 0.0;
            d3_vals.r = 0.0;
            t = ((vals.t + staticv->hy_dens_theta_int[k][0][0]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - staticv->hy_pressure_int[k][0][0]);
            flux->frho[k][j][i] = ((r * w) - (hv_coef * d3_vals.r));
            flux->fu[k][j][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            flux->fw[k][j][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            flux->fpt[k][j][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesZ_0_20_info ComputeTendenciesZ_0_20_para;
      ComputeTendenciesZ_0_20_para.lz = 201 - 1;
      ComputeTendenciesZ_0_20_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesZ_0_20_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesZ_0_20_para.oz = 1;
      ComputeTendenciesZ_0_20_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesZ_0_20_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesZ_0_20_para.hx = 2;
      ComputeTendenciesZ_0_20_para.hy = 0;
      ComputeTendenciesZ_0_20_para.hz = 2;
      ComputeTendenciesZ_0_20_para.bx = 256;
      ComputeTendenciesZ_0_20_para.by = 1;
      ComputeTendenciesZ_0_20_para.bz = 4;
      ComputeTendenciesZ_0_20_para.mx = 1;
      ComputeTendenciesZ_0_20_para.my = 8;
      ComputeTendenciesZ_0_20_para.mz = 8;
      ComputeTendenciesZ_0_20_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.state_forcingu = &state_forcing->u[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.state_forcingw = &state_forcing->w[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.state_forcingpt = &state_forcing->pt[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.fluxfrho = &flux->frho[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.fluxfu = &flux->fu[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.fluxfw = &flux->fw[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.fluxfpt = &flux->fpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_20_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_20_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_20_para.ox];
      ComputeTendenciesZ_0_20_para.hy_dens_int = &staticv->hy_dens_int[0][0][0];
      ComputeTendenciesZ_0_20_para.hy_dens_theta_int = &staticv->hy_dens_theta_int[0][0][0];
      ComputeTendenciesZ_0_20_para.hy_pressure_int = &staticv->hy_pressure_int[0][0][0];
      ComputeTendenciesZ_0_20_para.hv_coef = hv_coef;
      
      athread_spawn(ComputeTendenciesZ_0_20, &ComputeTendenciesZ_0_20_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            vals.r = ((((-state_forcing->rho[(k - 2)][j][i] / 12) + ((7 * state_forcing->rho[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->rho[k][j][i]) / 12)) - (state_forcing->rho[(k + 1)][j][i] / 12));
            d3_vals.r = (((-state_forcing->rho[(k - 2)][j][i] + (3 * state_forcing->rho[(k - 1)][j][i])) - (3 * state_forcing->rho[k][j][i])) + state_forcing->rho[(k + 1)][j][i]);
            vals.u = ((((-state_forcing->u[(k - 2)][j][i] / 12) + ((7 * state_forcing->u[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->u[k][j][i]) / 12)) - (state_forcing->u[(k + 1)][j][i] / 12));
            d3_vals.u = (((-state_forcing->u[(k - 2)][j][i] + (3 * state_forcing->u[(k - 1)][j][i])) - (3 * state_forcing->u[k][j][i])) + state_forcing->u[(k + 1)][j][i]);
            vals.w = ((((-state_forcing->w[(k - 2)][j][i] / 12) + ((7 * state_forcing->w[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->w[k][j][i]) / 12)) - (state_forcing->w[(k + 1)][j][i] / 12));
            d3_vals.w = (((-state_forcing->w[(k - 2)][j][i] + (3 * state_forcing->w[(k - 1)][j][i])) - (3 * state_forcing->w[k][j][i])) + state_forcing->w[(k + 1)][j][i]);
            vals.t = ((((-state_forcing->pt[(k - 2)][j][i] / 12) + ((7 * state_forcing->pt[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->pt[k][j][i]) / 12)) - (state_forcing->pt[(k + 1)][j][i] / 12));
            d3_vals.t = (((-state_forcing->pt[(k - 2)][j][i] + (3 * state_forcing->pt[(k - 1)][j][i])) - (3 * state_forcing->pt[k][j][i])) + state_forcing->pt[(k + 1)][j][i]);
            r = (vals.r + staticv->hy_dens_int[k][0][0]);
            u = (vals.u / r);
            w = 0.0;
            d3_vals.r = 0.0;
            t = ((vals.t + staticv->hy_dens_theta_int[k][0][0]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - staticv->hy_pressure_int[k][0][0]);
            flux->frho[k][j][i] = ((r * w) - (hv_coef * d3_vals.r));
            flux->fu[k][j][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            flux->fw[k][j][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            flux->fpt[k][j][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &flux->frho[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_frho], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fu[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fu], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fw[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fw], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fpt[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fpt], false, true, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_frho]);
      HaloWait(Proc,&Proc.FieldReq[async_fu]);
      HaloWait(Proc,&Proc.FieldReq[async_fw]);
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_fpt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesZ_0_21_info ComputeTendenciesZ_0_21_para;
      ComputeTendenciesZ_0_21_para.lz = 200 - 0;
      ComputeTendenciesZ_0_21_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesZ_0_21_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesZ_0_21_para.oz = 0;
      ComputeTendenciesZ_0_21_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesZ_0_21_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesZ_0_21_para.hx = 2;
      ComputeTendenciesZ_0_21_para.hy = 0;
      ComputeTendenciesZ_0_21_para.hz = 2;
      ComputeTendenciesZ_0_21_para.bx = 256;
      ComputeTendenciesZ_0_21_para.by = 1;
      ComputeTendenciesZ_0_21_para.bz = 4;
      ComputeTendenciesZ_0_21_para.mx = 1;
      ComputeTendenciesZ_0_21_para.my = 8;
      ComputeTendenciesZ_0_21_para.mz = 8;
      ComputeTendenciesZ_0_21_para.fluxfrho = &flux->frho[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.fluxfu = &flux->fu[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.fluxfw = &flux->fw[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.tenddw = &tend->dw[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.fluxfpt = &flux->fpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.tenddrho = &tend->drho[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.tenddu = &tend->du[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      ComputeTendenciesZ_0_21_para.tenddpt = &tend->dpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_21_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_21_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_21_para.ox];
      
      athread_spawn(ComputeTendenciesZ_0_21, &ComputeTendenciesZ_0_21_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZdrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  else
  {
    if (((DirPara->stepnum) < 3))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k + 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k + 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k - 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->rho[k][j][i] = state_forcing->rho[(k - 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->rho[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_rho], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k + 2)][j][i] / staticv->hy_dens_cell[(k + 2)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k + 1)][j][i] / staticv->hy_dens_cell[(k + 1)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k - 1)][j][i] / staticv->hy_dens_cell[(k - 1)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->u[k][j][i] = ((state_forcing->u[(k - 2)][j][i] / staticv->hy_dens_cell[(k - 2)][0][0]) * staticv->hy_dens_cell[k][0][0]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->u[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_u], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->w[k][j][i] = 0.0;
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZwTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &state_forcing->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->w[k][j][i] = 0.0;
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZwTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &state_forcing->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k + 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[-2 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k + 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[-1 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k - 1)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=201+Proc.lev_hw; k<202+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            state_forcing->pt[k][j][i] = state_forcing->pt[(k - 2)][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      SetBoundaryZptTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &state_forcing->pt[201 + Proc.lev_hw][0][0], &Proc.FieldReq[async_pt], true, true, true, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      dt = dt;
      hv_coef = ((-0.05 * 50.0) / (16 * dt));
      vals = (struct Vector4){0.0,0.0,0.0,0.0};
      d3_vals = (struct Vector4){0.0,0.0,0.0,0.0};
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            vals.r = ((((-state_forcing->rho[(k - 2)][j][i] / 12) + ((7 * state_forcing->rho[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->rho[k][j][i]) / 12)) - (state_forcing->rho[(k + 1)][j][i] / 12));
            d3_vals.r = (((-state_forcing->rho[(k - 2)][j][i] + (3 * state_forcing->rho[(k - 1)][j][i])) - (3 * state_forcing->rho[k][j][i])) + state_forcing->rho[(k + 1)][j][i]);
            vals.u = ((((-state_forcing->u[(k - 2)][j][i] / 12) + ((7 * state_forcing->u[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->u[k][j][i]) / 12)) - (state_forcing->u[(k + 1)][j][i] / 12));
            d3_vals.u = (((-state_forcing->u[(k - 2)][j][i] + (3 * state_forcing->u[(k - 1)][j][i])) - (3 * state_forcing->u[k][j][i])) + state_forcing->u[(k + 1)][j][i]);
            vals.w = ((((-state_forcing->w[(k - 2)][j][i] / 12) + ((7 * state_forcing->w[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->w[k][j][i]) / 12)) - (state_forcing->w[(k + 1)][j][i] / 12));
            d3_vals.w = (((-state_forcing->w[(k - 2)][j][i] + (3 * state_forcing->w[(k - 1)][j][i])) - (3 * state_forcing->w[k][j][i])) + state_forcing->w[(k + 1)][j][i]);
            vals.t = ((((-state_forcing->pt[(k - 2)][j][i] / 12) + ((7 * state_forcing->pt[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->pt[k][j][i]) / 12)) - (state_forcing->pt[(k + 1)][j][i] / 12));
            d3_vals.t = (((-state_forcing->pt[(k - 2)][j][i] + (3 * state_forcing->pt[(k - 1)][j][i])) - (3 * state_forcing->pt[k][j][i])) + state_forcing->pt[(k + 1)][j][i]);
            r = (vals.r + staticv->hy_dens_int[k][0][0]);
            u = (vals.u / r);
            w = 0.0;
            d3_vals.r = 0.0;
            t = ((vals.t + staticv->hy_dens_theta_int[k][0][0]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - staticv->hy_pressure_int[k][0][0]);
            flux->frho[k][j][i] = ((r * w) - (hv_coef * d3_vals.r));
            flux->fu[k][j][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            flux->fw[k][j][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            flux->fpt[k][j][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesZ_0_22_info ComputeTendenciesZ_0_22_para;
      ComputeTendenciesZ_0_22_para.lz = 201 - 1;
      ComputeTendenciesZ_0_22_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesZ_0_22_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesZ_0_22_para.oz = 1;
      ComputeTendenciesZ_0_22_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesZ_0_22_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesZ_0_22_para.hx = 2;
      ComputeTendenciesZ_0_22_para.hy = 0;
      ComputeTendenciesZ_0_22_para.hz = 2;
      ComputeTendenciesZ_0_22_para.bx = 256;
      ComputeTendenciesZ_0_22_para.by = 1;
      ComputeTendenciesZ_0_22_para.bz = 4;
      ComputeTendenciesZ_0_22_para.mx = 1;
      ComputeTendenciesZ_0_22_para.my = 8;
      ComputeTendenciesZ_0_22_para.mz = 8;
      ComputeTendenciesZ_0_22_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.state_forcingu = &state_forcing->u[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.state_forcingw = &state_forcing->w[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.state_forcingpt = &state_forcing->pt[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.fluxfrho = &flux->frho[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.fluxfu = &flux->fu[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.fluxfw = &flux->fw[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.fluxfpt = &flux->fpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_22_para.ox];
      ComputeTendenciesZ_0_22_para.hy_dens_int = &staticv->hy_dens_int[0][0][0];
      ComputeTendenciesZ_0_22_para.hy_dens_theta_int = &staticv->hy_dens_theta_int[0][0][0];
      ComputeTendenciesZ_0_22_para.hy_pressure_int = &staticv->hy_pressure_int[0][0][0];
      ComputeTendenciesZ_0_22_para.hv_coef = hv_coef;
      
      athread_spawn(ComputeTendenciesZ_0_22, &ComputeTendenciesZ_0_22_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=200+Proc.lev_hw; k<201+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            vals.r = ((((-state_forcing->rho[(k - 2)][j][i] / 12) + ((7 * state_forcing->rho[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->rho[k][j][i]) / 12)) - (state_forcing->rho[(k + 1)][j][i] / 12));
            d3_vals.r = (((-state_forcing->rho[(k - 2)][j][i] + (3 * state_forcing->rho[(k - 1)][j][i])) - (3 * state_forcing->rho[k][j][i])) + state_forcing->rho[(k + 1)][j][i]);
            vals.u = ((((-state_forcing->u[(k - 2)][j][i] / 12) + ((7 * state_forcing->u[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->u[k][j][i]) / 12)) - (state_forcing->u[(k + 1)][j][i] / 12));
            d3_vals.u = (((-state_forcing->u[(k - 2)][j][i] + (3 * state_forcing->u[(k - 1)][j][i])) - (3 * state_forcing->u[k][j][i])) + state_forcing->u[(k + 1)][j][i]);
            vals.w = ((((-state_forcing->w[(k - 2)][j][i] / 12) + ((7 * state_forcing->w[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->w[k][j][i]) / 12)) - (state_forcing->w[(k + 1)][j][i] / 12));
            d3_vals.w = (((-state_forcing->w[(k - 2)][j][i] + (3 * state_forcing->w[(k - 1)][j][i])) - (3 * state_forcing->w[k][j][i])) + state_forcing->w[(k + 1)][j][i]);
            vals.t = ((((-state_forcing->pt[(k - 2)][j][i] / 12) + ((7 * state_forcing->pt[(k - 1)][j][i]) / 12)) + ((7 * state_forcing->pt[k][j][i]) / 12)) - (state_forcing->pt[(k + 1)][j][i] / 12));
            d3_vals.t = (((-state_forcing->pt[(k - 2)][j][i] + (3 * state_forcing->pt[(k - 1)][j][i])) - (3 * state_forcing->pt[k][j][i])) + state_forcing->pt[(k + 1)][j][i]);
            r = (vals.r + staticv->hy_dens_int[k][0][0]);
            u = (vals.u / r);
            w = 0.0;
            d3_vals.r = 0.0;
            t = ((vals.t + staticv->hy_dens_theta_int[k][0][0]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - staticv->hy_pressure_int[k][0][0]);
            flux->frho[k][j][i] = ((r * w) - (hv_coef * d3_vals.r));
            flux->fu[k][j][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            flux->fw[k][j][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            flux->fpt[k][j][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_2d_D(Proc, &flux->frho[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_frho], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fu[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fu], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fw[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fw], false, true, false, true, false, false, true);
      UpdateHalo_2d_D(Proc, &flux->fpt[200 + Proc.lev_hw][0][0], &Proc.FieldReq[async_fpt], false, true, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_frho]);
      HaloWait(Proc,&Proc.FieldReq[async_fu]);
      HaloWait(Proc,&Proc.FieldReq[async_fw]);
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_fpt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesZ_0_23_info ComputeTendenciesZ_0_23_para;
      ComputeTendenciesZ_0_23_para.lz = 200 - 0;
      ComputeTendenciesZ_0_23_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesZ_0_23_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesZ_0_23_para.oz = 0;
      ComputeTendenciesZ_0_23_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesZ_0_23_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesZ_0_23_para.hx = 2;
      ComputeTendenciesZ_0_23_para.hy = 0;
      ComputeTendenciesZ_0_23_para.hz = 2;
      ComputeTendenciesZ_0_23_para.bx = 256;
      ComputeTendenciesZ_0_23_para.by = 1;
      ComputeTendenciesZ_0_23_para.bz = 4;
      ComputeTendenciesZ_0_23_para.mx = 1;
      ComputeTendenciesZ_0_23_para.my = 8;
      ComputeTendenciesZ_0_23_para.mz = 8;
      ComputeTendenciesZ_0_23_para.fluxfrho = &flux->frho[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.fluxfu = &flux->fu[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.fluxfw = &flux->fw[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.tenddw = &tend->dw[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.fluxfpt = &flux->fpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.tenddrho = &tend->drho[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.tenddu = &tend->du[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      ComputeTendenciesZ_0_23_para.tenddpt = &tend->dpt[Proc.lev_hw - 2+ComputeTendenciesZ_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesZ_0_23_para.oy][Proc.lon_hw -2+ComputeTendenciesZ_0_23_para.ox];
      
      athread_spawn(ComputeTendenciesZ_0_23, &ComputeTendenciesZ_0_23_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesZdrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      dt = dt;
      hv_coef = ((-0.05 * 50.0) / (16 * dt));
      vals = (struct Vector4){0.0,0.0,0.0,0.0};
      d3_vals = (struct Vector4){0.0,0.0,0.0,0.0};
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_rho]);
      HaloWait(Proc,&Proc.FieldReq[async_u]);
      HaloWait(Proc,&Proc.FieldReq[async_w]);
      HaloWait(Proc,&Proc.FieldReq[async_pt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesX_0_22_info ComputeTendenciesX_0_22_para;
      ComputeTendenciesX_0_22_para.lz = 200 - 0;
      ComputeTendenciesX_0_22_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesX_0_22_para.lx = MIN(Proc.lon_end+1, 401) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesX_0_22_para.oz = 0;
      ComputeTendenciesX_0_22_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesX_0_22_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesX_0_22_para.hx = 2;
      ComputeTendenciesX_0_22_para.hy = 0;
      ComputeTendenciesX_0_22_para.hz = 0;
      ComputeTendenciesX_0_22_para.bx = 256;
      ComputeTendenciesX_0_22_para.by = 8;
      ComputeTendenciesX_0_22_para.bz = 1;
      ComputeTendenciesX_0_22_para.mx = 1;
      ComputeTendenciesX_0_22_para.my = 2;
      ComputeTendenciesX_0_22_para.mz = 32;
      ComputeTendenciesX_0_22_para.state_forcingrho = &state_forcing->rho[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.state_forcingu = &state_forcing->u[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.state_forcingw = &state_forcing->w[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.state_forcingpt = &state_forcing->pt[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.fluxfrho = &flux->frho[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.fluxfu = &flux->fu[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.fluxfw = &flux->fw[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.fluxfpt = &flux->fpt[Proc.lev_hw - 0+ComputeTendenciesX_0_22_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_22_para.oy][Proc.lon_hw -2+ComputeTendenciesX_0_22_para.ox];
      ComputeTendenciesX_0_22_para.hy_dens_cell = &staticv->hy_dens_cell[0][0][0];
      ComputeTendenciesX_0_22_para.hy_dens_theta_cell = &staticv->hy_dens_theta_cell[0][0][0];
      ComputeTendenciesX_0_22_para.hv_coef = hv_coef;
      
      athread_spawn(ComputeTendenciesX_0_22, &ComputeTendenciesX_0_22_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesXfrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &flux->frho[0][0][0], &Proc.FieldReq[async_frho], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fu[0][0][0], &Proc.FieldReq[async_fu], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fw[0][0][0], &Proc.FieldReq[async_fw], false, true, false, false, true, false, false, true);
      UpdateHalo_3d_D(Proc, &flux->fpt[0][0][0], &Proc.FieldReq[async_fpt], false, true, false, false, true, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_frho]);
      HaloWait(Proc,&Proc.FieldReq[async_fu]);
      HaloWait(Proc,&Proc.FieldReq[async_fw]);
      HaloWait(Proc,&Proc.FieldReq[async_fpt]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      ComputeTendenciesX_0_23_info ComputeTendenciesX_0_23_para;
      ComputeTendenciesX_0_23_para.lz = 200 - 0;
      ComputeTendenciesX_0_23_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
      ComputeTendenciesX_0_23_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
      ComputeTendenciesX_0_23_para.oz = 0;
      ComputeTendenciesX_0_23_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
      ComputeTendenciesX_0_23_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
      ComputeTendenciesX_0_23_para.hx = 1;
      ComputeTendenciesX_0_23_para.hy = 0;
      ComputeTendenciesX_0_23_para.hz = 0;
      ComputeTendenciesX_0_23_para.bx = 256;
      ComputeTendenciesX_0_23_para.by = 8;
      ComputeTendenciesX_0_23_para.bz = 1;
      ComputeTendenciesX_0_23_para.mx = 1;
      ComputeTendenciesX_0_23_para.my = 2;
      ComputeTendenciesX_0_23_para.mz = 32;
      ComputeTendenciesX_0_23_para.fluxfrho = &flux->frho[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.fluxfu = &flux->fu[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.fluxfw = &flux->fw[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.fluxfpt = &flux->fpt[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.tenddrho = &tend->drho[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.tenddu = &tend->du[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.tenddw = &tend->dw[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      ComputeTendenciesX_0_23_para.tenddpt = &tend->dpt[Proc.lev_hw - 0+ComputeTendenciesX_0_23_para.oz][Proc.lat_hw -0+ComputeTendenciesX_0_23_para.oy][Proc.lon_hw -1+ComputeTendenciesX_0_23_para.ox];
      
      athread_spawn(ComputeTendenciesX_0_23, &ComputeTendenciesX_0_23_para);
      athread_join();
      
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      ComputeTendenciesXdrhoTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_rho]);
  HaloWait(Proc,&Proc.FieldReq[async_u]);
  HaloWait(Proc,&Proc.FieldReq[async_w]);
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateState_0_5_info UpdateState_0_5_para;
  UpdateState_0_5_para.lz = 200 - 0;
  UpdateState_0_5_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
  UpdateState_0_5_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
  UpdateState_0_5_para.oz = 0;
  UpdateState_0_5_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  UpdateState_0_5_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  UpdateState_0_5_para.hx = 2;
  UpdateState_0_5_para.hy = 0;
  UpdateState_0_5_para.hz = 0;
  UpdateState_0_5_para.bx = 256;
  UpdateState_0_5_para.by = 8;
  UpdateState_0_5_para.bz = 1;
  UpdateState_0_5_para.mx = 1;
  UpdateState_0_5_para.my = 2;
  UpdateState_0_5_para.mz = 32;
  UpdateState_0_5_para.state_initrho = &state_init->rho[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.tenddrho = &tend->drho[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_initu = &state_init->u[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.tenddu = &tend->du[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_initw = &state_init->w[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.tenddw = &tend->dw[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_initpt = &state_init->pt[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.tenddpt = &tend->dpt[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_outrho = &state_out->rho[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_outu = &state_out->u[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_outw = &state_out->w[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.state_outpt = &state_out->pt[Proc.lev_hw - 0+UpdateState_0_5_para.oz][Proc.lat_hw -0+UpdateState_0_5_para.oy][Proc.lon_hw -2+UpdateState_0_5_para.ox];
  UpdateState_0_5_para.dt = dt;
  
  athread_spawn(UpdateState_0_5, &UpdateState_0_5_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  UpdateStaterhoTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_out->rho[0][0][0], &Proc.FieldReq[async_rho], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_out->u[0][0][0], &Proc.FieldReq[async_u], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_out->w[0][0][0], &Proc.FieldReq[async_w], true, true, true, true, true, false, false, true);
  UpdateHalo_3d_D(Proc, &state_out->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (DirPara->stepnum) = ((DirPara->stepnum) + 1);
}

void SwitchDirection_0(struct HybridDirPara* DirPara)
{
  if ((DirPara->xd))
  {
    (DirPara->xd) = 0;
    (DirPara->zd) = 1;
  }
  else
  {
    (DirPara->xd) = 1;
    (DirPara->zd) = 0;
  }
  (DirPara->stepnum) = 0;
}

void OutputPrepare_0(struct HybridStateField* state, struct HybridOutField* stateout, struct HybridStaticField* staticv, struct HybridOutPara* OutPara)
{
  int k, j, i, nid;
  (OutPara->step) = ((OutPara->step) + 1);
  if (((OutPara->step) == (OutPara->interval)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_rho]);
    HaloWait(Proc,&Proc.FieldReq[async_u]);
    HaloWait(Proc,&Proc.FieldReq[async_w]);
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    OutputPrepare_0_0_info OutputPrepare_0_0_para;
    OutputPrepare_0_0_para.lz = 200 - 0;
    OutputPrepare_0_0_para.ly = MIN(Proc.lat_end+1, 1) - MAX(Proc.lat_beg, 0);
    OutputPrepare_0_0_para.lx = MIN(Proc.lon_end+1, 400) - MAX(Proc.lon_beg, 0);
    OutputPrepare_0_0_para.oz = 0;
    OutputPrepare_0_0_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
    OutputPrepare_0_0_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
    OutputPrepare_0_0_para.hx = 2;
    OutputPrepare_0_0_para.hy = 0;
    OutputPrepare_0_0_para.hz = 0;
    OutputPrepare_0_0_para.bx = 256;
    OutputPrepare_0_0_para.by = 8;
    OutputPrepare_0_0_para.bz = 1;
    OutputPrepare_0_0_para.mx = 1;
    OutputPrepare_0_0_para.my = 2;
    OutputPrepare_0_0_para.mz = 32;
    OutputPrepare_0_0_para.staterho = &state->rho[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.stateu = &state->u[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.statew = &state->w[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.statept = &state->pt[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.stateoutdens = &stateout->dens[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.stateoutuwnd = &stateout->uwnd[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.stateoutwwnd = &stateout->wwnd[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.stateouttheta = &stateout->theta[Proc.lev_hw - 0+OutputPrepare_0_0_para.oz][Proc.lat_hw -0+OutputPrepare_0_0_para.oy][Proc.lon_hw -2+OutputPrepare_0_0_para.ox];
    OutputPrepare_0_0_para.hy_dens_cell = &staticv->hy_dens_cell[0][0][0];
    OutputPrepare_0_0_para.hy_dens_theta_cell = &staticv->hy_dens_theta_cell[0][0][0];
    
    athread_spawn(OutputPrepare_0_0, &OutputPrepare_0_0_para);
    athread_join();
    
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    OutputPreparedensTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    nid = ncOpenFile("output.nc",1);
    ncDefDim(nid,200,1,400);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    ncDefDimTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    ncDefVar(nid, "dens", 1, 1, 1);
    ncPutVar(nid, "dens", 0, 200, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(stateout->dens[0][0][0]));
    
    ncDefVar(nid, "uwnd", 1, 1, 1);
    ncPutVar(nid, "uwnd", 0, 200, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(stateout->uwnd[0][0][0]));
    
    ncDefVar(nid, "wwnd", 1, 1, 1);
    ncPutVar(nid, "wwnd", 0, 200, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(stateout->wwnd[0][0][0]));
    
    ncDefVar(nid, "theta", 1, 1, 1);
    ncPutVar(nid, "theta", 0, 200, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(stateout->theta[0][0][0]));
    
    ncCloseFile(nid);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    ncCloseFileTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (OutPara->step) = 0;
  }
}

void MeshInit_0(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh)
{
  MeshInit(&global_mesh[0]);
  MeshInit_cp(&global_mesh[0], &mesh[0]);
}

void timeOperatorInit_0(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridMeshField* mesh, struct HybridDirPara* DirPara)
{
  spaceOperatorInit_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &mesh[0], DirPara);
}

void performTimestep_0(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridFluxField* flux, struct HybridTendField* tend, struct HybridMeshField* mesh, struct HybridDirPara* DirPara, struct HybridOutField* stateoutput, struct HybridOutPara* OutPara, double dt, double dtd2, double dtd3)
{
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &flux[0], &tend[0], DirPara, dtd3);
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[Timeinfo.newt], &staticv[0], &flux[0], &tend[0], DirPara, dtd2);
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[Timeinfo.oldt], &staticv[0], &flux[0], &tend[0], DirPara, dt);
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &flux[0], &tend[0], DirPara, dtd3);
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[Timeinfo.newt], &staticv[0], &flux[0], &tend[0], DirPara, dtd2);
  SemiDiscreteStep_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[Timeinfo.oldt], &staticv[0], &flux[0], &tend[0], DirPara, dt);
  SwitchDirection_0(DirPara);
  OutputPrepare_0(&state[Timeinfo.oldt], &stateoutput[Timeinfo.oldt], &staticv[0], OutPara);
}

int main(int argc, char **argv)
{
  int size,rank;
  double tottime_beg, tottime_end;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ProcInit_LonLat_Domain(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev);
  
  //Proc Ngb Init
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 0){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,2, 0, 2, 2, 2, DTYPE_DOUBLE, 26);
  }
  
  PhysicalVariableInit();
  
  TimeInit();
  
  athread_init();
  
  global_info gi;
  gi.fnx = Proc.full_nlon;
  gi.fny = Proc.full_nlat;
  gi.fnz = Proc.full_nlev;
  gi.hnx = Proc.half_nlon;
  gi.hny = Proc.half_nlat;
  gi.hnz = Proc.half_nlev;
  gi.ghx = Proc.lon_hw;
  gi.ghy = Proc.lat_hw;
  gi.ghz = Proc.lev_hw;
  athread_spawn(global_prepare, &gi);
  athread_join();
  
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 0)
  {
    MeshInit_0(global_mesh, mesh);
    timeOperatorInit_0(state, staticv, mesh, &DirPara);
    tottime_beg = MPI_Wtime();
    for (int t = 0; t < 12; t += 1){
      performTimestep_0(state, staticv, flux, tend, mesh, &DirPara, stateout, &OutPara, 0.16666666666666666, 0.08333333333333333, 0.05555555555555555);
    }
    tottime_end = MPI_Wtime();
  }
  
  ProfilingOutPut(Proc.id);
  if (Proc.id == 0) printf("Time is %.2f\n",tottime_end-tottime_beg);
}
