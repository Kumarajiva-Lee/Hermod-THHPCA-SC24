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
extern SLAVE_FUN(updateState_0_0)(void*);
extern SLAVE_FUN(updateState_0_1)(void*);
extern SLAVE_FUN(updateState_0_2)(void*);
extern SLAVE_FUN(updateState_0_3)(void*);
extern SLAVE_FUN(updateState_0_4)(void*);
extern SLAVE_FUN(updateState_0_5)(void*);
extern SLAVE_FUN(updateState_0_6)(void*);
extern SLAVE_FUN(updateState_0_7)(void*);
extern SLAVE_FUN(updateState_0_8)(void*);
extern SLAVE_FUN(update_rk1_0_0)(void*);
extern SLAVE_FUN(update_rk1_0_1)(void*);
extern SLAVE_FUN(update_rk1_0_2)(void*);
extern SLAVE_FUN(update_rk1_0_3)(void*);
extern SLAVE_FUN(update_rk1_0_4)(void*);
extern SLAVE_FUN(update_rk1_0_5)(void*);
extern SLAVE_FUN(update_rk1_0_6)(void*);
extern SLAVE_FUN(update_rk1_0_7)(void*);
extern SLAVE_FUN(update_rk1_0_8)(void*);
extern SLAVE_FUN(update_rk2_0_0)(void*);
extern SLAVE_FUN(update_rk2_0_1)(void*);
extern SLAVE_FUN(update_rk2_0_2)(void*);
extern SLAVE_FUN(update_rk2_0_3)(void*);
extern SLAVE_FUN(update_rk2_0_4)(void*);
extern SLAVE_FUN(update_rk2_0_5)(void*);
extern SLAVE_FUN(update_rk2_0_6)(void*);
extern SLAVE_FUN(update_rk2_0_7)(void*);
extern SLAVE_FUN(update_rk2_0_8)(void*);
extern SLAVE_FUN(updateX_0_30)(void*);
extern SLAVE_FUN(updateX_0_31)(void*);
extern SLAVE_FUN(updateX_0_32)(void*);
extern SLAVE_FUN(updateX_0_33)(void*);
extern SLAVE_FUN(updateX_0_34)(void*);
extern SLAVE_FUN(updateX_0_35)(void*);
extern SLAVE_FUN(updateX_0_36)(void*);
extern SLAVE_FUN(updateX_0_37)(void*);
extern SLAVE_FUN(updateX_0_38)(void*);
extern SLAVE_FUN(updateX_0_39)(void*);
extern SLAVE_FUN(updateX_0_40)(void*);
extern SLAVE_FUN(updateX_0_41)(void*);
extern SLAVE_FUN(updateX_0_42)(void*);
extern SLAVE_FUN(updateX_0_43)(void*);
extern SLAVE_FUN(updateX_0_44)(void*);
extern SLAVE_FUN(updateY_0_30)(void*);
extern SLAVE_FUN(updateY_0_31)(void*);
extern SLAVE_FUN(updateY_0_32)(void*);
extern SLAVE_FUN(updateY_0_33)(void*);
extern SLAVE_FUN(updateY_0_34)(void*);
extern SLAVE_FUN(updateY_0_35)(void*);
extern SLAVE_FUN(updateY_0_36)(void*);
extern SLAVE_FUN(updateY_0_37)(void*);
extern SLAVE_FUN(updateY_0_38)(void*);
extern SLAVE_FUN(updateY_0_39)(void*);
extern SLAVE_FUN(updateY_0_40)(void*);
extern SLAVE_FUN(updateY_0_41)(void*);
extern SLAVE_FUN(updateY_0_42)(void*);
extern SLAVE_FUN(updateY_0_43)(void*);
extern SLAVE_FUN(updateY_0_44)(void*);
extern SLAVE_FUN(update_rk3_0_0)(void*);
extern SLAVE_FUN(update_rk3_0_1)(void*);
extern SLAVE_FUN(update_rk3_0_2)(void*);
extern SLAVE_FUN(update_rk3_0_3)(void*);
extern SLAVE_FUN(update_rk3_0_4)(void*);
extern SLAVE_FUN(update_rk3_0_5)(void*);
extern SLAVE_FUN(update_rk3_0_6)(void*);
extern SLAVE_FUN(update_rk3_0_7)(void*);
extern SLAVE_FUN(update_rk3_0_8)(void*);

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

struct Vector4 matrixA(double k, double lamb, double theta)
{
  struct Vector4 ret;
  double alambda, atheta, a, b, c, d, temp;
  ret = (struct Vector4){0.0,0.0,0.0,0.0};
  if ((k <= 4))
  {
    alambda = (lamb - (((k - 1) * 3.141592653589793) / 2.0));
    atheta = theta;
    a = sin(alambda);
    b = cos(alambda);
    c = sin(atheta);
    d = cos(atheta);
    ret.a11 = d;
    ret.a12 = 0.0;
    ret.a21 = (((-c * d) * a) / b);
    ret.a22 = (((b * d) * d) + ((c * c) / b));
  }
  else
  {
    if ((k == 5))
    {
      alambda = lamb;
      atheta = theta;
      a = sin(alambda);
      b = cos(alambda);
      c = sin(atheta);
      d = cos(atheta);
      temp = (1.0 + (((((a * a) * d) * d) / c) / c));
      ret.a11 = ((b * c) * temp);
      ret.a21 = (((-c * c) * a) * temp);
      temp = (1.0 + (((((b * b) * d) * d) / c) / c));
      ret.a12 = ((a * c) * temp);
      ret.a22 = (((b * c) * c) * temp);
    }
    else
    {
      alambda = lamb;
      atheta = theta;
      a = sin(alambda);
      b = cos(alambda);
      c = sin(atheta);
      d = cos(atheta);
      temp = (1.0 + (((((a * a) * d) * d) / c) / c));
      ret.a11 = ((-b * c) * temp);
      ret.a21 = (((c * c) * a) * temp);
      temp = (1.0 + (((((b * b) * d) * d) / c) / c));
      ret.a12 = ((a * c) * temp);
      ret.a22 = (((b * c) * c) * temp);
    }
  }
  return ret;
}

struct Vector4 matrixIG(double x, double y)
{
  struct Vector4 ret;
  double rho;
  ret = (struct Vector4){0.0,0.0,0.0,0.0};
  rho = sqrt(((1.0 + pow(tan((x / 6371220.0)),2)) + pow(tan((y / 6371220.0)),2)));
  ret.a11 = ((1.0 + pow(tan((y / 6371220.0)),2)) * ((pow(rho,2) * pow(cos((x / 6371220.0)),2)) * pow(cos((y / 6371220.0)),2)));
  ret.a12 = ((tan((x / 6371220.0)) * tan((y / 6371220.0))) * ((pow(rho,2) * pow(cos((x / 6371220.0)),2)) * pow(cos((y / 6371220.0)),2)));
  ret.a21 = ret.a12;
  ret.a22 = ((1.0 + pow(tan((x / 6371220.0)),2)) * ((pow(rho,2) * pow(cos((x / 6371220.0)),2)) * pow(cos((y / 6371220.0)),2)));
  return ret;
}

struct Vector2 pprop2sp(int k, double x, double y)
{
  double x1, y1, a, b;
  struct Vector2 ret;
  x1 = (x / 6371220.0);
  y1 = (y / 6371220.0);
  ret = (struct Vector2){0.0,0.0};
  if ((k <= 4))
  {
    ret.x = (x1 + (((k - 1) * 3.141592653589793) / 2.0));
    ret.y = atan2((tan(y1) * cos(x1)),1.0);
  }
  else
  {
    if ((k == 5))
    {
      a = tan(x1);
      b = tan(y1);
      ret.x = atan2(a,-b);
      ret.y = atan2(1.0,sqrt(((a * a) + (b * b))));
    }
    else
    {
      a = tan(x1);
      b = tan(y1);
      ret.x = atan2(a,b);
      ret.y = -atan2(1.0,sqrt(((a * a) + (b * b))));
    }
  }
  return ret;
}

struct Vector2 covprosp2p(int k, double sv1, double sv2, double lamb, double theta)
{
  struct Vector4 a;
  struct Vector2 ret;
  a = (struct Vector4){0.0,0.0,0.0,0.0};
  a = matrixA(k,lamb,theta);
  ret = (struct Vector2){0.0,0.0};
  ret.x = ((a.a11 * sv1) + (a.a21 * sv2));
  ret.y = ((a.a12 * sv1) + (a.a22 * sv2));
  return ret;
}

struct Vector2 cov2contrav(double cov1, double cov2, double x, double y)
{
  struct Vector4 a;
  struct Vector2 ret;
  a = (struct Vector4){0.0,0.0,0.0,0.0};
  a = matrixIG(x,y);
  ret = (struct Vector2){0.0,0.0};
  ret.x = ((a.a11 * cov1) + (a.a12 * cov2));
  ret.y = ((a.a21 * cov1) + (a.a22 * cov2));
  return ret;
}

double computeh(double lamb, double theta)
{
  double h;
  h = (29400.0 - ((((6371220.0 * 7.292e-05) * 38.61068276698372) + ((0.5 * 38.61068276698372) * 38.61068276698372)) * pow((((-cos(lamb) * cos(theta)) * sin(0)) + (sin(theta) * cos(0))),2)));
  h = (h / 9.80616);
  return h;
}

double computejab(double x, double y)
{
  double alpha, beta, jab;
  alpha = (x / 6371220.0);
  beta = (y / 6371220.0);
  jab = (((((1.0 / pow(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))),1.5)) / cos(alpha)) / cos(alpha)) / cos(beta)) / cos(beta));
  return jab;
}

struct Vector2 computevel(int k, double lamb, double theta)
{
  double alpha, us, vs;
  struct Vector2 ret;
  alpha = 0.0;
  us = (38.61068276698372 * ((cos(alpha) * cos(theta)) + ((sin(alpha) * cos(lamb)) * sin(theta))));
  vs = ((-38.61068276698372 * sin(alpha)) * sin(lamb));
  ret = covprosp2p(k,us,vs,lamb,theta);
  return ret;
}

void MeshInit(struct HybridMeshField* mesh)
{
  int i, j;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=0; i<600+2*Proc.lon_hw; i+=1){
    mesh->x_1[0][0][i] = (((-6371220.0 * 3.141592653589793) / 4.0) + ((i - 1) * 16679.814955336966));
    mesh->x_2[0][0][i] = (mesh->x_1[0][0][i] + 8339.907477668483);
    mesh->x_3[0][0][i] = (mesh->x_2[0][0][i] + 8339.907477668483);
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
  
  for(j=0; j<600+2*Proc.lat_hw; j+=1){
    mesh->y_1[0][j][0] = (((-6371220.0 * 3.141592653589793) / 4.0) + ((j - 1) * 16679.814955336966));
    mesh->y_2[0][j][0] = (mesh->y_1[0][j][0] + 8339.907477668483);
    mesh->y_3[0][j][0] = (mesh->y_2[0][j][0] + 8339.907477668483);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  MeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void spaceOperatorInit_0(struct HybridStateField* state, struct HybridMeshField* mesh)
{
  int pypanel, j, i;
  struct Vector2 ltret, uv, uvc;
  double lamb, theta, h, jab, alpha, beta, rho;
  pypanel = Proc.panel;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_1[0][0][i],mesh->y_1[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_1[0][0][i],mesh->y_1[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_11[0][j][i] = (jab * h);
      state->u_11[0][j][i] = uv.x;
      state->v_11[0][j][i] = uv.y;
      state->jab_11[0][j][i] = jab;
      alpha = (mesh->x_1[0][0][i] / 6371220.0);
      beta = (mesh->y_1[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_11[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_11[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_11[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_11[0][j][i],state->v_11[0][j][i],mesh->x_1[0][0][i],mesh->y_1[0][j][0]);
      state->uc_11[0][j][i] = uvc.x;
      state->vc_11[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_11[0][0][0], &Proc.FieldReq[async_h_11], async_h_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_11[0][0][0], &Proc.FieldReq[async_u_11], async_u_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_11[0][0][0], &Proc.FieldReq[async_v_11], async_v_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_11[0][0][0], &Proc.FieldReq[async_jab_11], async_jab_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_11[0][0][0], &Proc.FieldReq[async_uc_11], async_uc_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_11[0][0][0], &Proc.FieldReq[async_vc_11], async_vc_11, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state->v_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_2[0][0][i],mesh->y_1[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_2[0][0][i],mesh->y_1[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_21[0][j][i] = (jab * h);
      state->u_21[0][j][i] = uv.x;
      state->v_21[0][j][i] = uv.y;
      state->jab_21[0][j][i] = jab;
      alpha = (mesh->x_2[0][0][i] / 6371220.0);
      beta = (mesh->y_1[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_21[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_21[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_21[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_21[0][j][i],state->v_21[0][j][i],mesh->x_2[0][0][i],mesh->y_1[0][j][0]);
      state->uc_21[0][j][i] = uvc.x;
      state->vc_21[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_21Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_21[0][0][0], &Proc.FieldReq[async_h_21], async_h_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_21[0][0][0], &Proc.FieldReq[async_u_21], async_u_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_21[0][0][0], &Proc.FieldReq[async_v_21], async_v_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_21[0][0][0], &Proc.FieldReq[async_jab_21], async_jab_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_21[0][0][0], &Proc.FieldReq[async_uc_21], async_uc_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_21[0][0][0], &Proc.FieldReq[async_vc_21], async_vc_21, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_3[0][0][i],mesh->y_1[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_3[0][0][i],mesh->y_1[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_31[0][j][i] = (jab * h);
      state->u_31[0][j][i] = uv.x;
      state->v_31[0][j][i] = uv.y;
      state->jab_31[0][j][i] = jab;
      alpha = (mesh->x_3[0][0][i] / 6371220.0);
      beta = (mesh->y_1[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_31[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_31[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_31[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_31[0][j][i],state->v_31[0][j][i],mesh->x_3[0][0][i],mesh->y_1[0][j][0]);
      state->uc_31[0][j][i] = uvc.x;
      state->vc_31[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_31Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_31[0][0][0], &Proc.FieldReq[async_h_31], async_h_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_31[0][0][0], &Proc.FieldReq[async_u_31], async_u_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_31[0][0][0], &Proc.FieldReq[async_v_31], async_v_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_31[0][0][0], &Proc.FieldReq[async_jab_31], async_jab_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab5_31[0][0][0], &Proc.FieldReq[async_jab5_31], async_jab5_31, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_31[0][0][0], &Proc.FieldReq[async_uc_31], async_uc_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_31[0][0][0], &Proc.FieldReq[async_vc_31], async_vc_31, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state->v_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_1[0][0][i],mesh->y_2[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_1[0][0][i],mesh->y_2[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_12[0][j][i] = (jab * h);
      state->u_12[0][j][i] = uv.x;
      state->v_12[0][j][i] = uv.y;
      state->jab_12[0][j][i] = jab;
      alpha = (mesh->x_1[0][0][i] / 6371220.0);
      beta = (mesh->y_2[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_12[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_12[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_12[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_12[0][j][i],state->v_12[0][j][i],mesh->x_1[0][0][i],mesh->y_2[0][j][0]);
      state->uc_12[0][j][i] = uvc.x;
      state->vc_12[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_12[0][0][0], &Proc.FieldReq[async_h_12], async_h_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_12[0][0][0], &Proc.FieldReq[async_u_12], async_u_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_12[0][0][0], &Proc.FieldReq[async_v_12], async_v_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_12[0][0][0], &Proc.FieldReq[async_jab_12], async_jab_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_12[0][0][0], &Proc.FieldReq[async_uc_12], async_uc_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_12[0][0][0], &Proc.FieldReq[async_vc_12], async_vc_12, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state->v_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_2[0][0][i],mesh->y_2[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_2[0][0][i],mesh->y_2[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_22[0][j][i] = (jab * h);
      state->u_22[0][j][i] = uv.x;
      state->v_22[0][j][i] = uv.y;
      state->jab_22[0][j][i] = jab;
      alpha = (mesh->x_2[0][0][i] / 6371220.0);
      beta = (mesh->y_2[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_22[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_22[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_22[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_22[0][j][i],state->v_22[0][j][i],mesh->x_2[0][0][i],mesh->y_2[0][j][0]);
      state->uc_22[0][j][i] = uvc.x;
      state->vc_22[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_22Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_22[0][0][0], &Proc.FieldReq[async_h_22], async_h_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_22[0][0][0], &Proc.FieldReq[async_u_22], async_u_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_22[0][0][0], &Proc.FieldReq[async_v_22], async_v_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_22[0][0][0], &Proc.FieldReq[async_jab_22], async_jab_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_22[0][0][0], &Proc.FieldReq[async_uc_22], async_uc_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_22[0][0][0], &Proc.FieldReq[async_vc_22], async_vc_22, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state->v_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_3[0][0][i],mesh->y_2[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_3[0][0][i],mesh->y_2[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_32[0][j][i] = (jab * h);
      state->u_32[0][j][i] = uv.x;
      state->v_32[0][j][i] = uv.y;
      state->jab_32[0][j][i] = jab;
      alpha = (mesh->x_3[0][0][i] / 6371220.0);
      beta = (mesh->y_2[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_32[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_32[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_32[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_32[0][j][i],state->v_32[0][j][i],mesh->x_3[0][0][i],mesh->y_2[0][j][0]);
      state->uc_32[0][j][i] = uvc.x;
      state->vc_32[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_32Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_32[0][0][0], &Proc.FieldReq[async_h_32], async_h_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_32[0][0][0], &Proc.FieldReq[async_u_32], async_u_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_32[0][0][0], &Proc.FieldReq[async_v_32], async_v_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_32[0][0][0], &Proc.FieldReq[async_jab_32], async_jab_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab5_32[0][0][0], &Proc.FieldReq[async_jab5_32], async_jab5_32, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_32[0][0][0], &Proc.FieldReq[async_uc_32], async_uc_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_32[0][0][0], &Proc.FieldReq[async_vc_32], async_vc_32, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_1[0][0][i],mesh->y_3[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_1[0][0][i],mesh->y_3[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_13[0][j][i] = (jab * h);
      state->u_13[0][j][i] = uv.x;
      state->v_13[0][j][i] = uv.y;
      state->jab_13[0][j][i] = jab;
      alpha = (mesh->x_1[0][0][i] / 6371220.0);
      beta = (mesh->y_3[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_13[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_13[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_13[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_13[0][j][i],state->v_13[0][j][i],mesh->x_1[0][0][i],mesh->y_3[0][j][0]);
      state->uc_13[0][j][i] = uvc.x;
      state->vc_13[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_13[0][0][0], &Proc.FieldReq[async_h_13], async_h_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_13[0][0][0], &Proc.FieldReq[async_u_13], async_u_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_13[0][0][0], &Proc.FieldReq[async_v_13], async_v_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_13[0][0][0], &Proc.FieldReq[async_jab_13], async_jab_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab7_13[0][0][0], &Proc.FieldReq[async_jab7_13], async_jab7_13, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_13[0][0][0], &Proc.FieldReq[async_uc_13], async_uc_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_13[0][0][0], &Proc.FieldReq[async_vc_13], async_vc_13, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state->v_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_2[0][0][i],mesh->y_3[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_2[0][0][i],mesh->y_3[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_23[0][j][i] = (jab * h);
      state->u_23[0][j][i] = uv.x;
      state->v_23[0][j][i] = uv.y;
      state->jab_23[0][j][i] = jab;
      alpha = (mesh->x_2[0][0][i] / 6371220.0);
      beta = (mesh->y_3[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_23[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_23[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_23[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_23[0][j][i],state->v_23[0][j][i],mesh->x_2[0][0][i],mesh->y_3[0][j][0]);
      state->uc_23[0][j][i] = uvc.x;
      state->vc_23[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_23Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_23[0][0][0], &Proc.FieldReq[async_h_23], async_h_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_23[0][0][0], &Proc.FieldReq[async_u_23], async_u_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_23[0][0][0], &Proc.FieldReq[async_v_23], async_v_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_23[0][0][0], &Proc.FieldReq[async_jab_23], async_jab_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab7_23[0][0][0], &Proc.FieldReq[async_jab7_23], async_jab7_23, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_23[0][0][0], &Proc.FieldReq[async_uc_23], async_uc_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_23[0][0][0], &Proc.FieldReq[async_vc_23], async_vc_23, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ltret = (struct Vector2){0.0,0.0};
      ltret = pprop2sp(pypanel,mesh->x_3[0][0][i],mesh->y_3[0][j][0]);
      lamb = ltret.x;
      theta = ltret.y;
      h = computeh(lamb,theta);
      jab = computejab(mesh->x_3[0][0][i],mesh->y_3[0][j][0]);
      uv = (struct Vector2){0.0,0.0};
      uv = computevel(pypanel,lamb,theta);
      state->h_33[0][j][i] = (jab * h);
      state->u_33[0][j][i] = uv.x;
      state->v_33[0][j][i] = uv.y;
      state->jab_33[0][j][i] = jab;
      alpha = (mesh->x_3[0][0][i] / 6371220.0);
      beta = (mesh->y_3[0][j][0] / 6371220.0);
      rho = sqrt(((1.0 + (tan(alpha) * tan(alpha))) + (tan(beta) * tan(beta))));
      state->jab5_33[0][j][i] = ((1.0 + pow(tan(beta),2)) * ((pow(rho,2.0) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab6_33[0][j][i] = ((tan(alpha) * tan(beta)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      state->jab7_33[0][j][i] = ((1.0 + pow(tan(alpha),2.0)) * ((pow(rho,2) * pow(cos(alpha),2)) * pow(cos(beta),2)));
      uvc = (struct Vector2){0.0,0.0};
      uvc = cov2contrav(state->u_33[0][j][i],state->v_33[0][j][i],mesh->x_3[0][0][i],mesh->y_3[0][j][0]);
      state->uc_33[0][j][i] = uvc.x;
      state->vc_33[0][j][i] = uvc.y;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  StateInitialh_33Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state->h_33[0][0][0], &Proc.FieldReq[async_h_33], async_h_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->u_33[0][0][0], &Proc.FieldReq[async_u_33], async_u_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->v_33[0][0][0], &Proc.FieldReq[async_v_33], async_v_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab_33[0][0][0], &Proc.FieldReq[async_jab_33], async_jab_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->jab5_33[0][0][0], &Proc.FieldReq[async_jab5_33], async_jab5_33, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state->jab7_33[0][0][0], &Proc.FieldReq[async_jab7_33], async_jab7_33, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->uc_33[0][0][0], &Proc.FieldReq[async_uc_33], async_uc_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state->vc_33[0][0][0], &Proc.FieldReq[async_vc_33], async_vc_33, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void updateState_0(struct HybridStateField* state_old, struct HybridStateField* state_new)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state_old->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state_old->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state_old->v_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_11], async_jab_11, state_old->jab_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state_old->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state_old->vc_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_0_info updateState_0_0_para;
  updateState_0_0_para.lz = 1 - 0;
  updateState_0_0_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_0_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_0_para.oz = 0;
  updateState_0_0_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_0_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_0_para.hx = 1;
  updateState_0_0_para.hy = 1;
  updateState_0_0_para.hz = 0;
  updateState_0_0_para.bx = 128;
  updateState_0_0_para.by = 4;
  updateState_0_0_para.bz = 1;
  updateState_0_0_para.mx = 1;
  updateState_0_0_para.my = 4;
  updateState_0_0_para.mz = 1;
  updateState_0_0_para.state_oldh_11 = &state_old->h_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldu_11 = &state_old->u_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldv_11 = &state_old->v_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldjab_11 = &state_old->jab_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldjab5_11 = &state_old->jab5_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldjab6_11 = &state_old->jab6_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldjab7_11 = &state_old->jab7_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_olduc_11 = &state_old->uc_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_oldvc_11 = &state_old->vc_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newh_11 = &state_new->h_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newu_11 = &state_new->u_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newv_11 = &state_new->v_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newjab_11 = &state_new->jab_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newjab5_11 = &state_new->jab5_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newjab6_11 = &state_new->jab6_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newjab7_11 = &state_new->jab7_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newuc_11 = &state_new->uc_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  updateState_0_0_para.state_newvc_11 = &state_new->vc_11[0][Proc.lat_hw -1+updateState_0_0_para.oy][Proc.lon_hw -1+updateState_0_0_para.ox];
  
  athread_spawn(updateState_0_0, &updateState_0_0_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_11[0][0][0], &Proc.FieldReq[async_h_11], async_h_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_11[0][0][0], &Proc.FieldReq[async_u_11], async_u_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_11[0][0][0], &Proc.FieldReq[async_v_11], async_v_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_11[0][0][0], &Proc.FieldReq[async_jab_11], async_jab_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_11[0][0][0], &Proc.FieldReq[async_uc_11], async_uc_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_11[0][0][0], &Proc.FieldReq[async_vc_11], async_vc_11, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state_old->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state_old->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state_old->v_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_12], async_jab_12, state_old->jab_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state_old->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state_old->vc_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_1_info updateState_0_1_para;
  updateState_0_1_para.lz = 1 - 0;
  updateState_0_1_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_1_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_1_para.oz = 0;
  updateState_0_1_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_1_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_1_para.hx = 1;
  updateState_0_1_para.hy = 1;
  updateState_0_1_para.hz = 0;
  updateState_0_1_para.bx = 128;
  updateState_0_1_para.by = 4;
  updateState_0_1_para.bz = 1;
  updateState_0_1_para.mx = 1;
  updateState_0_1_para.my = 4;
  updateState_0_1_para.mz = 1;
  updateState_0_1_para.state_oldh_12 = &state_old->h_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldu_12 = &state_old->u_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldv_12 = &state_old->v_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldjab_12 = &state_old->jab_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldjab5_12 = &state_old->jab5_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldjab6_12 = &state_old->jab6_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldjab7_12 = &state_old->jab7_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_olduc_12 = &state_old->uc_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_oldvc_12 = &state_old->vc_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newh_12 = &state_new->h_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newu_12 = &state_new->u_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newv_12 = &state_new->v_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newjab_12 = &state_new->jab_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newjab5_12 = &state_new->jab5_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newjab6_12 = &state_new->jab6_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newjab7_12 = &state_new->jab7_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newuc_12 = &state_new->uc_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  updateState_0_1_para.state_newvc_12 = &state_new->vc_12[0][Proc.lat_hw -1+updateState_0_1_para.oy][Proc.lon_hw -1+updateState_0_1_para.ox];
  
  athread_spawn(updateState_0_1, &updateState_0_1_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_12[0][0][0], &Proc.FieldReq[async_h_12], async_h_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_12[0][0][0], &Proc.FieldReq[async_u_12], async_u_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_12[0][0][0], &Proc.FieldReq[async_v_12], async_v_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_12[0][0][0], &Proc.FieldReq[async_jab_12], async_jab_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_12[0][0][0], &Proc.FieldReq[async_uc_12], async_uc_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_12[0][0][0], &Proc.FieldReq[async_vc_12], async_vc_12, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state_old->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state_old->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state_old->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state_old->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_13], async_jab7_13, state_old->jab7_13, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state_old->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state_old->vc_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_2_info updateState_0_2_para;
  updateState_0_2_para.lz = 1 - 0;
  updateState_0_2_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_2_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_2_para.oz = 0;
  updateState_0_2_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_2_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_2_para.hx = 1;
  updateState_0_2_para.hy = 1;
  updateState_0_2_para.hz = 0;
  updateState_0_2_para.bx = 128;
  updateState_0_2_para.by = 4;
  updateState_0_2_para.bz = 1;
  updateState_0_2_para.mx = 1;
  updateState_0_2_para.my = 4;
  updateState_0_2_para.mz = 1;
  updateState_0_2_para.state_oldh_13 = &state_old->h_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldu_13 = &state_old->u_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldv_13 = &state_old->v_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldjab_13 = &state_old->jab_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldjab5_13 = &state_old->jab5_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldjab6_13 = &state_old->jab6_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldjab7_13 = &state_old->jab7_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_olduc_13 = &state_old->uc_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_oldvc_13 = &state_old->vc_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newh_13 = &state_new->h_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newu_13 = &state_new->u_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newv_13 = &state_new->v_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newjab_13 = &state_new->jab_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newjab5_13 = &state_new->jab5_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newjab6_13 = &state_new->jab6_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newjab7_13 = &state_new->jab7_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newuc_13 = &state_new->uc_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  updateState_0_2_para.state_newvc_13 = &state_new->vc_13[0][Proc.lat_hw -1+updateState_0_2_para.oy][Proc.lon_hw -1+updateState_0_2_para.ox];
  
  athread_spawn(updateState_0_2, &updateState_0_2_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_13[0][0][0], &Proc.FieldReq[async_h_13], async_h_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_13[0][0][0], &Proc.FieldReq[async_u_13], async_u_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_13[0][0][0], &Proc.FieldReq[async_v_13], async_v_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_13[0][0][0], &Proc.FieldReq[async_jab_13], async_jab_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab7_13[0][0][0], &Proc.FieldReq[async_jab7_13], async_jab7_13, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_13[0][0][0], &Proc.FieldReq[async_uc_13], async_uc_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_13[0][0][0], &Proc.FieldReq[async_vc_13], async_vc_13, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state_old->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state_old->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state_old->v_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_21], async_jab_21, state_old->jab_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state_old->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state_old->vc_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_3_info updateState_0_3_para;
  updateState_0_3_para.lz = 1 - 0;
  updateState_0_3_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_3_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_3_para.oz = 0;
  updateState_0_3_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_3_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_3_para.hx = 1;
  updateState_0_3_para.hy = 1;
  updateState_0_3_para.hz = 0;
  updateState_0_3_para.bx = 128;
  updateState_0_3_para.by = 4;
  updateState_0_3_para.bz = 1;
  updateState_0_3_para.mx = 1;
  updateState_0_3_para.my = 4;
  updateState_0_3_para.mz = 1;
  updateState_0_3_para.state_oldh_21 = &state_old->h_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldu_21 = &state_old->u_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldv_21 = &state_old->v_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldjab_21 = &state_old->jab_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldjab5_21 = &state_old->jab5_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldjab6_21 = &state_old->jab6_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldjab7_21 = &state_old->jab7_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_olduc_21 = &state_old->uc_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_oldvc_21 = &state_old->vc_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newh_21 = &state_new->h_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newu_21 = &state_new->u_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newv_21 = &state_new->v_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newjab_21 = &state_new->jab_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newjab5_21 = &state_new->jab5_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newjab6_21 = &state_new->jab6_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newjab7_21 = &state_new->jab7_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newuc_21 = &state_new->uc_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  updateState_0_3_para.state_newvc_21 = &state_new->vc_21[0][Proc.lat_hw -1+updateState_0_3_para.oy][Proc.lon_hw -1+updateState_0_3_para.ox];
  
  athread_spawn(updateState_0_3, &updateState_0_3_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_21Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_21[0][0][0], &Proc.FieldReq[async_h_21], async_h_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_21[0][0][0], &Proc.FieldReq[async_u_21], async_u_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_21[0][0][0], &Proc.FieldReq[async_v_21], async_v_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_21[0][0][0], &Proc.FieldReq[async_jab_21], async_jab_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_21[0][0][0], &Proc.FieldReq[async_uc_21], async_uc_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_21[0][0][0], &Proc.FieldReq[async_vc_21], async_vc_21, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state_old->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state_old->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state_old->v_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_22], async_jab_22, state_old->jab_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state_old->uc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state_old->vc_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_4_info updateState_0_4_para;
  updateState_0_4_para.lz = 1 - 0;
  updateState_0_4_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_4_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_4_para.oz = 0;
  updateState_0_4_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_4_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_4_para.hx = 1;
  updateState_0_4_para.hy = 1;
  updateState_0_4_para.hz = 0;
  updateState_0_4_para.bx = 128;
  updateState_0_4_para.by = 4;
  updateState_0_4_para.bz = 1;
  updateState_0_4_para.mx = 1;
  updateState_0_4_para.my = 4;
  updateState_0_4_para.mz = 1;
  updateState_0_4_para.state_oldh_22 = &state_old->h_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldu_22 = &state_old->u_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldv_22 = &state_old->v_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldjab_22 = &state_old->jab_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldjab5_22 = &state_old->jab5_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldjab6_22 = &state_old->jab6_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldjab7_22 = &state_old->jab7_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_olduc_22 = &state_old->uc_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_oldvc_22 = &state_old->vc_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newh_22 = &state_new->h_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newu_22 = &state_new->u_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newv_22 = &state_new->v_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newjab_22 = &state_new->jab_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newjab5_22 = &state_new->jab5_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newjab6_22 = &state_new->jab6_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newjab7_22 = &state_new->jab7_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newuc_22 = &state_new->uc_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  updateState_0_4_para.state_newvc_22 = &state_new->vc_22[0][Proc.lat_hw -1+updateState_0_4_para.oy][Proc.lon_hw -1+updateState_0_4_para.ox];
  
  athread_spawn(updateState_0_4, &updateState_0_4_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_22Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_22[0][0][0], &Proc.FieldReq[async_h_22], async_h_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_22[0][0][0], &Proc.FieldReq[async_u_22], async_u_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_22[0][0][0], &Proc.FieldReq[async_v_22], async_v_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_22[0][0][0], &Proc.FieldReq[async_jab_22], async_jab_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_22[0][0][0], &Proc.FieldReq[async_uc_22], async_uc_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_22[0][0][0], &Proc.FieldReq[async_vc_22], async_vc_22, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state_old->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state_old->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state_old->v_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state_old->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_23], async_jab7_23, state_old->jab7_23, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state_old->uc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state_old->vc_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_5_info updateState_0_5_para;
  updateState_0_5_para.lz = 1 - 0;
  updateState_0_5_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_5_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_5_para.oz = 0;
  updateState_0_5_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_5_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_5_para.hx = 1;
  updateState_0_5_para.hy = 1;
  updateState_0_5_para.hz = 0;
  updateState_0_5_para.bx = 128;
  updateState_0_5_para.by = 4;
  updateState_0_5_para.bz = 1;
  updateState_0_5_para.mx = 1;
  updateState_0_5_para.my = 4;
  updateState_0_5_para.mz = 1;
  updateState_0_5_para.state_oldh_23 = &state_old->h_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldu_23 = &state_old->u_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldv_23 = &state_old->v_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldjab_23 = &state_old->jab_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldjab5_23 = &state_old->jab5_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldjab6_23 = &state_old->jab6_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldjab7_23 = &state_old->jab7_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_olduc_23 = &state_old->uc_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_oldvc_23 = &state_old->vc_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newh_23 = &state_new->h_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newu_23 = &state_new->u_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newv_23 = &state_new->v_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newjab_23 = &state_new->jab_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newjab5_23 = &state_new->jab5_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newjab6_23 = &state_new->jab6_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newjab7_23 = &state_new->jab7_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newuc_23 = &state_new->uc_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  updateState_0_5_para.state_newvc_23 = &state_new->vc_23[0][Proc.lat_hw -1+updateState_0_5_para.oy][Proc.lon_hw -1+updateState_0_5_para.ox];
  
  athread_spawn(updateState_0_5, &updateState_0_5_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_23Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_23[0][0][0], &Proc.FieldReq[async_h_23], async_h_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_23[0][0][0], &Proc.FieldReq[async_u_23], async_u_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_23[0][0][0], &Proc.FieldReq[async_v_23], async_v_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_23[0][0][0], &Proc.FieldReq[async_jab_23], async_jab_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab7_23[0][0][0], &Proc.FieldReq[async_jab7_23], async_jab7_23, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_23[0][0][0], &Proc.FieldReq[async_uc_23], async_uc_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_23[0][0][0], &Proc.FieldReq[async_vc_23], async_vc_23, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state_old->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state_old->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state_old->v_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state_old->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_31], async_jab5_31, state_old->jab5_31, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state_old->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state_old->vc_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_6_info updateState_0_6_para;
  updateState_0_6_para.lz = 1 - 0;
  updateState_0_6_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_6_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_6_para.oz = 0;
  updateState_0_6_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_6_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_6_para.hx = 1;
  updateState_0_6_para.hy = 1;
  updateState_0_6_para.hz = 0;
  updateState_0_6_para.bx = 128;
  updateState_0_6_para.by = 4;
  updateState_0_6_para.bz = 1;
  updateState_0_6_para.mx = 1;
  updateState_0_6_para.my = 4;
  updateState_0_6_para.mz = 1;
  updateState_0_6_para.state_oldh_31 = &state_old->h_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldu_31 = &state_old->u_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldv_31 = &state_old->v_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldjab_31 = &state_old->jab_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldjab5_31 = &state_old->jab5_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldjab6_31 = &state_old->jab6_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldjab7_31 = &state_old->jab7_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_olduc_31 = &state_old->uc_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_oldvc_31 = &state_old->vc_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newh_31 = &state_new->h_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newu_31 = &state_new->u_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newv_31 = &state_new->v_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newjab_31 = &state_new->jab_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newjab5_31 = &state_new->jab5_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newjab6_31 = &state_new->jab6_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newjab7_31 = &state_new->jab7_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newuc_31 = &state_new->uc_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  updateState_0_6_para.state_newvc_31 = &state_new->vc_31[0][Proc.lat_hw -1+updateState_0_6_para.oy][Proc.lon_hw -1+updateState_0_6_para.ox];
  
  athread_spawn(updateState_0_6, &updateState_0_6_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_31Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_31[0][0][0], &Proc.FieldReq[async_h_31], async_h_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_31[0][0][0], &Proc.FieldReq[async_u_31], async_u_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_31[0][0][0], &Proc.FieldReq[async_v_31], async_v_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_31[0][0][0], &Proc.FieldReq[async_jab_31], async_jab_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab5_31[0][0][0], &Proc.FieldReq[async_jab5_31], async_jab5_31, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_31[0][0][0], &Proc.FieldReq[async_uc_31], async_uc_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_31[0][0][0], &Proc.FieldReq[async_vc_31], async_vc_31, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state_old->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state_old->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state_old->v_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state_old->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_32], async_jab5_32, state_old->jab5_32, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state_old->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state_old->vc_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_7_info updateState_0_7_para;
  updateState_0_7_para.lz = 1 - 0;
  updateState_0_7_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_7_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_7_para.oz = 0;
  updateState_0_7_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_7_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_7_para.hx = 1;
  updateState_0_7_para.hy = 1;
  updateState_0_7_para.hz = 0;
  updateState_0_7_para.bx = 128;
  updateState_0_7_para.by = 4;
  updateState_0_7_para.bz = 1;
  updateState_0_7_para.mx = 1;
  updateState_0_7_para.my = 4;
  updateState_0_7_para.mz = 1;
  updateState_0_7_para.state_oldh_32 = &state_old->h_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldu_32 = &state_old->u_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldv_32 = &state_old->v_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldjab_32 = &state_old->jab_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldjab5_32 = &state_old->jab5_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldjab6_32 = &state_old->jab6_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldjab7_32 = &state_old->jab7_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_olduc_32 = &state_old->uc_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_oldvc_32 = &state_old->vc_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newh_32 = &state_new->h_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newu_32 = &state_new->u_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newv_32 = &state_new->v_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newjab_32 = &state_new->jab_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newjab5_32 = &state_new->jab5_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newjab6_32 = &state_new->jab6_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newjab7_32 = &state_new->jab7_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newuc_32 = &state_new->uc_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  updateState_0_7_para.state_newvc_32 = &state_new->vc_32[0][Proc.lat_hw -1+updateState_0_7_para.oy][Proc.lon_hw -1+updateState_0_7_para.ox];
  
  athread_spawn(updateState_0_7, &updateState_0_7_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_32Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_32[0][0][0], &Proc.FieldReq[async_h_32], async_h_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_32[0][0][0], &Proc.FieldReq[async_u_32], async_u_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_32[0][0][0], &Proc.FieldReq[async_v_32], async_v_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_32[0][0][0], &Proc.FieldReq[async_jab_32], async_jab_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab5_32[0][0][0], &Proc.FieldReq[async_jab5_32], async_jab5_32, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_32[0][0][0], &Proc.FieldReq[async_uc_32], async_uc_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_32[0][0][0], &Proc.FieldReq[async_vc_32], async_vc_32, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state_old->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state_old->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state_old->v_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state_old->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_33], async_jab5_33, state_old->jab5_33, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_33], async_jab7_33, state_old->jab7_33, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state_old->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state_old->vc_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateState_0_8_info updateState_0_8_para;
  updateState_0_8_para.lz = 1 - 0;
  updateState_0_8_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 1);
  updateState_0_8_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 1);
  updateState_0_8_para.oz = 0;
  updateState_0_8_para.oy = MAX(Proc.lat_beg, 1) - Proc.lat_beg;
  updateState_0_8_para.ox = MAX(Proc.lon_beg, 1) - Proc.lon_beg;
  updateState_0_8_para.hx = 1;
  updateState_0_8_para.hy = 1;
  updateState_0_8_para.hz = 0;
  updateState_0_8_para.bx = 128;
  updateState_0_8_para.by = 4;
  updateState_0_8_para.bz = 1;
  updateState_0_8_para.mx = 1;
  updateState_0_8_para.my = 4;
  updateState_0_8_para.mz = 1;
  updateState_0_8_para.state_oldh_33 = &state_old->h_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldu_33 = &state_old->u_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldv_33 = &state_old->v_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldjab_33 = &state_old->jab_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldjab5_33 = &state_old->jab5_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldjab6_33 = &state_old->jab6_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldjab7_33 = &state_old->jab7_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_olduc_33 = &state_old->uc_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_oldvc_33 = &state_old->vc_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newh_33 = &state_new->h_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newu_33 = &state_new->u_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newv_33 = &state_new->v_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newjab_33 = &state_new->jab_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newjab5_33 = &state_new->jab5_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newjab6_33 = &state_new->jab6_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newjab7_33 = &state_new->jab7_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newuc_33 = &state_new->uc_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  updateState_0_8_para.state_newvc_33 = &state_new->vc_33[0][Proc.lat_hw -1+updateState_0_8_para.oy][Proc.lon_hw -1+updateState_0_8_para.ox];
  
  athread_spawn(updateState_0_8, &updateState_0_8_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateStateh_33Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_33[0][0][0], &Proc.FieldReq[async_h_33], async_h_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_33[0][0][0], &Proc.FieldReq[async_u_33], async_u_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_33[0][0][0], &Proc.FieldReq[async_v_33], async_v_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab_33[0][0][0], &Proc.FieldReq[async_jab_33], async_jab_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab5_33[0][0][0], &Proc.FieldReq[async_jab5_33], async_jab5_33, true, true, true, false, false, false);
  UpdateHaloCS_2d_D(Proc, &state_new->jab7_33[0][0][0], &Proc.FieldReq[async_jab7_33], async_jab7_33, true, true, false, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->uc_33[0][0][0], &Proc.FieldReq[async_uc_33], async_uc_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->vc_33[0][0][0], &Proc.FieldReq[async_vc_33], async_vc_33, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void mcvupdateXY_0(struct HybridStateField* state, struct HybridTendField* tend, struct HybridMeshField* mesh)
{
  double ulocal, dql, dqr, dfl, dfr, sl, sr, sc, vlocal;
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_31], async_jab5_31, state->jab5_31, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_11], async_jab_11, state->jab_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_30_info updateX_0_30_para;
  updateX_0_30_para.lz = 1 - 0;
  updateX_0_30_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_30_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_30_para.oz = 0;
  updateX_0_30_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_30_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_30_para.hx = 1;
  updateX_0_30_para.hy = 1;
  updateX_0_30_para.hz = 0;
  updateX_0_30_para.bx = 128;
  updateX_0_30_para.by = 4;
  updateX_0_30_para.bz = 1;
  updateX_0_30_para.mx = 1;
  updateX_0_30_para.my = 4;
  updateX_0_30_para.mz = 1;
  updateX_0_30_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statejab5_31 = &state->jab5_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.stateu_31 = &state->u_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statev_31 = &state->v_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statejab_11 = &state->jab_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.stateu_11 = &state->u_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statev_11 = &state->v_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  updateX_0_30_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_30_para.oy][Proc.lon_hw -1+updateX_0_30_para.ox];
  
  athread_spawn(updateX_0_30, &updateX_0_30_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_31], async_jab5_31, state->jab5_31, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_31_info updateX_0_31_para;
  updateX_0_31_para.lz = 1 - 0;
  updateX_0_31_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_31_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_31_para.oz = 0;
  updateX_0_31_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_31_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_31_para.hx = 1;
  updateX_0_31_para.hy = 1;
  updateX_0_31_para.hz = 0;
  updateX_0_31_para.bx = 128;
  updateX_0_31_para.by = 8;
  updateX_0_31_para.bz = 1;
  updateX_0_31_para.mx = 1;
  updateX_0_31_para.my = 4;
  updateX_0_31_para.mz = 1;
  updateX_0_31_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.statejab5_31 = &state->jab5_31[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.stateh_21 = &state->h_21[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  updateX_0_31_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateX_0_31_para.oy][Proc.lon_hw -1+updateX_0_31_para.ox];
  
  athread_spawn(updateX_0_31, &updateX_0_31_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_31], async_jab5_31, state->jab5_31, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_11], async_jab_11, state->jab_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_21], async_jab_21, state->jab_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state->v_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_32_info updateX_0_32_para;
  updateX_0_32_para.lz = 1 - 0;
  updateX_0_32_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_32_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_32_para.oz = 0;
  updateX_0_32_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_32_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_32_para.hx = 1;
  updateX_0_32_para.hy = 1;
  updateX_0_32_para.hz = 0;
  updateX_0_32_para.bx = 128;
  updateX_0_32_para.by = 4;
  updateX_0_32_para.bz = 1;
  updateX_0_32_para.mx = 1;
  updateX_0_32_para.my = 4;
  updateX_0_32_para.mz = 1;
  updateX_0_32_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statejab5_31 = &state->jab5_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateu_11 = &state->u_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateu_21 = &state->u_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateu_31 = &state->u_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statejab_11 = &state->jab_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statev_11 = &state->v_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateh_21 = &state->h_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statejab_21 = &state->jab_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statev_21 = &state->v_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statev_31 = &state->v_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  updateX_0_32_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateX_0_32_para.oy][Proc.lon_hw -1+updateX_0_32_para.ox];
  
  athread_spawn(updateX_0_32, &updateX_0_32_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_31], async_jab5_31, state->jab5_31, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state->v_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_33_info updateX_0_33_para;
  updateX_0_33_para.lz = 1 - 0;
  updateX_0_33_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_33_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_33_para.oz = 0;
  updateX_0_33_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_33_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_33_para.hx = 1;
  updateX_0_33_para.hy = 1;
  updateX_0_33_para.hz = 0;
  updateX_0_33_para.bx = 128;
  updateX_0_33_para.by = 8;
  updateX_0_33_para.bz = 1;
  updateX_0_33_para.mx = 1;
  updateX_0_33_para.my = 4;
  updateX_0_33_para.mz = 1;
  updateX_0_33_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.statejab5_31 = &state->jab5_31[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.statev_11 = &state->v_11[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.statev_21 = &state->v_21[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.statev_31 = &state->v_31[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  updateX_0_33_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_33_para.oy][Proc.lon_hw -1+updateX_0_33_para.ox];
  
  athread_spawn(updateX_0_33, &updateX_0_33_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_2], async_qv_2, tend->qv_2, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_1], async_qv_1, tend->qv_1, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_34_info updateX_0_34_para;
  updateX_0_34_para.lz = 1 - 0;
  updateX_0_34_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_34_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_34_para.oz = 0;
  updateX_0_34_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_34_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_34_para.hx = 1;
  updateX_0_34_para.hy = 1;
  updateX_0_34_para.hz = 0;
  updateX_0_34_para.bx = 128;
  updateX_0_34_para.by = 4;
  updateX_0_34_para.bz = 1;
  updateX_0_34_para.mx = 1;
  updateX_0_34_para.my = 4;
  updateX_0_34_para.mz = 1;
  updateX_0_34_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxh_11 = &tend->dxh_11[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxh_21 = &tend->dxh_21[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxh_31 = &tend->dxh_31[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxu_11 = &tend->dxu_11[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxu_21 = &tend->dxu_21[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxu_31 = &tend->dxu_31[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxv_11 = &tend->dxv_11[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxv_21 = &tend->dxv_21[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  updateX_0_34_para.tenddxv_31 = &tend->dxv_31[0][Proc.lat_hw -1+updateX_0_34_para.oy][Proc.lon_hw -1+updateX_0_34_para.ox];
  
  athread_spawn(updateX_0_34, &updateX_0_34_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXdxh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_32], async_jab5_32, state->jab5_32, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state->v_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_12], async_jab_12, state->jab_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state->v_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_35_info updateX_0_35_para;
  updateX_0_35_para.lz = 1 - 0;
  updateX_0_35_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_35_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_35_para.oz = 0;
  updateX_0_35_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_35_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_35_para.hx = 1;
  updateX_0_35_para.hy = 1;
  updateX_0_35_para.hz = 0;
  updateX_0_35_para.bx = 128;
  updateX_0_35_para.by = 4;
  updateX_0_35_para.bz = 1;
  updateX_0_35_para.mx = 1;
  updateX_0_35_para.my = 4;
  updateX_0_35_para.mz = 1;
  updateX_0_35_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statejab5_32 = &state->jab5_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statejab_32 = &state->jab_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.stateh_12 = &state->h_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.stateu_32 = &state->u_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statev_32 = &state->v_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statejab_12 = &state->jab_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.stateu_12 = &state->u_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statev_12 = &state->v_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  updateX_0_35_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_35_para.oy][Proc.lon_hw -1+updateX_0_35_para.ox];
  
  athread_spawn(updateX_0_35, &updateX_0_35_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_32], async_jab5_32, state->jab5_32, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state->uc_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_36_info updateX_0_36_para;
  updateX_0_36_para.lz = 1 - 0;
  updateX_0_36_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_36_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_36_para.oz = 0;
  updateX_0_36_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_36_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_36_para.hx = 1;
  updateX_0_36_para.hy = 1;
  updateX_0_36_para.hz = 0;
  updateX_0_36_para.bx = 128;
  updateX_0_36_para.by = 8;
  updateX_0_36_para.bz = 1;
  updateX_0_36_para.mx = 1;
  updateX_0_36_para.my = 4;
  updateX_0_36_para.mz = 1;
  updateX_0_36_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.statejab5_32 = &state->jab5_32[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.statejab_32 = &state->jab_32[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.stateh_12 = &state->h_12[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.stateh_22 = &state->h_22[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.stateuc_22 = &state->uc_22[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  updateX_0_36_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateX_0_36_para.oy][Proc.lon_hw -1+updateX_0_36_para.ox];
  
  athread_spawn(updateX_0_36, &updateX_0_36_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_32], async_jab5_32, state->jab5_32, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_12], async_jab_12, state->jab_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state->v_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_22], async_jab_22, state->jab_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state->uc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state->v_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state->vc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state->v_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_37_info updateX_0_37_para;
  updateX_0_37_para.lz = 1 - 0;
  updateX_0_37_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_37_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_37_para.oz = 0;
  updateX_0_37_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_37_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_37_para.hx = 1;
  updateX_0_37_para.hy = 1;
  updateX_0_37_para.hz = 0;
  updateX_0_37_para.bx = 128;
  updateX_0_37_para.by = 4;
  updateX_0_37_para.bz = 1;
  updateX_0_37_para.mx = 1;
  updateX_0_37_para.my = 4;
  updateX_0_37_para.mz = 1;
  updateX_0_37_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statejab5_32 = &state->jab5_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statejab_32 = &state->jab_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateu_12 = &state->u_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateu_22 = &state->u_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateu_32 = &state->u_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateh_12 = &state->h_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statejab_12 = &state->jab_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statev_12 = &state->v_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateh_22 = &state->h_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statejab_22 = &state->jab_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.stateuc_22 = &state->uc_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statev_22 = &state->v_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statevc_22 = &state->vc_22[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statev_32 = &state->v_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  updateX_0_37_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateX_0_37_para.oy][Proc.lon_hw -1+updateX_0_37_para.ox];
  
  athread_spawn(updateX_0_37, &updateX_0_37_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_32], async_jab5_32, state->jab5_32, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state->v_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state->v_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state->v_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_38_info updateX_0_38_para;
  updateX_0_38_para.lz = 1 - 0;
  updateX_0_38_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_38_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_38_para.oz = 0;
  updateX_0_38_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_38_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_38_para.hx = 1;
  updateX_0_38_para.hy = 1;
  updateX_0_38_para.hz = 0;
  updateX_0_38_para.bx = 128;
  updateX_0_38_para.by = 8;
  updateX_0_38_para.bz = 1;
  updateX_0_38_para.mx = 1;
  updateX_0_38_para.my = 4;
  updateX_0_38_para.mz = 1;
  updateX_0_38_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.statejab5_32 = &state->jab5_32[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.statejab_32 = &state->jab_32[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.statev_12 = &state->v_12[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.statev_22 = &state->v_22[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.statev_32 = &state->v_32[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  updateX_0_38_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_38_para.oy][Proc.lon_hw -1+updateX_0_38_para.ox];
  
  athread_spawn(updateX_0_38, &updateX_0_38_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_2], async_qv_2, tend->qv_2, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_1], async_qv_1, tend->qv_1, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state->vc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state->uc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_39_info updateX_0_39_para;
  updateX_0_39_para.lz = 1 - 0;
  updateX_0_39_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_39_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_39_para.oz = 0;
  updateX_0_39_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_39_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_39_para.hx = 1;
  updateX_0_39_para.hy = 1;
  updateX_0_39_para.hz = 0;
  updateX_0_39_para.bx = 128;
  updateX_0_39_para.by = 4;
  updateX_0_39_para.bz = 1;
  updateX_0_39_para.mx = 1;
  updateX_0_39_para.my = 4;
  updateX_0_39_para.mz = 1;
  updateX_0_39_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.statevc_22 = &state->vc_22[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.stateuc_22 = &state->uc_22[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxh_12 = &tend->dxh_12[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxh_22 = &tend->dxh_22[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxh_32 = &tend->dxh_32[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxu_12 = &tend->dxu_12[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxu_22 = &tend->dxu_22[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxu_32 = &tend->dxu_32[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxv_12 = &tend->dxv_12[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxv_22 = &tend->dxv_22[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  updateX_0_39_para.tenddxv_32 = &tend->dxv_32[0][Proc.lat_hw -1+updateX_0_39_para.oy][Proc.lon_hw -1+updateX_0_39_para.ox];
  
  athread_spawn(updateX_0_39, &updateX_0_39_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXdxh_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_33], async_jab5_33, state->jab5_33, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_40_info updateX_0_40_para;
  updateX_0_40_para.lz = 1 - 0;
  updateX_0_40_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_40_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_40_para.oz = 0;
  updateX_0_40_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_40_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_40_para.hx = 1;
  updateX_0_40_para.hy = 1;
  updateX_0_40_para.hz = 0;
  updateX_0_40_para.bx = 128;
  updateX_0_40_para.by = 4;
  updateX_0_40_para.bz = 1;
  updateX_0_40_para.mx = 1;
  updateX_0_40_para.my = 4;
  updateX_0_40_para.mz = 1;
  updateX_0_40_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statejab5_33 = &state->jab5_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.stateu_33 = &state->u_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statev_33 = &state->v_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.stateu_13 = &state->u_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statev_13 = &state->v_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  updateX_0_40_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_40_para.oy][Proc.lon_hw -1+updateX_0_40_para.ox];
  
  athread_spawn(updateX_0_40, &updateX_0_40_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fh_1[0][0][0], &Proc.FieldReq[async_fh_1], async_fh_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->fu_1[0][0][0], &Proc.FieldReq[async_fu_1], async_fu_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->fv_1[0][0][0], &Proc.FieldReq[async_fv_1], async_fv_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qu_1[0][0][0], &Proc.FieldReq[async_qu_1], async_qu_1, true, true, false, false, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qv_1[0][0][0], &Proc.FieldReq[async_qv_1], async_qv_1, true, true, false, true, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_33], async_jab5_33, state->jab5_33, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_41_info updateX_0_41_para;
  updateX_0_41_para.lz = 1 - 0;
  updateX_0_41_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_41_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_41_para.oz = 0;
  updateX_0_41_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_41_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_41_para.hx = 1;
  updateX_0_41_para.hy = 1;
  updateX_0_41_para.hz = 0;
  updateX_0_41_para.bx = 128;
  updateX_0_41_para.by = 8;
  updateX_0_41_para.bz = 1;
  updateX_0_41_para.mx = 1;
  updateX_0_41_para.my = 4;
  updateX_0_41_para.mz = 1;
  updateX_0_41_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.statejab5_33 = &state->jab5_33[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  updateX_0_41_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateX_0_41_para.oy][Proc.lon_hw -1+updateX_0_41_para.ox];
  
  athread_spawn(updateX_0_41, &updateX_0_41_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fh_2[0][0][0], &Proc.FieldReq[async_fh_2], async_fh_2, true, true, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_33], async_jab5_33, state->jab5_33, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state->v_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_42_info updateX_0_42_para;
  updateX_0_42_para.lz = 1 - 0;
  updateX_0_42_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_42_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_42_para.oz = 0;
  updateX_0_42_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_42_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_42_para.hx = 1;
  updateX_0_42_para.hy = 1;
  updateX_0_42_para.hz = 0;
  updateX_0_42_para.bx = 128;
  updateX_0_42_para.by = 4;
  updateX_0_42_para.bz = 1;
  updateX_0_42_para.mx = 1;
  updateX_0_42_para.my = 4;
  updateX_0_42_para.mz = 1;
  updateX_0_42_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statejab5_33 = &state->jab5_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateu_13 = &state->u_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateu_23 = &state->u_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateu_33 = &state->u_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statev_13 = &state->v_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statejab_23 = &state->jab_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statev_23 = &state->v_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statev_33 = &state->v_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  updateX_0_42_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateX_0_42_para.oy][Proc.lon_hw -1+updateX_0_42_para.ox];
  
  athread_spawn(updateX_0_42, &updateX_0_42_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fu_2[0][0][0], &Proc.FieldReq[async_fu_2], async_fu_2, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qu_2[0][0][0], &Proc.FieldReq[async_qu_2], async_qu_2, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab5_33], async_jab5_33, state->jab5_33, true, true, true, true, false, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state->v_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_43_info updateX_0_43_para;
  updateX_0_43_para.lz = 1 - 0;
  updateX_0_43_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_43_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_43_para.oz = 0;
  updateX_0_43_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_43_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_43_para.hx = 1;
  updateX_0_43_para.hy = 1;
  updateX_0_43_para.hz = 0;
  updateX_0_43_para.bx = 128;
  updateX_0_43_para.by = 8;
  updateX_0_43_para.bz = 1;
  updateX_0_43_para.mx = 1;
  updateX_0_43_para.my = 4;
  updateX_0_43_para.mz = 1;
  updateX_0_43_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.statejab5_33 = &state->jab5_33[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.statev_13 = &state->v_13[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.statev_23 = &state->v_23[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.statev_33 = &state->v_33[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  updateX_0_43_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_43_para.oy][Proc.lon_hw -1+updateX_0_43_para.ox];
  
  athread_spawn(updateX_0_43, &updateX_0_43_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fv_2[0][0][0], &Proc.FieldReq[async_fv_2], async_fv_2, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qv_2[0][0][0], &Proc.FieldReq[async_qv_2], async_qv_2, true, true, false, true, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_2], async_qv_2, tend->qv_2, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qv_1], async_qv_1, tend->qv_1, true, true, true, false, true, false, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateX_0_44_info updateX_0_44_para;
  updateX_0_44_para.lz = 1 - 0;
  updateX_0_44_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateX_0_44_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateX_0_44_para.oz = 0;
  updateX_0_44_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateX_0_44_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateX_0_44_para.hx = 1;
  updateX_0_44_para.hy = 1;
  updateX_0_44_para.hz = 0;
  updateX_0_44_para.bx = 128;
  updateX_0_44_para.by = 4;
  updateX_0_44_para.bz = 1;
  updateX_0_44_para.mx = 1;
  updateX_0_44_para.my = 4;
  updateX_0_44_para.mz = 1;
  updateX_0_44_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxh_13 = &tend->dxh_13[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxh_23 = &tend->dxh_23[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxh_33 = &tend->dxh_33[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxu_13 = &tend->dxu_13[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxu_23 = &tend->dxu_23[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxu_33 = &tend->dxu_33[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxv_13 = &tend->dxv_13[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxv_23 = &tend->dxv_23[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  updateX_0_44_para.tenddxv_33 = &tend->dxv_33[0][Proc.lat_hw -1+updateX_0_44_para.oy][Proc.lon_hw -1+updateX_0_44_para.ox];
  
  athread_spawn(updateX_0_44, &updateX_0_44_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateXdxh_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_13], async_jab7_13, state->jab7_13, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_11], async_jab_11, state->jab_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_30_info updateY_0_30_para;
  updateY_0_30_para.lz = 1 - 0;
  updateY_0_30_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_30_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_30_para.oz = 0;
  updateY_0_30_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_30_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_30_para.hx = 1;
  updateY_0_30_para.hy = 1;
  updateY_0_30_para.hz = 0;
  updateY_0_30_para.bx = 128;
  updateY_0_30_para.by = 4;
  updateY_0_30_para.bz = 1;
  updateY_0_30_para.mx = 1;
  updateY_0_30_para.my = 4;
  updateY_0_30_para.mz = 1;
  updateY_0_30_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statejab7_13 = &state->jab7_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateu_13 = &state->u_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateu_11 = &state->u_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statev_13 = &state->v_13[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statejab_11 = &state->jab_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.statev_11 = &state->v_11[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  updateY_0_30_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateY_0_30_para.oy][Proc.lon_hw -1+updateY_0_30_para.ox];
  
  athread_spawn(updateY_0_30, &updateY_0_30_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_13], async_jab7_13, state->jab7_13, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_31_info updateY_0_31_para;
  updateY_0_31_para.lz = 1 - 0;
  updateY_0_31_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_31_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_31_para.oz = 0;
  updateY_0_31_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_31_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_31_para.hx = 1;
  updateY_0_31_para.hy = 1;
  updateY_0_31_para.hz = 0;
  updateY_0_31_para.bx = 128;
  updateY_0_31_para.by = 8;
  updateY_0_31_para.bz = 1;
  updateY_0_31_para.mx = 1;
  updateY_0_31_para.my = 4;
  updateY_0_31_para.mz = 1;
  updateY_0_31_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.statejab7_13 = &state->jab7_13[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.stateh_12 = &state->h_12[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  updateY_0_31_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateY_0_31_para.oy][Proc.lon_hw -1+updateY_0_31_para.ox];
  
  athread_spawn(updateY_0_31, &updateY_0_31_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_13], async_jab7_13, state->jab7_13, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_32_info updateY_0_32_para;
  updateY_0_32_para.lz = 1 - 0;
  updateY_0_32_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_32_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_32_para.oz = 0;
  updateY_0_32_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_32_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_32_para.hx = 1;
  updateY_0_32_para.hy = 1;
  updateY_0_32_para.hz = 0;
  updateY_0_32_para.bx = 128;
  updateY_0_32_para.by = 8;
  updateY_0_32_para.bz = 1;
  updateY_0_32_para.mx = 1;
  updateY_0_32_para.my = 4;
  updateY_0_32_para.mz = 1;
  updateY_0_32_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.statejab7_13 = &state->jab7_13[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.stateu_11 = &state->u_11[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.stateu_12 = &state->u_12[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.stateu_13 = &state->u_13[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  updateY_0_32_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_32_para.oy][Proc.lon_hw -1+updateY_0_32_para.ox];
  
  athread_spawn(updateY_0_32, &updateY_0_32_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_13], async_jab7_13, state->jab7_13, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_13], async_jab_13, state->jab_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state->v_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state->v_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state->v_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_11], async_jab_11, state->jab_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_12], async_jab_12, state->jab_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_33_info updateY_0_33_para;
  updateY_0_33_para.lz = 1 - 0;
  updateY_0_33_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_33_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_33_para.oz = 0;
  updateY_0_33_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_33_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_33_para.hx = 1;
  updateY_0_33_para.hy = 1;
  updateY_0_33_para.hz = 0;
  updateY_0_33_para.bx = 128;
  updateY_0_33_para.by = 4;
  updateY_0_33_para.bz = 1;
  updateY_0_33_para.mx = 1;
  updateY_0_33_para.my = 4;
  updateY_0_33_para.mz = 1;
  updateY_0_33_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statejab7_13 = &state->jab7_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateh_13 = &state->h_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statejab_13 = &state->jab_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statev_11 = &state->v_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statev_12 = &state->v_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statev_13 = &state->v_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateh_11 = &state->h_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statejab_11 = &state->jab_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateu_11 = &state->u_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateh_12 = &state->h_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statejab_12 = &state->jab_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateu_12 = &state->u_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateu_13 = &state->u_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  updateY_0_33_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateY_0_33_para.oy][Proc.lon_hw -1+updateY_0_33_para.ox];
  
  athread_spawn(updateY_0_33, &updateY_0_33_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_2], async_qu_2, tend->qu_2, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_1], async_qu_1, tend->qu_1, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_11], async_vc_11, state->vc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_12], async_vc_12, state->vc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_13], async_vc_13, state->vc_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_11], async_uc_11, state->uc_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_12], async_uc_12, state->uc_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_13], async_uc_13, state->uc_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_34_info updateY_0_34_para;
  updateY_0_34_para.lz = 1 - 0;
  updateY_0_34_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_34_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_34_para.oz = 0;
  updateY_0_34_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_34_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_34_para.hx = 1;
  updateY_0_34_para.hy = 1;
  updateY_0_34_para.hz = 0;
  updateY_0_34_para.bx = 128;
  updateY_0_34_para.by = 4;
  updateY_0_34_para.bz = 1;
  updateY_0_34_para.mx = 1;
  updateY_0_34_para.my = 4;
  updateY_0_34_para.mz = 1;
  updateY_0_34_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.statevc_11 = &state->vc_11[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.statevc_12 = &state->vc_12[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.statevc_13 = &state->vc_13[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.stateuc_11 = &state->uc_11[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.stateuc_12 = &state->uc_12[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.stateuc_13 = &state->uc_13[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyh_11 = &tend->dyh_11[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyh_12 = &tend->dyh_12[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyh_13 = &tend->dyh_13[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyu_11 = &tend->dyu_11[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyu_12 = &tend->dyu_12[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyu_13 = &tend->dyu_13[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyv_11 = &tend->dyv_11[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyv_12 = &tend->dyv_12[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  updateY_0_34_para.tenddyv_13 = &tend->dyv_13[0][Proc.lat_hw -1+updateY_0_34_para.oy][Proc.lon_hw -1+updateY_0_34_para.ox];
  
  athread_spawn(updateY_0_34, &updateY_0_34_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYdyh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_23], async_jab7_23, state->jab7_23, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state->v_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_21], async_jab_21, state->jab_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state->v_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_35_info updateY_0_35_para;
  updateY_0_35_para.lz = 1 - 0;
  updateY_0_35_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_35_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_35_para.oz = 0;
  updateY_0_35_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_35_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_35_para.hx = 1;
  updateY_0_35_para.hy = 1;
  updateY_0_35_para.hz = 0;
  updateY_0_35_para.bx = 128;
  updateY_0_35_para.by = 4;
  updateY_0_35_para.bz = 1;
  updateY_0_35_para.mx = 1;
  updateY_0_35_para.my = 4;
  updateY_0_35_para.mz = 1;
  updateY_0_35_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statejab7_23 = &state->jab7_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statejab_23 = &state->jab_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateh_21 = &state->h_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateu_23 = &state->u_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateu_21 = &state->u_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statev_23 = &state->v_23[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statejab_21 = &state->jab_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.statev_21 = &state->v_21[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  updateY_0_35_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateY_0_35_para.oy][Proc.lon_hw -1+updateY_0_35_para.ox];
  
  athread_spawn(updateY_0_35, &updateY_0_35_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_23], async_jab7_23, state->jab7_23, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state->vc_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_36_info updateY_0_36_para;
  updateY_0_36_para.lz = 1 - 0;
  updateY_0_36_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_36_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_36_para.oz = 0;
  updateY_0_36_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_36_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_36_para.hx = 1;
  updateY_0_36_para.hy = 1;
  updateY_0_36_para.hz = 0;
  updateY_0_36_para.bx = 128;
  updateY_0_36_para.by = 8;
  updateY_0_36_para.bz = 1;
  updateY_0_36_para.mx = 1;
  updateY_0_36_para.my = 4;
  updateY_0_36_para.mz = 1;
  updateY_0_36_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.statejab7_23 = &state->jab7_23[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.statejab_23 = &state->jab_23[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.stateh_21 = &state->h_21[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.stateh_22 = &state->h_22[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.statevc_22 = &state->vc_22[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  updateY_0_36_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateY_0_36_para.oy][Proc.lon_hw -1+updateY_0_36_para.ox];
  
  athread_spawn(updateY_0_36, &updateY_0_36_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_23], async_jab7_23, state->jab7_23, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state->u_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_37_info updateY_0_37_para;
  updateY_0_37_para.lz = 1 - 0;
  updateY_0_37_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_37_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_37_para.oz = 0;
  updateY_0_37_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_37_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_37_para.hx = 1;
  updateY_0_37_para.hy = 1;
  updateY_0_37_para.hz = 0;
  updateY_0_37_para.bx = 128;
  updateY_0_37_para.by = 8;
  updateY_0_37_para.bz = 1;
  updateY_0_37_para.mx = 1;
  updateY_0_37_para.my = 4;
  updateY_0_37_para.mz = 1;
  updateY_0_37_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.statejab7_23 = &state->jab7_23[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.statejab_23 = &state->jab_23[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.stateu_21 = &state->u_21[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.stateu_22 = &state->u_22[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.stateu_23 = &state->u_23[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  updateY_0_37_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_37_para.oy][Proc.lon_hw -1+updateY_0_37_para.ox];
  
  athread_spawn(updateY_0_37, &updateY_0_37_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_23], async_jab7_23, state->jab7_23, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_23], async_jab_23, state->jab_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state->v_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state->v_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state->v_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_21], async_jab_21, state->jab_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_22], async_jab_22, state->jab_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state->uc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state->vc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_38_info updateY_0_38_para;
  updateY_0_38_para.lz = 1 - 0;
  updateY_0_38_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_38_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_38_para.oz = 0;
  updateY_0_38_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_38_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_38_para.hx = 1;
  updateY_0_38_para.hy = 1;
  updateY_0_38_para.hz = 0;
  updateY_0_38_para.bx = 128;
  updateY_0_38_para.by = 4;
  updateY_0_38_para.bz = 1;
  updateY_0_38_para.mx = 1;
  updateY_0_38_para.my = 4;
  updateY_0_38_para.mz = 1;
  updateY_0_38_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statejab7_23 = &state->jab7_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateh_23 = &state->h_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statejab_23 = &state->jab_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statev_21 = &state->v_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statev_22 = &state->v_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statev_23 = &state->v_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateh_21 = &state->h_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statejab_21 = &state->jab_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateu_21 = &state->u_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateh_22 = &state->h_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statejab_22 = &state->jab_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateu_22 = &state->u_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateuc_22 = &state->uc_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.statevc_22 = &state->vc_22[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateu_23 = &state->u_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  updateY_0_38_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateY_0_38_para.oy][Proc.lon_hw -1+updateY_0_38_para.ox];
  
  athread_spawn(updateY_0_38, &updateY_0_38_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_2], async_qu_2, tend->qu_2, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_1], async_qu_1, tend->qu_1, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_21], async_vc_21, state->vc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_22], async_vc_22, state->vc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_23], async_vc_23, state->vc_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_21], async_uc_21, state->uc_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_22], async_uc_22, state->uc_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_23], async_uc_23, state->uc_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_39_info updateY_0_39_para;
  updateY_0_39_para.lz = 1 - 0;
  updateY_0_39_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_39_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_39_para.oz = 0;
  updateY_0_39_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_39_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_39_para.hx = 1;
  updateY_0_39_para.hy = 1;
  updateY_0_39_para.hz = 0;
  updateY_0_39_para.bx = 128;
  updateY_0_39_para.by = 4;
  updateY_0_39_para.bz = 1;
  updateY_0_39_para.mx = 1;
  updateY_0_39_para.my = 4;
  updateY_0_39_para.mz = 1;
  updateY_0_39_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.statevc_21 = &state->vc_21[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.statevc_22 = &state->vc_22[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.statevc_23 = &state->vc_23[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.stateuc_21 = &state->uc_21[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.stateuc_22 = &state->uc_22[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.stateuc_23 = &state->uc_23[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyh_11 = &tend->dyh_11[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyh_12 = &tend->dyh_12[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyh_13 = &tend->dyh_13[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyu_11 = &tend->dyu_11[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyu_12 = &tend->dyu_12[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyu_13 = &tend->dyu_13[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyv_11 = &tend->dyv_11[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyv_12 = &tend->dyv_12[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  updateY_0_39_para.tenddyv_13 = &tend->dyv_13[0][Proc.lat_hw -1+updateY_0_39_para.oy][Proc.lon_hw -1+updateY_0_39_para.ox];
  
  athread_spawn(updateY_0_39, &updateY_0_39_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYdyh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_33], async_jab7_33, state->jab7_33, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_40_info updateY_0_40_para;
  updateY_0_40_para.lz = 1 - 0;
  updateY_0_40_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_40_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_40_para.oz = 0;
  updateY_0_40_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_40_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_40_para.hx = 1;
  updateY_0_40_para.hy = 1;
  updateY_0_40_para.hz = 0;
  updateY_0_40_para.bx = 128;
  updateY_0_40_para.by = 4;
  updateY_0_40_para.bz = 1;
  updateY_0_40_para.mx = 1;
  updateY_0_40_para.my = 4;
  updateY_0_40_para.mz = 1;
  updateY_0_40_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statejab7_33 = &state->jab7_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateu_33 = &state->u_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateu_31 = &state->u_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statev_33 = &state->v_33[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.statev_31 = &state->v_31[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendqh_1 = &tend->qh_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  updateY_0_40_para.tendqv_1 = &tend->qv_1[0][Proc.lat_hw -1+updateY_0_40_para.oy][Proc.lon_hw -1+updateY_0_40_para.ox];
  
  athread_spawn(updateY_0_40, &updateY_0_40_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_1Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fh_1[0][0][0], &Proc.FieldReq[async_fh_1], async_fh_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->fu_1[0][0][0], &Proc.FieldReq[async_fu_1], async_fu_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->fv_1[0][0][0], &Proc.FieldReq[async_fv_1], async_fv_1, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qu_1[0][0][0], &Proc.FieldReq[async_qu_1], async_qu_1, true, true, false, false, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qv_1[0][0][0], &Proc.FieldReq[async_qv_1], async_qv_1, true, true, false, true, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_33], async_jab7_33, state->jab7_33, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_41_info updateY_0_41_para;
  updateY_0_41_para.lz = 1 - 0;
  updateY_0_41_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_41_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_41_para.oz = 0;
  updateY_0_41_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_41_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_41_para.hx = 1;
  updateY_0_41_para.hy = 1;
  updateY_0_41_para.hz = 0;
  updateY_0_41_para.bx = 128;
  updateY_0_41_para.by = 8;
  updateY_0_41_para.bz = 1;
  updateY_0_41_para.mx = 1;
  updateY_0_41_para.my = 4;
  updateY_0_41_para.mz = 1;
  updateY_0_41_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.statejab7_33 = &state->jab7_33[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  updateY_0_41_para.tendqh_2 = &tend->qh_2[0][Proc.lat_hw -1+updateY_0_41_para.oy][Proc.lon_hw -1+updateY_0_41_para.ox];
  
  athread_spawn(updateY_0_41, &updateY_0_41_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfh_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fh_2[0][0][0], &Proc.FieldReq[async_fh_2], async_fh_2, true, true, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_33], async_jab7_33, state->jab7_33, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_42_info updateY_0_42_para;
  updateY_0_42_para.lz = 1 - 0;
  updateY_0_42_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_42_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_42_para.oz = 0;
  updateY_0_42_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_42_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_42_para.hx = 1;
  updateY_0_42_para.hy = 1;
  updateY_0_42_para.hz = 0;
  updateY_0_42_para.bx = 128;
  updateY_0_42_para.by = 8;
  updateY_0_42_para.bz = 1;
  updateY_0_42_para.mx = 1;
  updateY_0_42_para.my = 4;
  updateY_0_42_para.mz = 1;
  updateY_0_42_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.statejab7_33 = &state->jab7_33[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.stateu_31 = &state->u_31[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.stateu_32 = &state->u_32[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.stateu_33 = &state->u_33[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  updateY_0_42_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_42_para.oy][Proc.lon_hw -1+updateY_0_42_para.ox];
  
  athread_spawn(updateY_0_42, &updateY_0_42_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfu_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fu_2[0][0][0], &Proc.FieldReq[async_fu_2], async_fu_2, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qu_2[0][0][0], &Proc.FieldReq[async_qu_2], async_qu_2, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab7_33], async_jab7_33, state->jab7_33, true, true, true, false, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_33], async_jab_33, state->jab_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state->v_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state->v_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state->v_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_31], async_jab_31, state->jab_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_jab_32], async_jab_32, state->jab_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_43_info updateY_0_43_para;
  updateY_0_43_para.lz = 1 - 0;
  updateY_0_43_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_43_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_43_para.oz = 0;
  updateY_0_43_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_43_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_43_para.hx = 1;
  updateY_0_43_para.hy = 1;
  updateY_0_43_para.hz = 0;
  updateY_0_43_para.bx = 128;
  updateY_0_43_para.by = 4;
  updateY_0_43_para.bz = 1;
  updateY_0_43_para.mx = 1;
  updateY_0_43_para.my = 4;
  updateY_0_43_para.mz = 1;
  updateY_0_43_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statejab7_33 = &state->jab7_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateh_33 = &state->h_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statejab_33 = &state->jab_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statev_31 = &state->v_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statev_32 = &state->v_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statev_33 = &state->v_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateh_31 = &state->h_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statejab_31 = &state->jab_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateu_31 = &state->u_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateh_32 = &state->h_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statejab_32 = &state->jab_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateu_32 = &state->u_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateu_33 = &state->u_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  updateY_0_43_para.tendqv_2 = &tend->qv_2[0][Proc.lat_hw -1+updateY_0_43_para.oy][Proc.lon_hw -1+updateY_0_43_para.ox];
  
  athread_spawn(updateY_0_43, &updateY_0_43_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYfv_2Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &tend->fv_2[0][0][0], &Proc.FieldReq[async_fv_2], async_fv_2, true, true, false, true, false, true);
  UpdateHaloCS_2d_D(Proc, &tend->qv_2[0][0][0], &Proc.FieldReq[async_qv_2], async_qv_2, true, true, false, true, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_2], async_qu_2, tend->qu_2, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_qu_1], async_qu_1, tend->qu_1, true, true, true, false, false, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_2], async_fh_2, tend->fh_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fh_1], async_fh_1, tend->fh_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_2], async_fu_2, tend->fu_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_31], async_vc_31, state->vc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fu_1], async_fu_1, tend->fu_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_32], async_vc_32, state->vc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_vc_33], async_vc_33, state->vc_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_2], async_fv_2, tend->fv_2, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_31], async_uc_31, state->uc_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_fv_1], async_fv_1, tend->fv_1, true, true, true, false, true, false, true);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_32], async_uc_32, state->uc_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_uc_33], async_uc_33, state->uc_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  updateY_0_44_info updateY_0_44_para;
  updateY_0_44_para.lz = 1 - 0;
  updateY_0_44_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  updateY_0_44_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  updateY_0_44_para.oz = 0;
  updateY_0_44_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  updateY_0_44_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  updateY_0_44_para.hx = 1;
  updateY_0_44_para.hy = 1;
  updateY_0_44_para.hz = 0;
  updateY_0_44_para.bx = 128;
  updateY_0_44_para.by = 4;
  updateY_0_44_para.bz = 1;
  updateY_0_44_para.mx = 1;
  updateY_0_44_para.my = 4;
  updateY_0_44_para.mz = 1;
  updateY_0_44_para.tendqu_2 = &tend->qu_2[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendqu_1 = &tend->qu_1[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfh_2 = &tend->fh_2[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfh_1 = &tend->fh_1[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfu_2 = &tend->fu_2[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.statevc_31 = &state->vc_31[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfu_1 = &tend->fu_1[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.statevc_32 = &state->vc_32[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.statevc_33 = &state->vc_33[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfv_2 = &tend->fv_2[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.stateuc_31 = &state->uc_31[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tendfv_1 = &tend->fv_1[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.stateuc_32 = &state->uc_32[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.stateuc_33 = &state->uc_33[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyh_11 = &tend->dyh_11[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyh_12 = &tend->dyh_12[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyh_13 = &tend->dyh_13[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyu_11 = &tend->dyu_11[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyu_12 = &tend->dyu_12[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyu_13 = &tend->dyu_13[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyv_11 = &tend->dyv_11[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyv_12 = &tend->dyv_12[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  updateY_0_44_para.tenddyv_13 = &tend->dyv_13[0][Proc.lat_hw -1+updateY_0_44_para.oy][Proc.lon_hw -1+updateY_0_44_para.ox];
  
  athread_spawn(updateY_0_44, &updateY_0_44_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  updateYdyh_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void update_rk1_0(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridTendField* tend1, double dt)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state_old->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state_old->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state_old->v_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_0_info update_rk1_0_0_para;
  update_rk1_0_0_para.lz = 1 - 0;
  update_rk1_0_0_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_0_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_0_para.oz = 0;
  update_rk1_0_0_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_0_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_0_para.hx = 1;
  update_rk1_0_0_para.hy = 1;
  update_rk1_0_0_para.hz = 0;
  update_rk1_0_0_para.bx = 128;
  update_rk1_0_0_para.by = 8;
  update_rk1_0_0_para.bz = 1;
  update_rk1_0_0_para.mx = 1;
  update_rk1_0_0_para.my = 4;
  update_rk1_0_0_para.mz = 1;
  update_rk1_0_0_para.state_oldh_11 = &state_old->h_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dxh_11 = &tend1->dxh_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dyh_11 = &tend1->dyh_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.state_oldu_11 = &state_old->u_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dxu_11 = &tend1->dxu_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dyu_11 = &tend1->dyu_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.state_oldv_11 = &state_old->v_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dxv_11 = &tend1->dxv_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.tend1dyv_11 = &tend1->dyv_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.state_newh_11 = &state_new->h_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.state_newu_11 = &state_new->u_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.state_newv_11 = &state_new->v_11[0][Proc.lat_hw -1+update_rk1_0_0_para.oy][Proc.lon_hw -1+update_rk1_0_0_para.ox];
  update_rk1_0_0_para.dt = dt;
  
  athread_spawn(update_rk1_0_0, &update_rk1_0_0_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_11[0][0][0], &Proc.FieldReq[async_h_11], async_h_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_11[0][0][0], &Proc.FieldReq[async_u_11], async_u_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_11[0][0][0], &Proc.FieldReq[async_v_11], async_v_11, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state_old->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state_old->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state_old->v_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_1_info update_rk1_0_1_para;
  update_rk1_0_1_para.lz = 1 - 0;
  update_rk1_0_1_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_1_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_1_para.oz = 0;
  update_rk1_0_1_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_1_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_1_para.hx = 1;
  update_rk1_0_1_para.hy = 1;
  update_rk1_0_1_para.hz = 0;
  update_rk1_0_1_para.bx = 128;
  update_rk1_0_1_para.by = 8;
  update_rk1_0_1_para.bz = 1;
  update_rk1_0_1_para.mx = 1;
  update_rk1_0_1_para.my = 4;
  update_rk1_0_1_para.mz = 1;
  update_rk1_0_1_para.state_oldh_12 = &state_old->h_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dxh_12 = &tend1->dxh_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dyh_12 = &tend1->dyh_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.state_oldu_12 = &state_old->u_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dxu_12 = &tend1->dxu_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dyu_12 = &tend1->dyu_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.state_oldv_12 = &state_old->v_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dxv_12 = &tend1->dxv_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.tend1dyv_12 = &tend1->dyv_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.state_newh_12 = &state_new->h_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.state_newu_12 = &state_new->u_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.state_newv_12 = &state_new->v_12[0][Proc.lat_hw -1+update_rk1_0_1_para.oy][Proc.lon_hw -1+update_rk1_0_1_para.ox];
  update_rk1_0_1_para.dt = dt;
  
  athread_spawn(update_rk1_0_1, &update_rk1_0_1_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_12[0][0][0], &Proc.FieldReq[async_h_12], async_h_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_12[0][0][0], &Proc.FieldReq[async_u_12], async_u_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_12[0][0][0], &Proc.FieldReq[async_v_12], async_v_12, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state_old->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state_old->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state_old->v_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_2_info update_rk1_0_2_para;
  update_rk1_0_2_para.lz = 1 - 0;
  update_rk1_0_2_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_2_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_2_para.oz = 0;
  update_rk1_0_2_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_2_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_2_para.hx = 1;
  update_rk1_0_2_para.hy = 1;
  update_rk1_0_2_para.hz = 0;
  update_rk1_0_2_para.bx = 128;
  update_rk1_0_2_para.by = 8;
  update_rk1_0_2_para.bz = 1;
  update_rk1_0_2_para.mx = 1;
  update_rk1_0_2_para.my = 4;
  update_rk1_0_2_para.mz = 1;
  update_rk1_0_2_para.state_oldh_13 = &state_old->h_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dxh_13 = &tend1->dxh_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dyh_13 = &tend1->dyh_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.state_oldu_13 = &state_old->u_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dxu_13 = &tend1->dxu_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dyu_13 = &tend1->dyu_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.state_oldv_13 = &state_old->v_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dxv_13 = &tend1->dxv_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.tend1dyv_13 = &tend1->dyv_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.state_newh_13 = &state_new->h_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.state_newu_13 = &state_new->u_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.state_newv_13 = &state_new->v_13[0][Proc.lat_hw -1+update_rk1_0_2_para.oy][Proc.lon_hw -1+update_rk1_0_2_para.ox];
  update_rk1_0_2_para.dt = dt;
  
  athread_spawn(update_rk1_0_2, &update_rk1_0_2_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_13[0][0][0], &Proc.FieldReq[async_h_13], async_h_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_13[0][0][0], &Proc.FieldReq[async_u_13], async_u_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_13[0][0][0], &Proc.FieldReq[async_v_13], async_v_13, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state_old->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state_old->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state_old->v_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_3_info update_rk1_0_3_para;
  update_rk1_0_3_para.lz = 1 - 0;
  update_rk1_0_3_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_3_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_3_para.oz = 0;
  update_rk1_0_3_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_3_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_3_para.hx = 1;
  update_rk1_0_3_para.hy = 1;
  update_rk1_0_3_para.hz = 0;
  update_rk1_0_3_para.bx = 128;
  update_rk1_0_3_para.by = 8;
  update_rk1_0_3_para.bz = 1;
  update_rk1_0_3_para.mx = 1;
  update_rk1_0_3_para.my = 4;
  update_rk1_0_3_para.mz = 1;
  update_rk1_0_3_para.state_oldh_21 = &state_old->h_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dxh_21 = &tend1->dxh_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dyh_21 = &tend1->dyh_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.state_oldu_21 = &state_old->u_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dxu_21 = &tend1->dxu_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dyu_21 = &tend1->dyu_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.state_oldv_21 = &state_old->v_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dxv_21 = &tend1->dxv_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.tend1dyv_21 = &tend1->dyv_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.state_newh_21 = &state_new->h_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.state_newu_21 = &state_new->u_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.state_newv_21 = &state_new->v_21[0][Proc.lat_hw -1+update_rk1_0_3_para.oy][Proc.lon_hw -1+update_rk1_0_3_para.ox];
  update_rk1_0_3_para.dt = dt;
  
  athread_spawn(update_rk1_0_3, &update_rk1_0_3_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_21Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_21[0][0][0], &Proc.FieldReq[async_h_21], async_h_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_21[0][0][0], &Proc.FieldReq[async_u_21], async_u_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_21[0][0][0], &Proc.FieldReq[async_v_21], async_v_21, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state_old->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state_old->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state_old->v_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_4_info update_rk1_0_4_para;
  update_rk1_0_4_para.lz = 1 - 0;
  update_rk1_0_4_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_4_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_4_para.oz = 0;
  update_rk1_0_4_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_4_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_4_para.hx = 1;
  update_rk1_0_4_para.hy = 1;
  update_rk1_0_4_para.hz = 0;
  update_rk1_0_4_para.bx = 128;
  update_rk1_0_4_para.by = 8;
  update_rk1_0_4_para.bz = 1;
  update_rk1_0_4_para.mx = 1;
  update_rk1_0_4_para.my = 4;
  update_rk1_0_4_para.mz = 1;
  update_rk1_0_4_para.state_oldh_22 = &state_old->h_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dxh_22 = &tend1->dxh_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dyh_22 = &tend1->dyh_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.state_oldu_22 = &state_old->u_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dxu_22 = &tend1->dxu_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dyu_22 = &tend1->dyu_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.state_oldv_22 = &state_old->v_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dxv_22 = &tend1->dxv_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.tend1dyv_22 = &tend1->dyv_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.state_newh_22 = &state_new->h_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.state_newu_22 = &state_new->u_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.state_newv_22 = &state_new->v_22[0][Proc.lat_hw -1+update_rk1_0_4_para.oy][Proc.lon_hw -1+update_rk1_0_4_para.ox];
  update_rk1_0_4_para.dt = dt;
  
  athread_spawn(update_rk1_0_4, &update_rk1_0_4_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_22Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_22[0][0][0], &Proc.FieldReq[async_h_22], async_h_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_22[0][0][0], &Proc.FieldReq[async_u_22], async_u_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_22[0][0][0], &Proc.FieldReq[async_v_22], async_v_22, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state_old->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state_old->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state_old->v_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_5_info update_rk1_0_5_para;
  update_rk1_0_5_para.lz = 1 - 0;
  update_rk1_0_5_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_5_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_5_para.oz = 0;
  update_rk1_0_5_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_5_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_5_para.hx = 1;
  update_rk1_0_5_para.hy = 1;
  update_rk1_0_5_para.hz = 0;
  update_rk1_0_5_para.bx = 128;
  update_rk1_0_5_para.by = 8;
  update_rk1_0_5_para.bz = 1;
  update_rk1_0_5_para.mx = 1;
  update_rk1_0_5_para.my = 4;
  update_rk1_0_5_para.mz = 1;
  update_rk1_0_5_para.state_oldh_23 = &state_old->h_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dxh_23 = &tend1->dxh_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dyh_23 = &tend1->dyh_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.state_oldu_23 = &state_old->u_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dxu_23 = &tend1->dxu_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dyu_23 = &tend1->dyu_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.state_oldv_23 = &state_old->v_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dxv_23 = &tend1->dxv_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.tend1dyv_23 = &tend1->dyv_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.state_newh_23 = &state_new->h_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.state_newu_23 = &state_new->u_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.state_newv_23 = &state_new->v_23[0][Proc.lat_hw -1+update_rk1_0_5_para.oy][Proc.lon_hw -1+update_rk1_0_5_para.ox];
  update_rk1_0_5_para.dt = dt;
  
  athread_spawn(update_rk1_0_5, &update_rk1_0_5_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_23Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_23[0][0][0], &Proc.FieldReq[async_h_23], async_h_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_23[0][0][0], &Proc.FieldReq[async_u_23], async_u_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_23[0][0][0], &Proc.FieldReq[async_v_23], async_v_23, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state_old->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state_old->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state_old->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_6_info update_rk1_0_6_para;
  update_rk1_0_6_para.lz = 1 - 0;
  update_rk1_0_6_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_6_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_6_para.oz = 0;
  update_rk1_0_6_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_6_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_6_para.hx = 1;
  update_rk1_0_6_para.hy = 1;
  update_rk1_0_6_para.hz = 0;
  update_rk1_0_6_para.bx = 128;
  update_rk1_0_6_para.by = 8;
  update_rk1_0_6_para.bz = 1;
  update_rk1_0_6_para.mx = 1;
  update_rk1_0_6_para.my = 4;
  update_rk1_0_6_para.mz = 1;
  update_rk1_0_6_para.state_oldh_31 = &state_old->h_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dxh_31 = &tend1->dxh_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dyh_31 = &tend1->dyh_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.state_oldu_31 = &state_old->u_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dxu_31 = &tend1->dxu_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dyu_31 = &tend1->dyu_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.state_oldv_31 = &state_old->v_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dxv_31 = &tend1->dxv_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.tend1dyv_31 = &tend1->dyv_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.state_newh_31 = &state_new->h_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.state_newu_31 = &state_new->u_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.state_newv_31 = &state_new->v_31[0][Proc.lat_hw -1+update_rk1_0_6_para.oy][Proc.lon_hw -1+update_rk1_0_6_para.ox];
  update_rk1_0_6_para.dt = dt;
  
  athread_spawn(update_rk1_0_6, &update_rk1_0_6_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_31Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_31[0][0][0], &Proc.FieldReq[async_h_31], async_h_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_31[0][0][0], &Proc.FieldReq[async_u_31], async_u_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_31[0][0][0], &Proc.FieldReq[async_v_31], async_v_31, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state_old->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state_old->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state_old->v_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_7_info update_rk1_0_7_para;
  update_rk1_0_7_para.lz = 1 - 0;
  update_rk1_0_7_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_7_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_7_para.oz = 0;
  update_rk1_0_7_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_7_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_7_para.hx = 1;
  update_rk1_0_7_para.hy = 1;
  update_rk1_0_7_para.hz = 0;
  update_rk1_0_7_para.bx = 128;
  update_rk1_0_7_para.by = 8;
  update_rk1_0_7_para.bz = 1;
  update_rk1_0_7_para.mx = 1;
  update_rk1_0_7_para.my = 4;
  update_rk1_0_7_para.mz = 1;
  update_rk1_0_7_para.state_oldh_32 = &state_old->h_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dxh_32 = &tend1->dxh_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dyh_32 = &tend1->dyh_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.state_oldu_32 = &state_old->u_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dxu_32 = &tend1->dxu_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dyu_32 = &tend1->dyu_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.state_oldv_32 = &state_old->v_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dxv_32 = &tend1->dxv_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.tend1dyv_32 = &tend1->dyv_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.state_newh_32 = &state_new->h_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.state_newu_32 = &state_new->u_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.state_newv_32 = &state_new->v_32[0][Proc.lat_hw -1+update_rk1_0_7_para.oy][Proc.lon_hw -1+update_rk1_0_7_para.ox];
  update_rk1_0_7_para.dt = dt;
  
  athread_spawn(update_rk1_0_7, &update_rk1_0_7_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_32Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_32[0][0][0], &Proc.FieldReq[async_h_32], async_h_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_32[0][0][0], &Proc.FieldReq[async_u_32], async_u_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_32[0][0][0], &Proc.FieldReq[async_v_32], async_v_32, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state_old->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state_old->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state_old->v_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk1_0_8_info update_rk1_0_8_para;
  update_rk1_0_8_para.lz = 1 - 0;
  update_rk1_0_8_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk1_0_8_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk1_0_8_para.oz = 0;
  update_rk1_0_8_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk1_0_8_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk1_0_8_para.hx = 1;
  update_rk1_0_8_para.hy = 1;
  update_rk1_0_8_para.hz = 0;
  update_rk1_0_8_para.bx = 128;
  update_rk1_0_8_para.by = 8;
  update_rk1_0_8_para.bz = 1;
  update_rk1_0_8_para.mx = 1;
  update_rk1_0_8_para.my = 4;
  update_rk1_0_8_para.mz = 1;
  update_rk1_0_8_para.state_oldh_33 = &state_old->h_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dxh_33 = &tend1->dxh_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dyh_33 = &tend1->dyh_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.state_oldu_33 = &state_old->u_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dxu_33 = &tend1->dxu_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dyu_33 = &tend1->dyu_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.state_oldv_33 = &state_old->v_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dxv_33 = &tend1->dxv_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.tend1dyv_33 = &tend1->dyv_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.state_newh_33 = &state_new->h_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.state_newu_33 = &state_new->u_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.state_newv_33 = &state_new->v_33[0][Proc.lat_hw -1+update_rk1_0_8_para.oy][Proc.lon_hw -1+update_rk1_0_8_para.ox];
  update_rk1_0_8_para.dt = dt;
  
  athread_spawn(update_rk1_0_8, &update_rk1_0_8_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk1h_33Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_33[0][0][0], &Proc.FieldReq[async_h_33], async_h_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_33[0][0][0], &Proc.FieldReq[async_u_33], async_u_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_33[0][0][0], &Proc.FieldReq[async_v_33], async_v_33, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void update_rk2_0(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridTendField* tend1, struct HybridTendField* tend2, double dt)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state_old->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state_old->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state_old->v_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_0_info update_rk2_0_0_para;
  update_rk2_0_0_para.lz = 1 - 0;
  update_rk2_0_0_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_0_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_0_para.oz = 0;
  update_rk2_0_0_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_0_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_0_para.hx = 1;
  update_rk2_0_0_para.hy = 1;
  update_rk2_0_0_para.hz = 0;
  update_rk2_0_0_para.bx = 128;
  update_rk2_0_0_para.by = 4;
  update_rk2_0_0_para.bz = 1;
  update_rk2_0_0_para.mx = 1;
  update_rk2_0_0_para.my = 4;
  update_rk2_0_0_para.mz = 1;
  update_rk2_0_0_para.state_oldh_11 = &state_old->h_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dxh_11 = &tend1->dxh_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dyh_11 = &tend1->dyh_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.state_oldu_11 = &state_old->u_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dxu_11 = &tend1->dxu_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dyu_11 = &tend1->dyu_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.state_oldv_11 = &state_old->v_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dxv_11 = &tend1->dxv_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend1dyv_11 = &tend1->dyv_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dxh_11 = &tend2->dxh_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dyh_11 = &tend2->dyh_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dxu_11 = &tend2->dxu_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dyu_11 = &tend2->dyu_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dxv_11 = &tend2->dxv_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.tend2dyv_11 = &tend2->dyv_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.state_newh_11 = &state_new->h_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.state_newu_11 = &state_new->u_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.state_newv_11 = &state_new->v_11[0][Proc.lat_hw -1+update_rk2_0_0_para.oy][Proc.lon_hw -1+update_rk2_0_0_para.ox];
  update_rk2_0_0_para.dt = dt;
  
  athread_spawn(update_rk2_0_0, &update_rk2_0_0_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_11[0][0][0], &Proc.FieldReq[async_h_11], async_h_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_11[0][0][0], &Proc.FieldReq[async_u_11], async_u_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_11[0][0][0], &Proc.FieldReq[async_v_11], async_v_11, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state_old->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state_old->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state_old->v_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_1_info update_rk2_0_1_para;
  update_rk2_0_1_para.lz = 1 - 0;
  update_rk2_0_1_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_1_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_1_para.oz = 0;
  update_rk2_0_1_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_1_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_1_para.hx = 1;
  update_rk2_0_1_para.hy = 1;
  update_rk2_0_1_para.hz = 0;
  update_rk2_0_1_para.bx = 128;
  update_rk2_0_1_para.by = 4;
  update_rk2_0_1_para.bz = 1;
  update_rk2_0_1_para.mx = 1;
  update_rk2_0_1_para.my = 4;
  update_rk2_0_1_para.mz = 1;
  update_rk2_0_1_para.state_oldh_12 = &state_old->h_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dxh_12 = &tend1->dxh_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dyh_12 = &tend1->dyh_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.state_oldu_12 = &state_old->u_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dxu_12 = &tend1->dxu_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dyu_12 = &tend1->dyu_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.state_oldv_12 = &state_old->v_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dxv_12 = &tend1->dxv_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend1dyv_12 = &tend1->dyv_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dxh_12 = &tend2->dxh_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dyh_12 = &tend2->dyh_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dxu_12 = &tend2->dxu_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dyu_12 = &tend2->dyu_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dxv_12 = &tend2->dxv_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.tend2dyv_12 = &tend2->dyv_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.state_newh_12 = &state_new->h_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.state_newu_12 = &state_new->u_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.state_newv_12 = &state_new->v_12[0][Proc.lat_hw -1+update_rk2_0_1_para.oy][Proc.lon_hw -1+update_rk2_0_1_para.ox];
  update_rk2_0_1_para.dt = dt;
  
  athread_spawn(update_rk2_0_1, &update_rk2_0_1_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_12[0][0][0], &Proc.FieldReq[async_h_12], async_h_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_12[0][0][0], &Proc.FieldReq[async_u_12], async_u_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_12[0][0][0], &Proc.FieldReq[async_v_12], async_v_12, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state_old->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state_old->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state_old->v_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_2_info update_rk2_0_2_para;
  update_rk2_0_2_para.lz = 1 - 0;
  update_rk2_0_2_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_2_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_2_para.oz = 0;
  update_rk2_0_2_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_2_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_2_para.hx = 1;
  update_rk2_0_2_para.hy = 1;
  update_rk2_0_2_para.hz = 0;
  update_rk2_0_2_para.bx = 128;
  update_rk2_0_2_para.by = 4;
  update_rk2_0_2_para.bz = 1;
  update_rk2_0_2_para.mx = 1;
  update_rk2_0_2_para.my = 4;
  update_rk2_0_2_para.mz = 1;
  update_rk2_0_2_para.state_oldh_13 = &state_old->h_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dxh_13 = &tend1->dxh_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dyh_13 = &tend1->dyh_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.state_oldu_13 = &state_old->u_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dxu_13 = &tend1->dxu_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dyu_13 = &tend1->dyu_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.state_oldv_13 = &state_old->v_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dxv_13 = &tend1->dxv_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend1dyv_13 = &tend1->dyv_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dxh_13 = &tend2->dxh_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dyh_13 = &tend2->dyh_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dxu_13 = &tend2->dxu_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dyu_13 = &tend2->dyu_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dxv_13 = &tend2->dxv_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.tend2dyv_13 = &tend2->dyv_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.state_newh_13 = &state_new->h_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.state_newu_13 = &state_new->u_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.state_newv_13 = &state_new->v_13[0][Proc.lat_hw -1+update_rk2_0_2_para.oy][Proc.lon_hw -1+update_rk2_0_2_para.ox];
  update_rk2_0_2_para.dt = dt;
  
  athread_spawn(update_rk2_0_2, &update_rk2_0_2_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_13[0][0][0], &Proc.FieldReq[async_h_13], async_h_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_13[0][0][0], &Proc.FieldReq[async_u_13], async_u_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_13[0][0][0], &Proc.FieldReq[async_v_13], async_v_13, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state_old->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state_old->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state_old->v_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_3_info update_rk2_0_3_para;
  update_rk2_0_3_para.lz = 1 - 0;
  update_rk2_0_3_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_3_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_3_para.oz = 0;
  update_rk2_0_3_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_3_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_3_para.hx = 1;
  update_rk2_0_3_para.hy = 1;
  update_rk2_0_3_para.hz = 0;
  update_rk2_0_3_para.bx = 128;
  update_rk2_0_3_para.by = 4;
  update_rk2_0_3_para.bz = 1;
  update_rk2_0_3_para.mx = 1;
  update_rk2_0_3_para.my = 4;
  update_rk2_0_3_para.mz = 1;
  update_rk2_0_3_para.state_oldh_21 = &state_old->h_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dxh_21 = &tend1->dxh_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dyh_21 = &tend1->dyh_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.state_oldu_21 = &state_old->u_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dxu_21 = &tend1->dxu_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dyu_21 = &tend1->dyu_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.state_oldv_21 = &state_old->v_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dxv_21 = &tend1->dxv_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend1dyv_21 = &tend1->dyv_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dxh_21 = &tend2->dxh_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dyh_21 = &tend2->dyh_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dxu_21 = &tend2->dxu_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dyu_21 = &tend2->dyu_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dxv_21 = &tend2->dxv_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.tend2dyv_21 = &tend2->dyv_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.state_newh_21 = &state_new->h_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.state_newu_21 = &state_new->u_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.state_newv_21 = &state_new->v_21[0][Proc.lat_hw -1+update_rk2_0_3_para.oy][Proc.lon_hw -1+update_rk2_0_3_para.ox];
  update_rk2_0_3_para.dt = dt;
  
  athread_spawn(update_rk2_0_3, &update_rk2_0_3_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_21Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_21[0][0][0], &Proc.FieldReq[async_h_21], async_h_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_21[0][0][0], &Proc.FieldReq[async_u_21], async_u_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_21[0][0][0], &Proc.FieldReq[async_v_21], async_v_21, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state_old->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state_old->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state_old->v_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_4_info update_rk2_0_4_para;
  update_rk2_0_4_para.lz = 1 - 0;
  update_rk2_0_4_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_4_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_4_para.oz = 0;
  update_rk2_0_4_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_4_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_4_para.hx = 1;
  update_rk2_0_4_para.hy = 1;
  update_rk2_0_4_para.hz = 0;
  update_rk2_0_4_para.bx = 128;
  update_rk2_0_4_para.by = 4;
  update_rk2_0_4_para.bz = 1;
  update_rk2_0_4_para.mx = 1;
  update_rk2_0_4_para.my = 4;
  update_rk2_0_4_para.mz = 1;
  update_rk2_0_4_para.state_oldh_22 = &state_old->h_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dxh_22 = &tend1->dxh_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dyh_22 = &tend1->dyh_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.state_oldu_22 = &state_old->u_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dxu_22 = &tend1->dxu_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dyu_22 = &tend1->dyu_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.state_oldv_22 = &state_old->v_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dxv_22 = &tend1->dxv_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend1dyv_22 = &tend1->dyv_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dxh_22 = &tend2->dxh_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dyh_22 = &tend2->dyh_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dxu_22 = &tend2->dxu_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dyu_22 = &tend2->dyu_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dxv_22 = &tend2->dxv_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.tend2dyv_22 = &tend2->dyv_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.state_newh_22 = &state_new->h_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.state_newu_22 = &state_new->u_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.state_newv_22 = &state_new->v_22[0][Proc.lat_hw -1+update_rk2_0_4_para.oy][Proc.lon_hw -1+update_rk2_0_4_para.ox];
  update_rk2_0_4_para.dt = dt;
  
  athread_spawn(update_rk2_0_4, &update_rk2_0_4_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_22Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_22[0][0][0], &Proc.FieldReq[async_h_22], async_h_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_22[0][0][0], &Proc.FieldReq[async_u_22], async_u_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_22[0][0][0], &Proc.FieldReq[async_v_22], async_v_22, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state_old->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state_old->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state_old->v_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_5_info update_rk2_0_5_para;
  update_rk2_0_5_para.lz = 1 - 0;
  update_rk2_0_5_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_5_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_5_para.oz = 0;
  update_rk2_0_5_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_5_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_5_para.hx = 1;
  update_rk2_0_5_para.hy = 1;
  update_rk2_0_5_para.hz = 0;
  update_rk2_0_5_para.bx = 128;
  update_rk2_0_5_para.by = 4;
  update_rk2_0_5_para.bz = 1;
  update_rk2_0_5_para.mx = 1;
  update_rk2_0_5_para.my = 4;
  update_rk2_0_5_para.mz = 1;
  update_rk2_0_5_para.state_oldh_23 = &state_old->h_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dxh_23 = &tend1->dxh_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dyh_23 = &tend1->dyh_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.state_oldu_23 = &state_old->u_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dxu_23 = &tend1->dxu_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dyu_23 = &tend1->dyu_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.state_oldv_23 = &state_old->v_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dxv_23 = &tend1->dxv_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend1dyv_23 = &tend1->dyv_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dxh_23 = &tend2->dxh_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dyh_23 = &tend2->dyh_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dxu_23 = &tend2->dxu_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dyu_23 = &tend2->dyu_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dxv_23 = &tend2->dxv_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.tend2dyv_23 = &tend2->dyv_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.state_newh_23 = &state_new->h_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.state_newu_23 = &state_new->u_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.state_newv_23 = &state_new->v_23[0][Proc.lat_hw -1+update_rk2_0_5_para.oy][Proc.lon_hw -1+update_rk2_0_5_para.ox];
  update_rk2_0_5_para.dt = dt;
  
  athread_spawn(update_rk2_0_5, &update_rk2_0_5_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_23Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_23[0][0][0], &Proc.FieldReq[async_h_23], async_h_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_23[0][0][0], &Proc.FieldReq[async_u_23], async_u_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_23[0][0][0], &Proc.FieldReq[async_v_23], async_v_23, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state_old->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state_old->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state_old->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_6_info update_rk2_0_6_para;
  update_rk2_0_6_para.lz = 1 - 0;
  update_rk2_0_6_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_6_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_6_para.oz = 0;
  update_rk2_0_6_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_6_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_6_para.hx = 1;
  update_rk2_0_6_para.hy = 1;
  update_rk2_0_6_para.hz = 0;
  update_rk2_0_6_para.bx = 128;
  update_rk2_0_6_para.by = 4;
  update_rk2_0_6_para.bz = 1;
  update_rk2_0_6_para.mx = 1;
  update_rk2_0_6_para.my = 4;
  update_rk2_0_6_para.mz = 1;
  update_rk2_0_6_para.state_oldh_31 = &state_old->h_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dxh_31 = &tend1->dxh_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dyh_31 = &tend1->dyh_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.state_oldu_31 = &state_old->u_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dxu_31 = &tend1->dxu_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dyu_31 = &tend1->dyu_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.state_oldv_31 = &state_old->v_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dxv_31 = &tend1->dxv_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend1dyv_31 = &tend1->dyv_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dxh_31 = &tend2->dxh_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dyh_31 = &tend2->dyh_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dxu_31 = &tend2->dxu_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dyu_31 = &tend2->dyu_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dxv_31 = &tend2->dxv_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.tend2dyv_31 = &tend2->dyv_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.state_newh_31 = &state_new->h_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.state_newu_31 = &state_new->u_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.state_newv_31 = &state_new->v_31[0][Proc.lat_hw -1+update_rk2_0_6_para.oy][Proc.lon_hw -1+update_rk2_0_6_para.ox];
  update_rk2_0_6_para.dt = dt;
  
  athread_spawn(update_rk2_0_6, &update_rk2_0_6_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_31Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_31[0][0][0], &Proc.FieldReq[async_h_31], async_h_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_31[0][0][0], &Proc.FieldReq[async_u_31], async_u_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_31[0][0][0], &Proc.FieldReq[async_v_31], async_v_31, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state_old->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state_old->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state_old->v_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_7_info update_rk2_0_7_para;
  update_rk2_0_7_para.lz = 1 - 0;
  update_rk2_0_7_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_7_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_7_para.oz = 0;
  update_rk2_0_7_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_7_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_7_para.hx = 1;
  update_rk2_0_7_para.hy = 1;
  update_rk2_0_7_para.hz = 0;
  update_rk2_0_7_para.bx = 128;
  update_rk2_0_7_para.by = 4;
  update_rk2_0_7_para.bz = 1;
  update_rk2_0_7_para.mx = 1;
  update_rk2_0_7_para.my = 4;
  update_rk2_0_7_para.mz = 1;
  update_rk2_0_7_para.state_oldh_32 = &state_old->h_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dxh_32 = &tend1->dxh_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dyh_32 = &tend1->dyh_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.state_oldu_32 = &state_old->u_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dxu_32 = &tend1->dxu_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dyu_32 = &tend1->dyu_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.state_oldv_32 = &state_old->v_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dxv_32 = &tend1->dxv_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend1dyv_32 = &tend1->dyv_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dxh_32 = &tend2->dxh_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dyh_32 = &tend2->dyh_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dxu_32 = &tend2->dxu_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dyu_32 = &tend2->dyu_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dxv_32 = &tend2->dxv_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.tend2dyv_32 = &tend2->dyv_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.state_newh_32 = &state_new->h_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.state_newu_32 = &state_new->u_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.state_newv_32 = &state_new->v_32[0][Proc.lat_hw -1+update_rk2_0_7_para.oy][Proc.lon_hw -1+update_rk2_0_7_para.ox];
  update_rk2_0_7_para.dt = dt;
  
  athread_spawn(update_rk2_0_7, &update_rk2_0_7_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_32Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_32[0][0][0], &Proc.FieldReq[async_h_32], async_h_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_32[0][0][0], &Proc.FieldReq[async_u_32], async_u_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_32[0][0][0], &Proc.FieldReq[async_v_32], async_v_32, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state_old->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state_old->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state_old->v_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk2_0_8_info update_rk2_0_8_para;
  update_rk2_0_8_para.lz = 1 - 0;
  update_rk2_0_8_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk2_0_8_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk2_0_8_para.oz = 0;
  update_rk2_0_8_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk2_0_8_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk2_0_8_para.hx = 1;
  update_rk2_0_8_para.hy = 1;
  update_rk2_0_8_para.hz = 0;
  update_rk2_0_8_para.bx = 128;
  update_rk2_0_8_para.by = 4;
  update_rk2_0_8_para.bz = 1;
  update_rk2_0_8_para.mx = 1;
  update_rk2_0_8_para.my = 4;
  update_rk2_0_8_para.mz = 1;
  update_rk2_0_8_para.state_oldh_33 = &state_old->h_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dxh_33 = &tend1->dxh_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dyh_33 = &tend1->dyh_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.state_oldu_33 = &state_old->u_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dxu_33 = &tend1->dxu_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dyu_33 = &tend1->dyu_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.state_oldv_33 = &state_old->v_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dxv_33 = &tend1->dxv_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend1dyv_33 = &tend1->dyv_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dxh_33 = &tend2->dxh_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dyh_33 = &tend2->dyh_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dxu_33 = &tend2->dxu_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dyu_33 = &tend2->dyu_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dxv_33 = &tend2->dxv_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.tend2dyv_33 = &tend2->dyv_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.state_newh_33 = &state_new->h_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.state_newu_33 = &state_new->u_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.state_newv_33 = &state_new->v_33[0][Proc.lat_hw -1+update_rk2_0_8_para.oy][Proc.lon_hw -1+update_rk2_0_8_para.ox];
  update_rk2_0_8_para.dt = dt;
  
  athread_spawn(update_rk2_0_8, &update_rk2_0_8_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk2h_33Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_33[0][0][0], &Proc.FieldReq[async_h_33], async_h_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_33[0][0][0], &Proc.FieldReq[async_u_33], async_u_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_33[0][0][0], &Proc.FieldReq[async_v_33], async_v_33, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void update_rk3_0(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridTendField* tend1, struct HybridTendField* tend2, struct HybridTendField* tend3, double dt)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_11], async_h_11, state_old->h_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_11], async_u_11, state_old->u_11, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_11], async_v_11, state_old->v_11, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_0_info update_rk3_0_0_para;
  update_rk3_0_0_para.lz = 1 - 0;
  update_rk3_0_0_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_0_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_0_para.oz = 0;
  update_rk3_0_0_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_0_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_0_para.hx = 1;
  update_rk3_0_0_para.hy = 1;
  update_rk3_0_0_para.hz = 0;
  update_rk3_0_0_para.bx = 128;
  update_rk3_0_0_para.by = 2;
  update_rk3_0_0_para.bz = 1;
  update_rk3_0_0_para.mx = 1;
  update_rk3_0_0_para.my = 4;
  update_rk3_0_0_para.mz = 1;
  update_rk3_0_0_para.state_oldh_11 = &state_old->h_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dxh_11 = &tend1->dxh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dyh_11 = &tend1->dyh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.state_oldu_11 = &state_old->u_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dxu_11 = &tend1->dxu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dyu_11 = &tend1->dyu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.state_oldv_11 = &state_old->v_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dxv_11 = &tend1->dxv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend1dyv_11 = &tend1->dyv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dxh_11 = &tend2->dxh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dyh_11 = &tend2->dyh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dxh_11 = &tend3->dxh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dyh_11 = &tend3->dyh_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dxu_11 = &tend2->dxu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dyu_11 = &tend2->dyu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dxu_11 = &tend3->dxu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dyu_11 = &tend3->dyu_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dxv_11 = &tend2->dxv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend2dyv_11 = &tend2->dyv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dxv_11 = &tend3->dxv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.tend3dyv_11 = &tend3->dyv_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.state_newh_11 = &state_new->h_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.state_newu_11 = &state_new->u_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.state_newv_11 = &state_new->v_11[0][Proc.lat_hw -1+update_rk3_0_0_para.oy][Proc.lon_hw -1+update_rk3_0_0_para.ox];
  update_rk3_0_0_para.dt = dt;
  
  athread_spawn(update_rk3_0_0, &update_rk3_0_0_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_11Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_11[0][0][0], &Proc.FieldReq[async_h_11], async_h_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_11[0][0][0], &Proc.FieldReq[async_u_11], async_u_11, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_11[0][0][0], &Proc.FieldReq[async_v_11], async_v_11, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_12], async_h_12, state_old->h_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_12], async_u_12, state_old->u_12, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_12], async_v_12, state_old->v_12, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_1_info update_rk3_0_1_para;
  update_rk3_0_1_para.lz = 1 - 0;
  update_rk3_0_1_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_1_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_1_para.oz = 0;
  update_rk3_0_1_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_1_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_1_para.hx = 1;
  update_rk3_0_1_para.hy = 1;
  update_rk3_0_1_para.hz = 0;
  update_rk3_0_1_para.bx = 128;
  update_rk3_0_1_para.by = 2;
  update_rk3_0_1_para.bz = 1;
  update_rk3_0_1_para.mx = 1;
  update_rk3_0_1_para.my = 4;
  update_rk3_0_1_para.mz = 1;
  update_rk3_0_1_para.state_oldh_12 = &state_old->h_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dxh_12 = &tend1->dxh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dyh_12 = &tend1->dyh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.state_oldu_12 = &state_old->u_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dxu_12 = &tend1->dxu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dyu_12 = &tend1->dyu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.state_oldv_12 = &state_old->v_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dxv_12 = &tend1->dxv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend1dyv_12 = &tend1->dyv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dxh_12 = &tend2->dxh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dyh_12 = &tend2->dyh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dxh_12 = &tend3->dxh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dyh_12 = &tend3->dyh_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dxu_12 = &tend2->dxu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dyu_12 = &tend2->dyu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dxu_12 = &tend3->dxu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dyu_12 = &tend3->dyu_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dxv_12 = &tend2->dxv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend2dyv_12 = &tend2->dyv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dxv_12 = &tend3->dxv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.tend3dyv_12 = &tend3->dyv_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.state_newh_12 = &state_new->h_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.state_newu_12 = &state_new->u_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.state_newv_12 = &state_new->v_12[0][Proc.lat_hw -1+update_rk3_0_1_para.oy][Proc.lon_hw -1+update_rk3_0_1_para.ox];
  update_rk3_0_1_para.dt = dt;
  
  athread_spawn(update_rk3_0_1, &update_rk3_0_1_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_12Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_12[0][0][0], &Proc.FieldReq[async_h_12], async_h_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_12[0][0][0], &Proc.FieldReq[async_u_12], async_u_12, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_12[0][0][0], &Proc.FieldReq[async_v_12], async_v_12, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_13], async_h_13, state_old->h_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_13], async_u_13, state_old->u_13, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_13], async_v_13, state_old->v_13, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_2_info update_rk3_0_2_para;
  update_rk3_0_2_para.lz = 1 - 0;
  update_rk3_0_2_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_2_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_2_para.oz = 0;
  update_rk3_0_2_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_2_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_2_para.hx = 1;
  update_rk3_0_2_para.hy = 1;
  update_rk3_0_2_para.hz = 0;
  update_rk3_0_2_para.bx = 128;
  update_rk3_0_2_para.by = 2;
  update_rk3_0_2_para.bz = 1;
  update_rk3_0_2_para.mx = 1;
  update_rk3_0_2_para.my = 4;
  update_rk3_0_2_para.mz = 1;
  update_rk3_0_2_para.state_oldh_13 = &state_old->h_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dxh_13 = &tend1->dxh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dyh_13 = &tend1->dyh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.state_oldu_13 = &state_old->u_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dxu_13 = &tend1->dxu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dyu_13 = &tend1->dyu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.state_oldv_13 = &state_old->v_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dxv_13 = &tend1->dxv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend1dyv_13 = &tend1->dyv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dxh_13 = &tend2->dxh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dyh_13 = &tend2->dyh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dxh_13 = &tend3->dxh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dyh_13 = &tend3->dyh_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dxu_13 = &tend2->dxu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dyu_13 = &tend2->dyu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dxu_13 = &tend3->dxu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dyu_13 = &tend3->dyu_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dxv_13 = &tend2->dxv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend2dyv_13 = &tend2->dyv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dxv_13 = &tend3->dxv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.tend3dyv_13 = &tend3->dyv_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.state_newh_13 = &state_new->h_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.state_newu_13 = &state_new->u_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.state_newv_13 = &state_new->v_13[0][Proc.lat_hw -1+update_rk3_0_2_para.oy][Proc.lon_hw -1+update_rk3_0_2_para.ox];
  update_rk3_0_2_para.dt = dt;
  
  athread_spawn(update_rk3_0_2, &update_rk3_0_2_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_13Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_13[0][0][0], &Proc.FieldReq[async_h_13], async_h_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_13[0][0][0], &Proc.FieldReq[async_u_13], async_u_13, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_13[0][0][0], &Proc.FieldReq[async_v_13], async_v_13, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_21], async_h_21, state_old->h_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_21], async_u_21, state_old->u_21, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_21], async_v_21, state_old->v_21, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_3_info update_rk3_0_3_para;
  update_rk3_0_3_para.lz = 1 - 0;
  update_rk3_0_3_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_3_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_3_para.oz = 0;
  update_rk3_0_3_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_3_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_3_para.hx = 1;
  update_rk3_0_3_para.hy = 1;
  update_rk3_0_3_para.hz = 0;
  update_rk3_0_3_para.bx = 128;
  update_rk3_0_3_para.by = 2;
  update_rk3_0_3_para.bz = 1;
  update_rk3_0_3_para.mx = 1;
  update_rk3_0_3_para.my = 4;
  update_rk3_0_3_para.mz = 1;
  update_rk3_0_3_para.state_oldh_21 = &state_old->h_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dxh_21 = &tend1->dxh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dyh_21 = &tend1->dyh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.state_oldu_21 = &state_old->u_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dxu_21 = &tend1->dxu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dyu_21 = &tend1->dyu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.state_oldv_21 = &state_old->v_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dxv_21 = &tend1->dxv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend1dyv_21 = &tend1->dyv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dxh_21 = &tend2->dxh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dyh_21 = &tend2->dyh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dxh_21 = &tend3->dxh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dyh_21 = &tend3->dyh_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dxu_21 = &tend2->dxu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dyu_21 = &tend2->dyu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dxu_21 = &tend3->dxu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dyu_21 = &tend3->dyu_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dxv_21 = &tend2->dxv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend2dyv_21 = &tend2->dyv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dxv_21 = &tend3->dxv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.tend3dyv_21 = &tend3->dyv_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.state_newh_21 = &state_new->h_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.state_newu_21 = &state_new->u_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.state_newv_21 = &state_new->v_21[0][Proc.lat_hw -1+update_rk3_0_3_para.oy][Proc.lon_hw -1+update_rk3_0_3_para.ox];
  update_rk3_0_3_para.dt = dt;
  
  athread_spawn(update_rk3_0_3, &update_rk3_0_3_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_21Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_21[0][0][0], &Proc.FieldReq[async_h_21], async_h_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_21[0][0][0], &Proc.FieldReq[async_u_21], async_u_21, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_21[0][0][0], &Proc.FieldReq[async_v_21], async_v_21, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_22], async_h_22, state_old->h_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_22], async_u_22, state_old->u_22, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_22], async_v_22, state_old->v_22, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_4_info update_rk3_0_4_para;
  update_rk3_0_4_para.lz = 1 - 0;
  update_rk3_0_4_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_4_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_4_para.oz = 0;
  update_rk3_0_4_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_4_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_4_para.hx = 1;
  update_rk3_0_4_para.hy = 1;
  update_rk3_0_4_para.hz = 0;
  update_rk3_0_4_para.bx = 128;
  update_rk3_0_4_para.by = 2;
  update_rk3_0_4_para.bz = 1;
  update_rk3_0_4_para.mx = 1;
  update_rk3_0_4_para.my = 4;
  update_rk3_0_4_para.mz = 1;
  update_rk3_0_4_para.state_oldh_22 = &state_old->h_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dxh_22 = &tend1->dxh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dyh_22 = &tend1->dyh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.state_oldu_22 = &state_old->u_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dxu_22 = &tend1->dxu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dyu_22 = &tend1->dyu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.state_oldv_22 = &state_old->v_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dxv_22 = &tend1->dxv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend1dyv_22 = &tend1->dyv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dxh_22 = &tend2->dxh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dyh_22 = &tend2->dyh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dxh_22 = &tend3->dxh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dyh_22 = &tend3->dyh_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dxu_22 = &tend2->dxu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dyu_22 = &tend2->dyu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dxu_22 = &tend3->dxu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dyu_22 = &tend3->dyu_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dxv_22 = &tend2->dxv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend2dyv_22 = &tend2->dyv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dxv_22 = &tend3->dxv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.tend3dyv_22 = &tend3->dyv_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.state_newh_22 = &state_new->h_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.state_newu_22 = &state_new->u_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.state_newv_22 = &state_new->v_22[0][Proc.lat_hw -1+update_rk3_0_4_para.oy][Proc.lon_hw -1+update_rk3_0_4_para.ox];
  update_rk3_0_4_para.dt = dt;
  
  athread_spawn(update_rk3_0_4, &update_rk3_0_4_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_22Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_22[0][0][0], &Proc.FieldReq[async_h_22], async_h_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_22[0][0][0], &Proc.FieldReq[async_u_22], async_u_22, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_22[0][0][0], &Proc.FieldReq[async_v_22], async_v_22, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_23], async_h_23, state_old->h_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_23], async_u_23, state_old->u_23, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_23], async_v_23, state_old->v_23, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_5_info update_rk3_0_5_para;
  update_rk3_0_5_para.lz = 1 - 0;
  update_rk3_0_5_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_5_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_5_para.oz = 0;
  update_rk3_0_5_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_5_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_5_para.hx = 1;
  update_rk3_0_5_para.hy = 1;
  update_rk3_0_5_para.hz = 0;
  update_rk3_0_5_para.bx = 128;
  update_rk3_0_5_para.by = 2;
  update_rk3_0_5_para.bz = 1;
  update_rk3_0_5_para.mx = 1;
  update_rk3_0_5_para.my = 4;
  update_rk3_0_5_para.mz = 1;
  update_rk3_0_5_para.state_oldh_23 = &state_old->h_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dxh_23 = &tend1->dxh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dyh_23 = &tend1->dyh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.state_oldu_23 = &state_old->u_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dxu_23 = &tend1->dxu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dyu_23 = &tend1->dyu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.state_oldv_23 = &state_old->v_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dxv_23 = &tend1->dxv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend1dyv_23 = &tend1->dyv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dxh_23 = &tend2->dxh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dyh_23 = &tend2->dyh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dxh_23 = &tend3->dxh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dyh_23 = &tend3->dyh_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dxu_23 = &tend2->dxu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dyu_23 = &tend2->dyu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dxu_23 = &tend3->dxu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dyu_23 = &tend3->dyu_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dxv_23 = &tend2->dxv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend2dyv_23 = &tend2->dyv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dxv_23 = &tend3->dxv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.tend3dyv_23 = &tend3->dyv_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.state_newh_23 = &state_new->h_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.state_newu_23 = &state_new->u_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.state_newv_23 = &state_new->v_23[0][Proc.lat_hw -1+update_rk3_0_5_para.oy][Proc.lon_hw -1+update_rk3_0_5_para.ox];
  update_rk3_0_5_para.dt = dt;
  
  athread_spawn(update_rk3_0_5, &update_rk3_0_5_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_23Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_23[0][0][0], &Proc.FieldReq[async_h_23], async_h_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_23[0][0][0], &Proc.FieldReq[async_u_23], async_u_23, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_23[0][0][0], &Proc.FieldReq[async_v_23], async_v_23, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_31], async_h_31, state_old->h_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_31], async_u_31, state_old->u_31, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_31], async_v_31, state_old->v_31, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_6_info update_rk3_0_6_para;
  update_rk3_0_6_para.lz = 1 - 0;
  update_rk3_0_6_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_6_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_6_para.oz = 0;
  update_rk3_0_6_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_6_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_6_para.hx = 1;
  update_rk3_0_6_para.hy = 1;
  update_rk3_0_6_para.hz = 0;
  update_rk3_0_6_para.bx = 128;
  update_rk3_0_6_para.by = 2;
  update_rk3_0_6_para.bz = 1;
  update_rk3_0_6_para.mx = 1;
  update_rk3_0_6_para.my = 4;
  update_rk3_0_6_para.mz = 1;
  update_rk3_0_6_para.state_oldh_31 = &state_old->h_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dxh_31 = &tend1->dxh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dyh_31 = &tend1->dyh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.state_oldu_31 = &state_old->u_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dxu_31 = &tend1->dxu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dyu_31 = &tend1->dyu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.state_oldv_31 = &state_old->v_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dxv_31 = &tend1->dxv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend1dyv_31 = &tend1->dyv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dxh_31 = &tend2->dxh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dyh_31 = &tend2->dyh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dxh_31 = &tend3->dxh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dyh_31 = &tend3->dyh_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dxu_31 = &tend2->dxu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dyu_31 = &tend2->dyu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dxu_31 = &tend3->dxu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dyu_31 = &tend3->dyu_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dxv_31 = &tend2->dxv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend2dyv_31 = &tend2->dyv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dxv_31 = &tend3->dxv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.tend3dyv_31 = &tend3->dyv_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.state_newh_31 = &state_new->h_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.state_newu_31 = &state_new->u_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.state_newv_31 = &state_new->v_31[0][Proc.lat_hw -1+update_rk3_0_6_para.oy][Proc.lon_hw -1+update_rk3_0_6_para.ox];
  update_rk3_0_6_para.dt = dt;
  
  athread_spawn(update_rk3_0_6, &update_rk3_0_6_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_31Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_31[0][0][0], &Proc.FieldReq[async_h_31], async_h_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_31[0][0][0], &Proc.FieldReq[async_u_31], async_u_31, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_31[0][0][0], &Proc.FieldReq[async_v_31], async_v_31, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_32], async_h_32, state_old->h_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_32], async_u_32, state_old->u_32, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_32], async_v_32, state_old->v_32, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_7_info update_rk3_0_7_para;
  update_rk3_0_7_para.lz = 1 - 0;
  update_rk3_0_7_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_7_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_7_para.oz = 0;
  update_rk3_0_7_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_7_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_7_para.hx = 1;
  update_rk3_0_7_para.hy = 1;
  update_rk3_0_7_para.hz = 0;
  update_rk3_0_7_para.bx = 128;
  update_rk3_0_7_para.by = 2;
  update_rk3_0_7_para.bz = 1;
  update_rk3_0_7_para.mx = 1;
  update_rk3_0_7_para.my = 4;
  update_rk3_0_7_para.mz = 1;
  update_rk3_0_7_para.state_oldh_32 = &state_old->h_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dxh_32 = &tend1->dxh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dyh_32 = &tend1->dyh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.state_oldu_32 = &state_old->u_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dxu_32 = &tend1->dxu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dyu_32 = &tend1->dyu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.state_oldv_32 = &state_old->v_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dxv_32 = &tend1->dxv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend1dyv_32 = &tend1->dyv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dxh_32 = &tend2->dxh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dyh_32 = &tend2->dyh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dxh_32 = &tend3->dxh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dyh_32 = &tend3->dyh_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dxu_32 = &tend2->dxu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dyu_32 = &tend2->dyu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dxu_32 = &tend3->dxu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dyu_32 = &tend3->dyu_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dxv_32 = &tend2->dxv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend2dyv_32 = &tend2->dyv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dxv_32 = &tend3->dxv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.tend3dyv_32 = &tend3->dyv_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.state_newh_32 = &state_new->h_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.state_newu_32 = &state_new->u_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.state_newv_32 = &state_new->v_32[0][Proc.lat_hw -1+update_rk3_0_7_para.oy][Proc.lon_hw -1+update_rk3_0_7_para.ox];
  update_rk3_0_7_para.dt = dt;
  
  athread_spawn(update_rk3_0_7, &update_rk3_0_7_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_32Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_32[0][0][0], &Proc.FieldReq[async_h_32], async_h_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_32[0][0][0], &Proc.FieldReq[async_u_32], async_u_32, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_32[0][0][0], &Proc.FieldReq[async_v_32], async_v_32, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWaitCS(Proc,&Proc.FieldReq[async_h_33], async_h_33, state_old->h_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_u_33], async_u_33, state_old->u_33, true, true, true, true, false, true, false);
  HaloWaitCS(Proc,&Proc.FieldReq[async_v_33], async_v_33, state_old->v_33, true, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  update_rk3_0_8_info update_rk3_0_8_para;
  update_rk3_0_8_para.lz = 1 - 0;
  update_rk3_0_8_para.ly = MIN(Proc.lat_end+1, 600) - MAX(Proc.lat_beg, 0);
  update_rk3_0_8_para.lx = MIN(Proc.lon_end+1, 600) - MAX(Proc.lon_beg, 0);
  update_rk3_0_8_para.oz = 0;
  update_rk3_0_8_para.oy = MAX(Proc.lat_beg, 0) - Proc.lat_beg;
  update_rk3_0_8_para.ox = MAX(Proc.lon_beg, 0) - Proc.lon_beg;
  update_rk3_0_8_para.hx = 1;
  update_rk3_0_8_para.hy = 1;
  update_rk3_0_8_para.hz = 0;
  update_rk3_0_8_para.bx = 128;
  update_rk3_0_8_para.by = 2;
  update_rk3_0_8_para.bz = 1;
  update_rk3_0_8_para.mx = 1;
  update_rk3_0_8_para.my = 4;
  update_rk3_0_8_para.mz = 1;
  update_rk3_0_8_para.state_oldh_33 = &state_old->h_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dxh_33 = &tend1->dxh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dyh_33 = &tend1->dyh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.state_oldu_33 = &state_old->u_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dxu_33 = &tend1->dxu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dyu_33 = &tend1->dyu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.state_oldv_33 = &state_old->v_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dxv_33 = &tend1->dxv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend1dyv_33 = &tend1->dyv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dxh_33 = &tend2->dxh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dyh_33 = &tend2->dyh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dxh_33 = &tend3->dxh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dyh_33 = &tend3->dyh_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dxu_33 = &tend2->dxu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dyu_33 = &tend2->dyu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dxu_33 = &tend3->dxu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dyu_33 = &tend3->dyu_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dxv_33 = &tend2->dxv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend2dyv_33 = &tend2->dyv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dxv_33 = &tend3->dxv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.tend3dyv_33 = &tend3->dyv_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.state_newh_33 = &state_new->h_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.state_newu_33 = &state_new->u_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.state_newv_33 = &state_new->v_33[0][Proc.lat_hw -1+update_rk3_0_8_para.oy][Proc.lon_hw -1+update_rk3_0_8_para.ox];
  update_rk3_0_8_para.dt = dt;
  
  athread_spawn(update_rk3_0_8, &update_rk3_0_8_para);
  athread_join();
  
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  update_rk3h_33Time += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHaloCS_2d_D(Proc, &state_new->h_33[0][0][0], &Proc.FieldReq[async_h_33], async_h_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->u_33[0][0][0], &Proc.FieldReq[async_u_33], async_u_33, true, true, true, false, true, false);
  UpdateHaloCS_2d_D(Proc, &state_new->v_33[0][0][0], &Proc.FieldReq[async_v_33], async_v_33, true, true, true, false, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void MeshInit_0(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh)
{
  MeshInit(&global_mesh[0]);
  MeshInit_cp(&global_mesh[0], &mesh[0]);
}

void OperatorInit_0(struct HybridStateField* state, struct HybridMeshField* mesh)
{
  spaceOperatorInit_0(&state[Timeinfo.oldt], &mesh[0]);
}

void Rk3_0(struct HybridStateField* state, struct HybridTendField* tend, struct HybridMeshField* mesh, double dt)
{
  updateState_0(&state[Timeinfo.oldt], &state[2]);
  mcvupdateXY_0(&state[Timeinfo.oldt], &tend[Timeinfo.oldt], &mesh[0]);
  update_rk1_0(&state[Timeinfo.oldt], &state[2], &tend[Timeinfo.oldt], dt);
  mcvupdateXY_0(&state[Timeinfo.oldt], &tend[Timeinfo.newt], &mesh[0]);
  update_rk2_0(&state[Timeinfo.oldt], &state[2], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], dt);
  mcvupdateXY_0(&state[Timeinfo.oldt], &tend[2], &mesh[0]);
  update_rk3_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], &tend[2], dt);
}

int main(int argc, char **argv)
{
  int size,rank;
  double tottime_beg, tottime_end;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ProcInit_CubedSphere_Domain(&Proc,MPI_COMM_WORLD,size,rank,ProcLatNum,ProcLon,ProcLat,ncell,nlev);
  
  //Proc Ngb Init
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 599){
    ProcInit_CubedSphere_Ngb(&Proc,MPI_COMM_WORLD,size,rank,ProcLon[0],ProcLat[0],ncell,nlev,1, 1, 0, 2, DTYPE_DOUBLE, 157);
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
  
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 599)
  {
    MeshInit_0(global_mesh, mesh);
    OperatorInit_0(state, mesh);
    tottime_beg = MPI_Wtime();
    for (int t = 0; t < 30; t += 1){
      Rk3_0(state, tend, mesh, 7.5);
    }
    tottime_end = MPI_Wtime();
  }
  
  ProfilingOutPut(Proc.id);
  if (Proc.id == 0) printf("Time is %.8f\n",tottime_end-tottime_beg);
}