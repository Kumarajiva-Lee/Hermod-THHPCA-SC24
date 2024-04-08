#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>


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

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_11[0][j][i] = state_old->h_11[0][j][i];
      state_new->u_11[0][j][i] = state_old->u_11[0][j][i];
      state_new->v_11[0][j][i] = state_old->v_11[0][j][i];
      state_new->jab_11[0][j][i] = state_old->jab_11[0][j][i];
      state_new->jab5_11[0][j][i] = state_old->jab5_11[0][j][i];
      state_new->jab6_11[0][j][i] = state_old->jab6_11[0][j][i];
      state_new->jab7_11[0][j][i] = state_old->jab7_11[0][j][i];
      state_new->uc_11[0][j][i] = state_old->uc_11[0][j][i];
      state_new->vc_11[0][j][i] = state_old->vc_11[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_12[0][j][i] = state_old->h_12[0][j][i];
      state_new->u_12[0][j][i] = state_old->u_12[0][j][i];
      state_new->v_12[0][j][i] = state_old->v_12[0][j][i];
      state_new->jab_12[0][j][i] = state_old->jab_12[0][j][i];
      state_new->jab5_12[0][j][i] = state_old->jab5_12[0][j][i];
      state_new->jab6_12[0][j][i] = state_old->jab6_12[0][j][i];
      state_new->jab7_12[0][j][i] = state_old->jab7_12[0][j][i];
      state_new->uc_12[0][j][i] = state_old->uc_12[0][j][i];
      state_new->vc_12[0][j][i] = state_old->vc_12[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_13[0][j][i] = state_old->h_13[0][j][i];
      state_new->u_13[0][j][i] = state_old->u_13[0][j][i];
      state_new->v_13[0][j][i] = state_old->v_13[0][j][i];
      state_new->jab_13[0][j][i] = state_old->jab_13[0][j][i];
      state_new->jab5_13[0][j][i] = state_old->jab5_13[0][j][i];
      state_new->jab6_13[0][j][i] = state_old->jab6_13[0][j][i];
      state_new->jab7_13[0][j][i] = state_old->jab7_13[0][j][i];
      state_new->uc_13[0][j][i] = state_old->uc_13[0][j][i];
      state_new->vc_13[0][j][i] = state_old->vc_13[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_21[0][j][i] = state_old->h_21[0][j][i];
      state_new->u_21[0][j][i] = state_old->u_21[0][j][i];
      state_new->v_21[0][j][i] = state_old->v_21[0][j][i];
      state_new->jab_21[0][j][i] = state_old->jab_21[0][j][i];
      state_new->jab5_21[0][j][i] = state_old->jab5_21[0][j][i];
      state_new->jab6_21[0][j][i] = state_old->jab6_21[0][j][i];
      state_new->jab7_21[0][j][i] = state_old->jab7_21[0][j][i];
      state_new->uc_21[0][j][i] = state_old->uc_21[0][j][i];
      state_new->vc_21[0][j][i] = state_old->vc_21[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_22[0][j][i] = state_old->h_22[0][j][i];
      state_new->u_22[0][j][i] = state_old->u_22[0][j][i];
      state_new->v_22[0][j][i] = state_old->v_22[0][j][i];
      state_new->jab_22[0][j][i] = state_old->jab_22[0][j][i];
      state_new->jab5_22[0][j][i] = state_old->jab5_22[0][j][i];
      state_new->jab6_22[0][j][i] = state_old->jab6_22[0][j][i];
      state_new->jab7_22[0][j][i] = state_old->jab7_22[0][j][i];
      state_new->uc_22[0][j][i] = state_old->uc_22[0][j][i];
      state_new->vc_22[0][j][i] = state_old->vc_22[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_23[0][j][i] = state_old->h_23[0][j][i];
      state_new->u_23[0][j][i] = state_old->u_23[0][j][i];
      state_new->v_23[0][j][i] = state_old->v_23[0][j][i];
      state_new->jab_23[0][j][i] = state_old->jab_23[0][j][i];
      state_new->jab5_23[0][j][i] = state_old->jab5_23[0][j][i];
      state_new->jab6_23[0][j][i] = state_old->jab6_23[0][j][i];
      state_new->jab7_23[0][j][i] = state_old->jab7_23[0][j][i];
      state_new->uc_23[0][j][i] = state_old->uc_23[0][j][i];
      state_new->vc_23[0][j][i] = state_old->vc_23[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_31[0][j][i] = state_old->h_31[0][j][i];
      state_new->u_31[0][j][i] = state_old->u_31[0][j][i];
      state_new->v_31[0][j][i] = state_old->v_31[0][j][i];
      state_new->jab_31[0][j][i] = state_old->jab_31[0][j][i];
      state_new->jab5_31[0][j][i] = state_old->jab5_31[0][j][i];
      state_new->jab6_31[0][j][i] = state_old->jab6_31[0][j][i];
      state_new->jab7_31[0][j][i] = state_old->jab7_31[0][j][i];
      state_new->uc_31[0][j][i] = state_old->uc_31[0][j][i];
      state_new->vc_31[0][j][i] = state_old->vc_31[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_32[0][j][i] = state_old->h_32[0][j][i];
      state_new->u_32[0][j][i] = state_old->u_32[0][j][i];
      state_new->v_32[0][j][i] = state_old->v_32[0][j][i];
      state_new->jab_32[0][j][i] = state_old->jab_32[0][j][i];
      state_new->jab5_32[0][j][i] = state_old->jab5_32[0][j][i];
      state_new->jab6_32[0][j][i] = state_old->jab6_32[0][j][i];
      state_new->jab7_32[0][j][i] = state_old->jab7_32[0][j][i];
      state_new->uc_32[0][j][i] = state_old->uc_32[0][j][i];
      state_new->vc_32[0][j][i] = state_old->vc_32[0][j][i];
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 1)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_33[0][j][i] = state_old->h_33[0][j][i];
      state_new->u_33[0][j][i] = state_old->u_33[0][j][i];
      state_new->v_33[0][j][i] = state_old->v_33[0][j][i];
      state_new->jab_33[0][j][i] = state_old->jab_33[0][j][i];
      state_new->jab5_33[0][j][i] = state_old->jab5_33[0][j][i];
      state_new->jab6_33[0][j][i] = state_old->jab6_33[0][j][i];
      state_new->jab7_33[0][j][i] = state_old->jab7_33[0][j][i];
      state_new->uc_33[0][j][i] = state_old->uc_33[0][j][i];
      state_new->vc_33[0][j][i] = state_old->vc_33[0][j][i];
    }
  }
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
  int j, i;
  double ulocal, dql, dqr, dfl, dfr, sl, sr, sc, vlocal;
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_31[0][j][(i - 1)]) + sqrt((((state->jab5_31[0][j][(i - 1)] * 9.80616) * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)])));
      tend->fh_1[0][j][i] = ((((state->uc_31[0][j][(i - 1)] * state->h_31[0][j][(i - 1)]) + (state->uc_11[0][j][i] * state->h_11[0][j][i])) + (ulocal * (state->h_31[0][j][(i - 1)] - state->h_11[0][j][i]))) / 2);
      tend->fu_1[0][j][i] = ((((((9.80616 * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)]) + (0.5 * ((state->u_31[0][j][(i - 1)] * state->uc_31[0][j][(i - 1)]) + (state->v_31[0][j][(i - 1)] * state->vc_31[0][j][(i - 1)])))) + (((9.80616 * state->h_11[0][j][i]) / state->jab_11[0][j][i]) + (0.5 * ((state->u_11[0][j][i] * state->uc_11[0][j][i]) + (state->v_11[0][j][i] * state->vc_11[0][j][i]))))) + (ulocal * (state->u_31[0][j][(i - 1)] - state->u_11[0][j][i]))) / 2);
      tend->fv_1[0][j][i] = (((0.0 + 0.0) + (ulocal * (state->v_31[0][j][(i - 1)] - state->v_11[0][j][i]))) / 2);
      tend->qh_1[0][j][i] = ((state->h_31[0][j][(i - 1)] + state->h_11[0][j][i]) / 2);
      tend->qu_1[0][j][i] = ((state->u_31[0][j][(i - 1)] + state->u_11[0][j][i]) / 2);
      tend->qv_1[0][j][i] = ((state->v_31[0][j][(i - 1)] + state->v_11[0][j][i]) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_31[0][j][(i - 1)]) + sqrt((((state->jab5_31[0][j][(i - 1)] * 9.80616) * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)])));
      dql = ((state->h_11[0][j][(i - 1)] - (4 * state->h_21[0][j][(i - 1)])) + (3 * state->h_31[0][j][(i - 1)]));
      dqr = -(((3 * state->h_11[0][j][i]) - (4 * state->h_21[0][j][i])) + state->h_31[0][j][i]);
      dfl = (((state->uc_11[0][j][(i - 1)] * state->h_11[0][j][(i - 1)]) - (4 * (state->uc_21[0][j][(i - 1)] * state->h_21[0][j][(i - 1)]))) + (3 * (state->uc_31[0][j][(i - 1)] * state->h_31[0][j][(i - 1)])));
      dfr = -(((3 * (state->uc_11[0][j][i] * state->h_11[0][j][i])) - (4 * (state->uc_21[0][j][i] * state->h_21[0][j][i]))) + (state->uc_31[0][j][i] * state->h_31[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_31[0][j][(i - 1)]) + sqrt((((state->jab5_31[0][j][(i - 1)] * 9.80616) * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)])));
      dql = ((state->u_11[0][j][(i - 1)] - (4 * state->u_21[0][j][(i - 1)])) + (3 * state->u_31[0][j][(i - 1)]));
      dqr = -(((3 * state->u_11[0][j][i]) - (4 * state->u_21[0][j][i])) + state->u_31[0][j][i]);
      dfl = (((((9.80616 * state->h_11[0][j][(i - 1)]) / state->jab_11[0][j][(i - 1)]) + (0.5 * ((state->u_11[0][j][(i - 1)] * state->uc_11[0][j][(i - 1)]) + (state->v_11[0][j][(i - 1)] * state->vc_11[0][j][(i - 1)])))) - (4 * (((9.80616 * state->h_21[0][j][(i - 1)]) / state->jab_21[0][j][(i - 1)]) + (0.5 * ((state->u_21[0][j][(i - 1)] * state->uc_21[0][j][(i - 1)]) + (state->v_21[0][j][(i - 1)] * state->vc_21[0][j][(i - 1)])))))) + (3 * (((9.80616 * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)]) + (0.5 * ((state->u_31[0][j][(i - 1)] * state->uc_31[0][j][(i - 1)]) + (state->v_31[0][j][(i - 1)] * state->vc_31[0][j][(i - 1)]))))));
      dfr = -(((3 * (((9.80616 * state->h_11[0][j][i]) / state->jab_11[0][j][i]) + (0.5 * ((state->u_11[0][j][i] * state->uc_11[0][j][i]) + (state->v_11[0][j][i] * state->vc_11[0][j][i]))))) - (4 * (((9.80616 * state->h_21[0][j][i]) / state->jab_21[0][j][i]) + (0.5 * ((state->u_21[0][j][i] * state->uc_21[0][j][i]) + (state->v_21[0][j][i] * state->vc_21[0][j][i])))))) + (((9.80616 * state->h_31[0][j][i]) / state->jab_31[0][j][i]) + (0.5 * ((state->u_31[0][j][i] * state->uc_31[0][j][i]) + (state->v_31[0][j][i] * state->vc_31[0][j][i])))));
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_31[0][j][(i - 1)]) + sqrt((((state->jab5_31[0][j][(i - 1)] * 9.80616) * state->h_31[0][j][(i - 1)]) / state->jab_31[0][j][(i - 1)])));
      dql = ((state->v_11[0][j][(i - 1)] - (4 * state->v_21[0][j][(i - 1)])) + (3 * state->v_31[0][j][(i - 1)]));
      dqr = -(((3 * state->v_11[0][j][i]) - (4 * state->v_21[0][j][i])) + state->v_31[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qv_2[0][j][i] / 16679.814955336966);
      sr = (tend->qv_2[0][j][(i + 1)] / 16679.814955336966);
      sc = (((((3.0 * (tend->qv_1[0][j][(i + 1)] - tend->qv_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qv_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_11[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dxh_21[0][j][i] = (((((-3.0 * (tend->fh_1[0][j][(i + 1)] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_31[0][j][i] = (-tend->fh_2[0][j][(i + 1)] / 16679.814955336966);
      tend->dxu_11[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) + (state->vc_11[0][j][i] * sl));
      tend->dxu_21[0][j][i] = ((((((3.0 * (tend->fu_1[0][j][(i + 1)] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) + (state->vc_21[0][j][i] * sc));
      tend->dxu_31[0][j][i] = ((-tend->fu_2[0][j][(i + 1)] / 16679.814955336966) + (state->vc_31[0][j][i] * sr));
      tend->dxv_11[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) - (state->uc_11[0][j][i] * sl));
      tend->dxv_21[0][j][i] = ((((((-3.0 * (tend->fv_1[0][j][(i + 1)] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) - (state->uc_21[0][j][i] * sc));
      tend->dxv_31[0][j][i] = ((-tend->fv_2[0][j][(i + 1)] / 16679.814955336966) - (state->uc_31[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_32[0][j][(i - 1)]) + sqrt((((state->jab5_32[0][j][(i - 1)] * 9.80616) * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)])));
      tend->fh_1[0][j][i] = ((((state->uc_32[0][j][(i - 1)] * state->h_32[0][j][(i - 1)]) + (state->uc_12[0][j][i] * state->h_12[0][j][i])) + (ulocal * (state->h_32[0][j][(i - 1)] - state->h_12[0][j][i]))) / 2);
      tend->fu_1[0][j][i] = ((((((9.80616 * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)]) + (0.5 * ((state->u_32[0][j][(i - 1)] * state->uc_32[0][j][(i - 1)]) + (state->v_32[0][j][(i - 1)] * state->vc_32[0][j][(i - 1)])))) + (((9.80616 * state->h_12[0][j][i]) / state->jab_12[0][j][i]) + (0.5 * ((state->u_12[0][j][i] * state->uc_12[0][j][i]) + (state->v_12[0][j][i] * state->vc_12[0][j][i]))))) + (ulocal * (state->u_32[0][j][(i - 1)] - state->u_12[0][j][i]))) / 2);
      tend->fv_1[0][j][i] = (((0.0 + 0.0) + (ulocal * (state->v_32[0][j][(i - 1)] - state->v_12[0][j][i]))) / 2);
      tend->qh_1[0][j][i] = ((state->h_32[0][j][(i - 1)] + state->h_12[0][j][i]) / 2);
      tend->qu_1[0][j][i] = ((state->u_32[0][j][(i - 1)] + state->u_12[0][j][i]) / 2);
      tend->qv_1[0][j][i] = ((state->v_32[0][j][(i - 1)] + state->v_12[0][j][i]) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_32[0][j][(i - 1)]) + sqrt((((state->jab5_32[0][j][(i - 1)] * 9.80616) * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)])));
      dql = ((state->h_12[0][j][(i - 1)] - (4 * state->h_22[0][j][(i - 1)])) + (3 * state->h_32[0][j][(i - 1)]));
      dqr = -(((3 * state->h_12[0][j][i]) - (4 * state->h_22[0][j][i])) + state->h_32[0][j][i]);
      dfl = (((state->uc_12[0][j][(i - 1)] * state->h_12[0][j][(i - 1)]) - (4 * (state->uc_22[0][j][(i - 1)] * state->h_22[0][j][(i - 1)]))) + (3 * (state->uc_32[0][j][(i - 1)] * state->h_32[0][j][(i - 1)])));
      dfr = -(((3 * (state->uc_12[0][j][i] * state->h_12[0][j][i])) - (4 * (state->uc_22[0][j][i] * state->h_22[0][j][i]))) + (state->uc_32[0][j][i] * state->h_32[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_32[0][j][(i - 1)]) + sqrt((((state->jab5_32[0][j][(i - 1)] * 9.80616) * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)])));
      dql = ((state->u_12[0][j][(i - 1)] - (4 * state->u_22[0][j][(i - 1)])) + (3 * state->u_32[0][j][(i - 1)]));
      dqr = -(((3 * state->u_12[0][j][i]) - (4 * state->u_22[0][j][i])) + state->u_32[0][j][i]);
      dfl = (((((9.80616 * state->h_12[0][j][(i - 1)]) / state->jab_12[0][j][(i - 1)]) + (0.5 * ((state->u_12[0][j][(i - 1)] * state->uc_12[0][j][(i - 1)]) + (state->v_12[0][j][(i - 1)] * state->vc_12[0][j][(i - 1)])))) - (4 * (((9.80616 * state->h_22[0][j][(i - 1)]) / state->jab_22[0][j][(i - 1)]) + (0.5 * ((state->u_22[0][j][(i - 1)] * state->uc_22[0][j][(i - 1)]) + (state->v_22[0][j][(i - 1)] * state->vc_22[0][j][(i - 1)])))))) + (3 * (((9.80616 * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)]) + (0.5 * ((state->u_32[0][j][(i - 1)] * state->uc_32[0][j][(i - 1)]) + (state->v_32[0][j][(i - 1)] * state->vc_32[0][j][(i - 1)]))))));
      dfr = -(((3 * (((9.80616 * state->h_12[0][j][i]) / state->jab_12[0][j][i]) + (0.5 * ((state->u_12[0][j][i] * state->uc_12[0][j][i]) + (state->v_12[0][j][i] * state->vc_12[0][j][i]))))) - (4 * (((9.80616 * state->h_22[0][j][i]) / state->jab_22[0][j][i]) + (0.5 * ((state->u_22[0][j][i] * state->uc_22[0][j][i]) + (state->v_22[0][j][i] * state->vc_22[0][j][i])))))) + (((9.80616 * state->h_32[0][j][i]) / state->jab_32[0][j][i]) + (0.5 * ((state->u_32[0][j][i] * state->uc_32[0][j][i]) + (state->v_32[0][j][i] * state->vc_32[0][j][i])))));
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_32[0][j][(i - 1)]) + sqrt((((state->jab5_32[0][j][(i - 1)] * 9.80616) * state->h_32[0][j][(i - 1)]) / state->jab_32[0][j][(i - 1)])));
      dql = ((state->v_12[0][j][(i - 1)] - (4 * state->v_22[0][j][(i - 1)])) + (3 * state->v_32[0][j][(i - 1)]));
      dqr = -(((3 * state->v_12[0][j][i]) - (4 * state->v_22[0][j][i])) + state->v_32[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qv_2[0][j][i] / 16679.814955336966);
      sr = (tend->qv_2[0][j][(i + 1)] / 16679.814955336966);
      sc = (((((3.0 * (tend->qv_1[0][j][(i + 1)] - tend->qv_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qv_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_12[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dxh_22[0][j][i] = (((((-3.0 * (tend->fh_1[0][j][(i + 1)] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_32[0][j][i] = (-tend->fh_2[0][j][(i + 1)] / 16679.814955336966);
      tend->dxu_12[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) + (state->vc_12[0][j][i] * sl));
      tend->dxu_22[0][j][i] = ((((((3.0 * (tend->fu_1[0][j][(i + 1)] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) + (state->vc_22[0][j][i] * sc));
      tend->dxu_32[0][j][i] = ((-tend->fu_2[0][j][(i + 1)] / 16679.814955336966) + (state->vc_32[0][j][i] * sr));
      tend->dxv_12[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) - (state->uc_12[0][j][i] * sl));
      tend->dxv_22[0][j][i] = ((((((-3.0 * (tend->fv_1[0][j][(i + 1)] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) - (state->uc_22[0][j][i] * sc));
      tend->dxv_32[0][j][i] = ((-tend->fv_2[0][j][(i + 1)] / 16679.814955336966) - (state->uc_32[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_33[0][j][(i - 1)]) + sqrt((((state->jab5_33[0][j][(i - 1)] * 9.80616) * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)])));
      tend->fh_1[0][j][i] = ((((state->uc_33[0][j][(i - 1)] * state->h_33[0][j][(i - 1)]) + (state->uc_13[0][j][i] * state->h_13[0][j][i])) + (ulocal * (state->h_33[0][j][(i - 1)] - state->h_13[0][j][i]))) / 2);
      tend->fu_1[0][j][i] = ((((((9.80616 * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)]) + (0.5 * ((state->u_33[0][j][(i - 1)] * state->uc_33[0][j][(i - 1)]) + (state->v_33[0][j][(i - 1)] * state->vc_33[0][j][(i - 1)])))) + (((9.80616 * state->h_13[0][j][i]) / state->jab_13[0][j][i]) + (0.5 * ((state->u_13[0][j][i] * state->uc_13[0][j][i]) + (state->v_13[0][j][i] * state->vc_13[0][j][i]))))) + (ulocal * (state->u_33[0][j][(i - 1)] - state->u_13[0][j][i]))) / 2);
      tend->fv_1[0][j][i] = (((0.0 + 0.0) + (ulocal * (state->v_33[0][j][(i - 1)] - state->v_13[0][j][i]))) / 2);
      tend->qh_1[0][j][i] = ((state->h_33[0][j][(i - 1)] + state->h_13[0][j][i]) / 2);
      tend->qu_1[0][j][i] = ((state->u_33[0][j][(i - 1)] + state->u_13[0][j][i]) / 2);
      tend->qv_1[0][j][i] = ((state->v_33[0][j][(i - 1)] + state->v_13[0][j][i]) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_33[0][j][(i - 1)]) + sqrt((((state->jab5_33[0][j][(i - 1)] * 9.80616) * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)])));
      dql = ((state->h_13[0][j][(i - 1)] - (4 * state->h_23[0][j][(i - 1)])) + (3 * state->h_33[0][j][(i - 1)]));
      dqr = -(((3 * state->h_13[0][j][i]) - (4 * state->h_23[0][j][i])) + state->h_33[0][j][i]);
      dfl = (((state->uc_13[0][j][(i - 1)] * state->h_13[0][j][(i - 1)]) - (4 * (state->uc_23[0][j][(i - 1)] * state->h_23[0][j][(i - 1)]))) + (3 * (state->uc_33[0][j][(i - 1)] * state->h_33[0][j][(i - 1)])));
      dfr = -(((3 * (state->uc_13[0][j][i] * state->h_13[0][j][i])) - (4 * (state->uc_23[0][j][i] * state->h_23[0][j][i]))) + (state->uc_33[0][j][i] * state->h_33[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_33[0][j][(i - 1)]) + sqrt((((state->jab5_33[0][j][(i - 1)] * 9.80616) * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)])));
      dql = ((state->u_13[0][j][(i - 1)] - (4 * state->u_23[0][j][(i - 1)])) + (3 * state->u_33[0][j][(i - 1)]));
      dqr = -(((3 * state->u_13[0][j][i]) - (4 * state->u_23[0][j][i])) + state->u_33[0][j][i]);
      dfl = (((((9.80616 * state->h_13[0][j][(i - 1)]) / state->jab_13[0][j][(i - 1)]) + (0.5 * ((state->u_13[0][j][(i - 1)] * state->uc_13[0][j][(i - 1)]) + (state->v_13[0][j][(i - 1)] * state->vc_13[0][j][(i - 1)])))) - (4 * (((9.80616 * state->h_23[0][j][(i - 1)]) / state->jab_23[0][j][(i - 1)]) + (0.5 * ((state->u_23[0][j][(i - 1)] * state->uc_23[0][j][(i - 1)]) + (state->v_23[0][j][(i - 1)] * state->vc_23[0][j][(i - 1)])))))) + (3 * (((9.80616 * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)]) + (0.5 * ((state->u_33[0][j][(i - 1)] * state->uc_33[0][j][(i - 1)]) + (state->v_33[0][j][(i - 1)] * state->vc_33[0][j][(i - 1)]))))));
      dfr = -(((3 * (((9.80616 * state->h_13[0][j][i]) / state->jab_13[0][j][i]) + (0.5 * ((state->u_13[0][j][i] * state->uc_13[0][j][i]) + (state->v_13[0][j][i] * state->vc_13[0][j][i]))))) - (4 * (((9.80616 * state->h_23[0][j][i]) / state->jab_23[0][j][i]) + (0.5 * ((state->u_23[0][j][i] * state->uc_23[0][j][i]) + (state->v_23[0][j][i] * state->vc_23[0][j][i])))))) + (((9.80616 * state->h_33[0][j][i]) / state->jab_33[0][j][i]) + (0.5 * ((state->u_33[0][j][i] * state->uc_33[0][j][i]) + (state->v_33[0][j][i] * state->vc_33[0][j][i])))));
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      ulocal = (fabs(state->uc_33[0][j][(i - 1)]) + sqrt((((state->jab5_33[0][j][(i - 1)] * 9.80616) * state->h_33[0][j][(i - 1)]) / state->jab_33[0][j][(i - 1)])));
      dql = ((state->v_13[0][j][(i - 1)] - (4 * state->v_23[0][j][(i - 1)])) + (3 * state->v_33[0][j][(i - 1)]));
      dqr = -(((3 * state->v_13[0][j][i]) - (4 * state->v_23[0][j][i])) + state->v_33[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qv_2[0][j][i] / 16679.814955336966);
      sr = (tend->qv_2[0][j][(i + 1)] / 16679.814955336966);
      sc = (((((3.0 * (tend->qv_1[0][j][(i + 1)] - tend->qv_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qv_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_13[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dxh_23[0][j][i] = (((((-3.0 * (tend->fh_1[0][j][(i + 1)] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][j][(i + 1)] / 16679.814955336966) / 4.0));
      tend->dxh_33[0][j][i] = (-tend->fh_2[0][j][(i + 1)] / 16679.814955336966);
      tend->dxu_13[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) + (state->vc_13[0][j][i] * sl));
      tend->dxu_23[0][j][i] = ((((((3.0 * (tend->fu_1[0][j][(i + 1)] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) + (state->vc_23[0][j][i] * sc));
      tend->dxu_33[0][j][i] = ((-tend->fu_2[0][j][(i + 1)] / 16679.814955336966) + (state->vc_33[0][j][i] * sr));
      tend->dxv_13[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) - (state->uc_13[0][j][i] * sl));
      tend->dxv_23[0][j][i] = ((((((-3.0 * (tend->fv_1[0][j][(i + 1)] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][j][(i + 1)] / 16679.814955336966) / 4.0)) - (state->uc_23[0][j][i] * sc));
      tend->dxv_33[0][j][i] = ((-tend->fv_2[0][j][(i + 1)] / 16679.814955336966) - (state->uc_33[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_13[0][(j - 1)][i]) + sqrt((((state->jab7_13[0][(j - 1)][i] * 9.80616) * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i])));
      tend->fh_1[0][j][i] = (((state->vc_13[0][(j - 1)][i] * state->h_13[0][(j - 1)][i]) + (state->vc_11[0][j][i] * state->h_11[0][j][i])) + ((vlocal * (state->h_13[0][(j - 1)][i] - state->h_11[0][j][i])) / 2.0));
      tend->fu_1[0][j][i] = ((0.0 + 0.0) + ((vlocal * (state->u_13[0][(j - 1)][i] - state->u_11[0][j][i])) / 2.0));
      tend->fv_1[0][j][i] = (((((9.80616 * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i]) + (0.5 * ((state->u_13[0][(j - 1)][i] * state->uc_13[0][(j - 1)][i]) + (state->v_13[0][(j - 1)][i] * state->vc_13[0][(j - 1)][i])))) + (((9.80616 * state->h_11[0][j][i]) / state->jab_11[0][j][i]) + (0.5 * ((state->u_11[0][j][i] * state->uc_11[0][j][i]) + (state->v_11[0][j][i] * state->vc_11[0][j][i]))))) + ((vlocal * (state->v_13[0][(j - 1)][i] - state->v_11[0][j][i])) / 2.0));
      tend->qh_1[0][j][i] = ((state->h_13[0][(j - 1)][i] + state->h_11[0][j][i]) / 2.0);
      tend->qu_1[0][j][i] = ((state->u_13[0][(j - 1)][i] + state->u_11[0][j][i]) / 2.0);
      tend->qv_1[0][j][i] = ((state->v_13[0][(j - 1)][i] + state->v_11[0][j][i]) / 2.0);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_13[0][(j - 1)][i]) + sqrt((((state->jab7_13[0][(j - 1)][i] * 9.80616) * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i])));
      dql = ((state->h_11[0][(j - 1)][i] - (4 * state->h_12[0][(j - 1)][i])) + (3 * state->h_13[0][(j - 1)][i]));
      dqr = -(((3 * state->h_11[0][j][i]) - (4 * state->h_12[0][j][i])) + state->h_13[0][j][i]);
      dfl = (((state->vc_11[0][(j - 1)][i] * state->h_11[0][(j - 1)][i]) - (4 * (state->vc_12[0][(j - 1)][i] * state->h_12[0][(j - 1)][i]))) + (3 * (state->vc_13[0][(j - 1)][i] * state->h_13[0][(j - 1)][i])));
      dfr = -(((3 * (state->vc_11[0][j][i] * state->h_11[0][j][i])) - (4 * (state->vc_12[0][j][i] * state->h_12[0][j][i]))) + (state->vc_13[0][j][i] * state->h_13[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_13[0][(j - 1)][i]) + sqrt((((state->jab7_13[0][(j - 1)][i] * 9.80616) * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i])));
      dql = ((state->u_11[0][(j - 1)][i] - (4 * state->u_12[0][(j - 1)][i])) + (3 * state->u_13[0][(j - 1)][i]));
      dqr = -(((3 * state->u_11[0][j][i]) - (4 * state->u_12[0][j][i])) + state->u_13[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_13[0][(j - 1)][i]) + sqrt((((state->jab7_13[0][(j - 1)][i] * 9.80616) * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i])));
      dql = ((state->v_11[0][(j - 1)][i] - (4 * state->v_12[0][(j - 1)][i])) + (3 * state->v_13[0][(j - 1)][i]));
      dqr = -(((3 * state->v_11[0][j][i]) - (4 * state->v_12[0][j][i])) + state->v_13[0][j][i]);
      dfl = (((((9.80616 * state->h_11[0][(j - 1)][i]) / state->jab_11[0][(j - 1)][i]) + (0.5 * ((state->u_11[0][(j - 1)][i] * state->uc_11[0][(j - 1)][i]) + (state->v_11[0][(j - 1)][i] * state->vc_11[0][(j - 1)][i])))) - (4 * (((9.80616 * state->h_12[0][(j - 1)][i]) / state->jab_12[0][(j - 1)][i]) + (0.5 * ((state->u_12[0][(j - 1)][i] * state->uc_12[0][(j - 1)][i]) + (state->v_12[0][(j - 1)][i] * state->vc_12[0][(j - 1)][i])))))) + (3 * (((9.80616 * state->h_13[0][(j - 1)][i]) / state->jab_13[0][(j - 1)][i]) + (0.5 * ((state->u_13[0][(j - 1)][i] * state->uc_13[0][(j - 1)][i]) + (state->v_13[0][(j - 1)][i] * state->vc_13[0][(j - 1)][i]))))));
      dfr = -(((3 * (((9.80616 * state->h_11[0][j][i]) / state->jab_11[0][j][i]) + (0.5 * ((state->u_11[0][j][i] * state->uc_11[0][j][i]) + (state->v_11[0][j][i] * state->vc_11[0][j][i]))))) - (4 * (((9.80616 * state->h_12[0][j][i]) / state->jab_12[0][j][i]) + (0.5 * ((state->u_12[0][j][i] * state->uc_12[0][j][i]) + (state->v_12[0][j][i] * state->vc_12[0][j][i])))))) + (((9.80616 * state->h_13[0][j][i]) / state->jab_13[0][j][i]) + (0.5 * ((state->u_13[0][j][i] * state->uc_13[0][j][i]) + (state->v_13[0][j][i] * state->vc_13[0][j][i])))));
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qu_2[0][j][i] / 16679.814955336966);
      sr = (tend->qu_2[0][(j + 1)][i] / 16679.814955336966);
      sc = (((((3.0 * (tend->qu_1[0][(j + 1)][i] - tend->qu_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qu_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_11[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dyh_12[0][j][i] = (((((-3 * (tend->fh_1[0][(j + 1)][i] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_13[0][j][i] = (-tend->fh_2[0][(j + 1)][i] / 16679.814955336966);
      tend->dyu_11[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) - (state->vc_11[0][j][i] * sl));
      tend->dyu_12[0][j][i] = ((((((-3.0 * (tend->fu_1[0][(j + 1)][i] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->vc_12[0][j][i] * sc));
      tend->dyu_13[0][j][i] = ((-tend->fu_2[0][(j + 1)][i] / 16679.814955336966) - (state->vc_13[0][j][i] * sr));
      tend->dyv_11[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) + (state->uc_11[0][j][i] * sl));
      tend->dyv_12[0][j][i] = ((((((-3.0 * (tend->fv_1[0][(j + 1)][i] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->uc_12[0][j][i] * sc));
      tend->dyv_13[0][j][i] = ((-tend->fv_2[0][(j + 1)][i] / 16679.814955336966) - (state->uc_13[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_23[0][(j - 1)][i]) + sqrt((((state->jab7_23[0][(j - 1)][i] * 9.80616) * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i])));
      tend->fh_1[0][j][i] = (((state->vc_23[0][(j - 1)][i] * state->h_23[0][(j - 1)][i]) + (state->vc_21[0][j][i] * state->h_21[0][j][i])) + ((vlocal * (state->h_23[0][(j - 1)][i] - state->h_21[0][j][i])) / 2.0));
      tend->fu_1[0][j][i] = ((0.0 + 0.0) + ((vlocal * (state->u_23[0][(j - 1)][i] - state->u_21[0][j][i])) / 2.0));
      tend->fv_1[0][j][i] = (((((9.80616 * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i]) + (0.5 * ((state->u_23[0][(j - 1)][i] * state->uc_23[0][(j - 1)][i]) + (state->v_23[0][(j - 1)][i] * state->vc_23[0][(j - 1)][i])))) + (((9.80616 * state->h_21[0][j][i]) / state->jab_21[0][j][i]) + (0.5 * ((state->u_21[0][j][i] * state->uc_21[0][j][i]) + (state->v_21[0][j][i] * state->vc_21[0][j][i]))))) + ((vlocal * (state->v_23[0][(j - 1)][i] - state->v_21[0][j][i])) / 2.0));
      tend->qh_1[0][j][i] = ((state->h_23[0][(j - 1)][i] + state->h_21[0][j][i]) / 2.0);
      tend->qu_1[0][j][i] = ((state->u_23[0][(j - 1)][i] + state->u_21[0][j][i]) / 2.0);
      tend->qv_1[0][j][i] = ((state->v_23[0][(j - 1)][i] + state->v_21[0][j][i]) / 2.0);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_23[0][(j - 1)][i]) + sqrt((((state->jab7_23[0][(j - 1)][i] * 9.80616) * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i])));
      dql = ((state->h_21[0][(j - 1)][i] - (4 * state->h_22[0][(j - 1)][i])) + (3 * state->h_23[0][(j - 1)][i]));
      dqr = -(((3 * state->h_21[0][j][i]) - (4 * state->h_22[0][j][i])) + state->h_23[0][j][i]);
      dfl = (((state->vc_21[0][(j - 1)][i] * state->h_21[0][(j - 1)][i]) - (4 * (state->vc_22[0][(j - 1)][i] * state->h_22[0][(j - 1)][i]))) + (3 * (state->vc_23[0][(j - 1)][i] * state->h_23[0][(j - 1)][i])));
      dfr = -(((3 * (state->vc_21[0][j][i] * state->h_21[0][j][i])) - (4 * (state->vc_22[0][j][i] * state->h_22[0][j][i]))) + (state->vc_23[0][j][i] * state->h_23[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_23[0][(j - 1)][i]) + sqrt((((state->jab7_23[0][(j - 1)][i] * 9.80616) * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i])));
      dql = ((state->u_21[0][(j - 1)][i] - (4 * state->u_22[0][(j - 1)][i])) + (3 * state->u_23[0][(j - 1)][i]));
      dqr = -(((3 * state->u_21[0][j][i]) - (4 * state->u_22[0][j][i])) + state->u_23[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_23[0][(j - 1)][i]) + sqrt((((state->jab7_23[0][(j - 1)][i] * 9.80616) * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i])));
      dql = ((state->v_21[0][(j - 1)][i] - (4 * state->v_22[0][(j - 1)][i])) + (3 * state->v_23[0][(j - 1)][i]));
      dqr = -(((3 * state->v_21[0][j][i]) - (4 * state->v_22[0][j][i])) + state->v_23[0][j][i]);
      dfl = (((((9.80616 * state->h_21[0][(j - 1)][i]) / state->jab_21[0][(j - 1)][i]) + (0.5 * ((state->u_21[0][(j - 1)][i] * state->uc_21[0][(j - 1)][i]) + (state->v_21[0][(j - 1)][i] * state->vc_21[0][(j - 1)][i])))) - (4 * (((9.80616 * state->h_22[0][(j - 1)][i]) / state->jab_22[0][(j - 1)][i]) + (0.5 * ((state->u_22[0][(j - 1)][i] * state->uc_22[0][(j - 1)][i]) + (state->v_22[0][(j - 1)][i] * state->vc_22[0][(j - 1)][i])))))) + (3 * (((9.80616 * state->h_23[0][(j - 1)][i]) / state->jab_23[0][(j - 1)][i]) + (0.5 * ((state->u_23[0][(j - 1)][i] * state->uc_23[0][(j - 1)][i]) + (state->v_23[0][(j - 1)][i] * state->vc_23[0][(j - 1)][i]))))));
      dfr = -(((3 * (((9.80616 * state->h_21[0][j][i]) / state->jab_21[0][j][i]) + (0.5 * ((state->u_21[0][j][i] * state->uc_21[0][j][i]) + (state->v_21[0][j][i] * state->vc_21[0][j][i]))))) - (4 * (((9.80616 * state->h_22[0][j][i]) / state->jab_22[0][j][i]) + (0.5 * ((state->u_22[0][j][i] * state->uc_22[0][j][i]) + (state->v_22[0][j][i] * state->vc_22[0][j][i])))))) + (((9.80616 * state->h_23[0][j][i]) / state->jab_23[0][j][i]) + (0.5 * ((state->u_23[0][j][i] * state->uc_23[0][j][i]) + (state->v_23[0][j][i] * state->vc_23[0][j][i])))));
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qu_2[0][j][i] / 16679.814955336966);
      sr = (tend->qu_2[0][(j + 1)][i] / 16679.814955336966);
      sc = (((((3.0 * (tend->qu_1[0][(j + 1)][i] - tend->qu_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qu_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_11[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dyh_12[0][j][i] = (((((-3 * (tend->fh_1[0][(j + 1)][i] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_13[0][j][i] = (-tend->fh_2[0][(j + 1)][i] / 16679.814955336966);
      tend->dyu_11[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) - (state->vc_21[0][j][i] * sl));
      tend->dyu_12[0][j][i] = ((((((-3.0 * (tend->fu_1[0][(j + 1)][i] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->vc_22[0][j][i] * sc));
      tend->dyu_13[0][j][i] = ((-tend->fu_2[0][(j + 1)][i] / 16679.814955336966) - (state->vc_23[0][j][i] * sr));
      tend->dyv_11[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) + (state->uc_21[0][j][i] * sl));
      tend->dyv_12[0][j][i] = ((((((-3.0 * (tend->fv_1[0][(j + 1)][i] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->uc_22[0][j][i] * sc));
      tend->dyv_13[0][j][i] = ((-tend->fv_2[0][(j + 1)][i] / 16679.814955336966) - (state->uc_23[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_33[0][(j - 1)][i]) + sqrt((((state->jab7_33[0][(j - 1)][i] * 9.80616) * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i])));
      tend->fh_1[0][j][i] = (((state->vc_33[0][(j - 1)][i] * state->h_33[0][(j - 1)][i]) + (state->vc_31[0][j][i] * state->h_31[0][j][i])) + ((vlocal * (state->h_33[0][(j - 1)][i] - state->h_31[0][j][i])) / 2.0));
      tend->fu_1[0][j][i] = ((0.0 + 0.0) + ((vlocal * (state->u_33[0][(j - 1)][i] - state->u_31[0][j][i])) / 2.0));
      tend->fv_1[0][j][i] = (((((9.80616 * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i]) + (0.5 * ((state->u_33[0][(j - 1)][i] * state->uc_33[0][(j - 1)][i]) + (state->v_33[0][(j - 1)][i] * state->vc_33[0][(j - 1)][i])))) + (((9.80616 * state->h_31[0][j][i]) / state->jab_31[0][j][i]) + (0.5 * ((state->u_31[0][j][i] * state->uc_31[0][j][i]) + (state->v_31[0][j][i] * state->vc_31[0][j][i]))))) + ((vlocal * (state->v_33[0][(j - 1)][i] - state->v_31[0][j][i])) / 2.0));
      tend->qh_1[0][j][i] = ((state->h_33[0][(j - 1)][i] + state->h_31[0][j][i]) / 2.0);
      tend->qu_1[0][j][i] = ((state->u_33[0][(j - 1)][i] + state->u_31[0][j][i]) / 2.0);
      tend->qv_1[0][j][i] = ((state->v_33[0][(j - 1)][i] + state->v_31[0][j][i]) / 2.0);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_33[0][(j - 1)][i]) + sqrt((((state->jab7_33[0][(j - 1)][i] * 9.80616) * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i])));
      dql = ((state->h_31[0][(j - 1)][i] - (4 * state->h_32[0][(j - 1)][i])) + (3 * state->h_33[0][(j - 1)][i]));
      dqr = -(((3 * state->h_31[0][j][i]) - (4 * state->h_32[0][j][i])) + state->h_33[0][j][i]);
      dfl = (((state->vc_31[0][(j - 1)][i] * state->h_31[0][(j - 1)][i]) - (4 * (state->vc_32[0][(j - 1)][i] * state->h_32[0][(j - 1)][i]))) + (3 * (state->vc_33[0][(j - 1)][i] * state->h_33[0][(j - 1)][i])));
      dfr = -(((3 * (state->vc_31[0][j][i] * state->h_31[0][j][i])) - (4 * (state->vc_32[0][j][i] * state->h_32[0][j][i]))) + (state->vc_33[0][j][i] * state->h_33[0][j][i]));
      tend->fh_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qh_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_33[0][(j - 1)][i]) + sqrt((((state->jab7_33[0][(j - 1)][i] * 9.80616) * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i])));
      dql = ((state->u_31[0][(j - 1)][i] - (4 * state->u_32[0][(j - 1)][i])) + (3 * state->u_33[0][(j - 1)][i]));
      dqr = -(((3 * state->u_31[0][j][i]) - (4 * state->u_32[0][j][i])) + state->u_33[0][j][i]);
      dfl = 0.0;
      dfr = 0.0;
      tend->fu_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qu_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      vlocal = (fabs(state->vc_33[0][(j - 1)][i]) + sqrt((((state->jab7_33[0][(j - 1)][i] * 9.80616) * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i])));
      dql = ((state->v_31[0][(j - 1)][i] - (4 * state->v_32[0][(j - 1)][i])) + (3 * state->v_33[0][(j - 1)][i]));
      dqr = -(((3 * state->v_31[0][j][i]) - (4 * state->v_32[0][j][i])) + state->v_33[0][j][i]);
      dfl = (((((9.80616 * state->h_31[0][(j - 1)][i]) / state->jab_31[0][(j - 1)][i]) + (0.5 * ((state->u_31[0][(j - 1)][i] * state->uc_31[0][(j - 1)][i]) + (state->v_31[0][(j - 1)][i] * state->vc_31[0][(j - 1)][i])))) - (4 * (((9.80616 * state->h_32[0][(j - 1)][i]) / state->jab_32[0][(j - 1)][i]) + (0.5 * ((state->u_32[0][(j - 1)][i] * state->uc_32[0][(j - 1)][i]) + (state->v_32[0][(j - 1)][i] * state->vc_32[0][(j - 1)][i])))))) + (3 * (((9.80616 * state->h_33[0][(j - 1)][i]) / state->jab_33[0][(j - 1)][i]) + (0.5 * ((state->u_33[0][(j - 1)][i] * state->uc_33[0][(j - 1)][i]) + (state->v_33[0][(j - 1)][i] * state->vc_33[0][(j - 1)][i]))))));
      dfr = -(((3 * (((9.80616 * state->h_31[0][j][i]) / state->jab_31[0][j][i]) + (0.5 * ((state->u_31[0][j][i] * state->uc_31[0][j][i]) + (state->v_31[0][j][i] * state->vc_31[0][j][i]))))) - (4 * (((9.80616 * state->h_32[0][j][i]) / state->jab_32[0][j][i]) + (0.5 * ((state->u_32[0][j][i] * state->uc_32[0][j][i]) + (state->v_32[0][j][i] * state->vc_32[0][j][i])))))) + (((9.80616 * state->h_33[0][j][i]) / state->jab_33[0][j][i]) + (0.5 * ((state->u_33[0][j][i] * state->uc_33[0][j][i]) + (state->v_33[0][j][i] * state->vc_33[0][j][i])))));
      tend->fv_2[0][j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
      tend->qv_2[0][j][i] = ((dql + dqr) / 2);
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      sl = (tend->qu_2[0][j][i] / 16679.814955336966);
      sr = (tend->qu_2[0][(j + 1)][i] / 16679.814955336966);
      sc = (((((3.0 * (tend->qu_1[0][(j + 1)][i] - tend->qu_1[0][j][i])) / 16679.814955336966) / 2.0) - ((tend->qu_2[0][j][i] / 16679.814955336966) / 4.0)) - ((tend->qu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_11[0][j][i] = (-tend->fh_2[0][j][i] / 16679.814955336966);
      tend->dyh_12[0][j][i] = (((((-3 * (tend->fh_1[0][(j + 1)][i] - tend->fh_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fh_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fh_2[0][(j + 1)][i] / 16679.814955336966) / 4.0));
      tend->dyh_13[0][j][i] = (-tend->fh_2[0][(j + 1)][i] / 16679.814955336966);
      tend->dyu_11[0][j][i] = ((-tend->fu_2[0][j][i] / 16679.814955336966) - (state->vc_31[0][j][i] * sl));
      tend->dyu_12[0][j][i] = ((((((-3.0 * (tend->fu_1[0][(j + 1)][i] - tend->fu_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fu_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fu_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->vc_32[0][j][i] * sc));
      tend->dyu_13[0][j][i] = ((-tend->fu_2[0][(j + 1)][i] / 16679.814955336966) - (state->vc_33[0][j][i] * sr));
      tend->dyv_11[0][j][i] = ((-tend->fv_2[0][j][i] / 16679.814955336966) + (state->uc_31[0][j][i] * sl));
      tend->dyv_12[0][j][i] = ((((((-3.0 * (tend->fv_1[0][(j + 1)][i] - tend->fv_1[0][j][i])) / 16679.814955336966) / 2.0) + ((tend->fv_2[0][j][i] / 16679.814955336966) / 4.0)) + ((tend->fv_2[0][(j + 1)][i] / 16679.814955336966) / 4.0)) - (state->uc_32[0][j][i] * sc));
      tend->dyv_13[0][j][i] = ((-tend->fv_2[0][(j + 1)][i] / 16679.814955336966) - (state->uc_33[0][j][i] * sr));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_11[0][j][i] = (state_old->h_11[0][j][i] + ((tend1->dxh_11[0][j][i] + tend1->dyh_11[0][j][i]) * dt));
      state_new->u_11[0][j][i] = (state_old->u_11[0][j][i] + ((tend1->dxu_11[0][j][i] + tend1->dyu_11[0][j][i]) * dt));
      state_new->v_11[0][j][i] = (state_old->v_11[0][j][i] + ((tend1->dxv_11[0][j][i] + tend1->dyv_11[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_12[0][j][i] = (state_old->h_12[0][j][i] + ((tend1->dxh_12[0][j][i] + tend1->dyh_12[0][j][i]) * dt));
      state_new->u_12[0][j][i] = (state_old->u_12[0][j][i] + ((tend1->dxu_12[0][j][i] + tend1->dyu_12[0][j][i]) * dt));
      state_new->v_12[0][j][i] = (state_old->v_12[0][j][i] + ((tend1->dxv_12[0][j][i] + tend1->dyv_12[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_13[0][j][i] = (state_old->h_13[0][j][i] + ((tend1->dxh_13[0][j][i] + tend1->dyh_13[0][j][i]) * dt));
      state_new->u_13[0][j][i] = (state_old->u_13[0][j][i] + ((tend1->dxu_13[0][j][i] + tend1->dyu_13[0][j][i]) * dt));
      state_new->v_13[0][j][i] = (state_old->v_13[0][j][i] + ((tend1->dxv_13[0][j][i] + tend1->dyv_13[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_21[0][j][i] = (state_old->h_21[0][j][i] + ((tend1->dxh_21[0][j][i] + tend1->dyh_21[0][j][i]) * dt));
      state_new->u_21[0][j][i] = (state_old->u_21[0][j][i] + ((tend1->dxu_21[0][j][i] + tend1->dyu_21[0][j][i]) * dt));
      state_new->v_21[0][j][i] = (state_old->v_21[0][j][i] + ((tend1->dxv_21[0][j][i] + tend1->dyv_21[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_22[0][j][i] = (state_old->h_22[0][j][i] + ((tend1->dxh_22[0][j][i] + tend1->dyh_22[0][j][i]) * dt));
      state_new->u_22[0][j][i] = (state_old->u_22[0][j][i] + ((tend1->dxu_22[0][j][i] + tend1->dyu_22[0][j][i]) * dt));
      state_new->v_22[0][j][i] = (state_old->v_22[0][j][i] + ((tend1->dxv_22[0][j][i] + tend1->dyv_22[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_23[0][j][i] = (state_old->h_23[0][j][i] + ((tend1->dxh_23[0][j][i] + tend1->dyh_23[0][j][i]) * dt));
      state_new->u_23[0][j][i] = (state_old->u_23[0][j][i] + ((tend1->dxu_23[0][j][i] + tend1->dyu_23[0][j][i]) * dt));
      state_new->v_23[0][j][i] = (state_old->v_23[0][j][i] + ((tend1->dxv_23[0][j][i] + tend1->dyv_23[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_31[0][j][i] = (state_old->h_31[0][j][i] + ((tend1->dxh_31[0][j][i] + tend1->dyh_31[0][j][i]) * dt));
      state_new->u_31[0][j][i] = (state_old->u_31[0][j][i] + ((tend1->dxu_31[0][j][i] + tend1->dyu_31[0][j][i]) * dt));
      state_new->v_31[0][j][i] = (state_old->v_31[0][j][i] + ((tend1->dxv_31[0][j][i] + tend1->dyv_31[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_32[0][j][i] = (state_old->h_32[0][j][i] + ((tend1->dxh_32[0][j][i] + tend1->dyh_32[0][j][i]) * dt));
      state_new->u_32[0][j][i] = (state_old->u_32[0][j][i] + ((tend1->dxu_32[0][j][i] + tend1->dyu_32[0][j][i]) * dt));
      state_new->v_32[0][j][i] = (state_old->v_32[0][j][i] + ((tend1->dxv_32[0][j][i] + tend1->dyv_32[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_33[0][j][i] = (state_old->h_33[0][j][i] + ((tend1->dxh_33[0][j][i] + tend1->dyh_33[0][j][i]) * dt));
      state_new->u_33[0][j][i] = (state_old->u_33[0][j][i] + ((tend1->dxu_33[0][j][i] + tend1->dyu_33[0][j][i]) * dt));
      state_new->v_33[0][j][i] = (state_old->v_33[0][j][i] + ((tend1->dxv_33[0][j][i] + tend1->dyv_33[0][j][i]) * dt));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_11[0][j][i] = (state_old->h_11[0][j][i] + (((((tend1->dxh_11[0][j][i] + tend1->dyh_11[0][j][i]) + tend2->dxh_11[0][j][i]) + tend2->dyh_11[0][j][i]) * dt) / 4.0));
      state_new->u_11[0][j][i] = (state_old->u_11[0][j][i] + (((((tend1->dxu_11[0][j][i] + tend1->dyu_11[0][j][i]) + tend2->dxu_11[0][j][i]) + tend2->dyu_11[0][j][i]) * dt) / 4.0));
      state_new->v_11[0][j][i] = (state_old->v_11[0][j][i] + (((((tend1->dxv_11[0][j][i] + tend1->dyv_11[0][j][i]) + tend2->dxv_11[0][j][i]) + tend2->dyv_11[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_12[0][j][i] = (state_old->h_12[0][j][i] + (((((tend1->dxh_12[0][j][i] + tend1->dyh_12[0][j][i]) + tend2->dxh_12[0][j][i]) + tend2->dyh_12[0][j][i]) * dt) / 4.0));
      state_new->u_12[0][j][i] = (state_old->u_12[0][j][i] + (((((tend1->dxu_12[0][j][i] + tend1->dyu_12[0][j][i]) + tend2->dxu_12[0][j][i]) + tend2->dyu_12[0][j][i]) * dt) / 4.0));
      state_new->v_12[0][j][i] = (state_old->v_12[0][j][i] + (((((tend1->dxv_12[0][j][i] + tend1->dyv_12[0][j][i]) + tend2->dxv_12[0][j][i]) + tend2->dyv_12[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_13[0][j][i] = (state_old->h_13[0][j][i] + (((((tend1->dxh_13[0][j][i] + tend1->dyh_13[0][j][i]) + tend2->dxh_13[0][j][i]) + tend2->dyh_13[0][j][i]) * dt) / 4.0));
      state_new->u_13[0][j][i] = (state_old->u_13[0][j][i] + (((((tend1->dxu_13[0][j][i] + tend1->dyu_13[0][j][i]) + tend2->dxu_13[0][j][i]) + tend2->dyu_13[0][j][i]) * dt) / 4.0));
      state_new->v_13[0][j][i] = (state_old->v_13[0][j][i] + (((((tend1->dxv_13[0][j][i] + tend1->dyv_13[0][j][i]) + tend2->dxv_13[0][j][i]) + tend2->dyv_13[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_21[0][j][i] = (state_old->h_21[0][j][i] + (((((tend1->dxh_21[0][j][i] + tend1->dyh_21[0][j][i]) + tend2->dxh_21[0][j][i]) + tend2->dyh_21[0][j][i]) * dt) / 4.0));
      state_new->u_21[0][j][i] = (state_old->u_21[0][j][i] + (((((tend1->dxu_21[0][j][i] + tend1->dyu_21[0][j][i]) + tend2->dxu_21[0][j][i]) + tend2->dyu_21[0][j][i]) * dt) / 4.0));
      state_new->v_21[0][j][i] = (state_old->v_21[0][j][i] + (((((tend1->dxv_21[0][j][i] + tend1->dyv_21[0][j][i]) + tend2->dxv_21[0][j][i]) + tend2->dyv_21[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_22[0][j][i] = (state_old->h_22[0][j][i] + (((((tend1->dxh_22[0][j][i] + tend1->dyh_22[0][j][i]) + tend2->dxh_22[0][j][i]) + tend2->dyh_22[0][j][i]) * dt) / 4.0));
      state_new->u_22[0][j][i] = (state_old->u_22[0][j][i] + (((((tend1->dxu_22[0][j][i] + tend1->dyu_22[0][j][i]) + tend2->dxu_22[0][j][i]) + tend2->dyu_22[0][j][i]) * dt) / 4.0));
      state_new->v_22[0][j][i] = (state_old->v_22[0][j][i] + (((((tend1->dxv_22[0][j][i] + tend1->dyv_22[0][j][i]) + tend2->dxv_22[0][j][i]) + tend2->dyv_22[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_23[0][j][i] = (state_old->h_23[0][j][i] + (((((tend1->dxh_23[0][j][i] + tend1->dyh_23[0][j][i]) + tend2->dxh_23[0][j][i]) + tend2->dyh_23[0][j][i]) * dt) / 4.0));
      state_new->u_23[0][j][i] = (state_old->u_23[0][j][i] + (((((tend1->dxu_23[0][j][i] + tend1->dyu_23[0][j][i]) + tend2->dxu_23[0][j][i]) + tend2->dyu_23[0][j][i]) * dt) / 4.0));
      state_new->v_23[0][j][i] = (state_old->v_23[0][j][i] + (((((tend1->dxv_23[0][j][i] + tend1->dyv_23[0][j][i]) + tend2->dxv_23[0][j][i]) + tend2->dyv_23[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_31[0][j][i] = (state_old->h_31[0][j][i] + (((((tend1->dxh_31[0][j][i] + tend1->dyh_31[0][j][i]) + tend2->dxh_31[0][j][i]) + tend2->dyh_31[0][j][i]) * dt) / 4.0));
      state_new->u_31[0][j][i] = (state_old->u_31[0][j][i] + (((((tend1->dxu_31[0][j][i] + tend1->dyu_31[0][j][i]) + tend2->dxu_31[0][j][i]) + tend2->dyu_31[0][j][i]) * dt) / 4.0));
      state_new->v_31[0][j][i] = (state_old->v_31[0][j][i] + (((((tend1->dxv_31[0][j][i] + tend1->dyv_31[0][j][i]) + tend2->dxv_31[0][j][i]) + tend2->dyv_31[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_32[0][j][i] = (state_old->h_32[0][j][i] + (((((tend1->dxh_32[0][j][i] + tend1->dyh_32[0][j][i]) + tend2->dxh_32[0][j][i]) + tend2->dyh_32[0][j][i]) * dt) / 4.0));
      state_new->u_32[0][j][i] = (state_old->u_32[0][j][i] + (((((tend1->dxu_32[0][j][i] + tend1->dyu_32[0][j][i]) + tend2->dxu_32[0][j][i]) + tend2->dyu_32[0][j][i]) * dt) / 4.0));
      state_new->v_32[0][j][i] = (state_old->v_32[0][j][i] + (((((tend1->dxv_32[0][j][i] + tend1->dyv_32[0][j][i]) + tend2->dxv_32[0][j][i]) + tend2->dyv_32[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_33[0][j][i] = (state_old->h_33[0][j][i] + (((((tend1->dxh_33[0][j][i] + tend1->dyh_33[0][j][i]) + tend2->dxh_33[0][j][i]) + tend2->dyh_33[0][j][i]) * dt) / 4.0));
      state_new->u_33[0][j][i] = (state_old->u_33[0][j][i] + (((((tend1->dxu_33[0][j][i] + tend1->dyu_33[0][j][i]) + tend2->dxu_33[0][j][i]) + tend2->dyu_33[0][j][i]) * dt) / 4.0));
      state_new->v_33[0][j][i] = (state_old->v_33[0][j][i] + (((((tend1->dxv_33[0][j][i] + tend1->dyv_33[0][j][i]) + tend2->dxv_33[0][j][i]) + tend2->dyv_33[0][j][i]) * dt) / 4.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_11[0][j][i] = (state_old->h_11[0][j][i] + (((((((tend1->dxh_11[0][j][i] + tend1->dyh_11[0][j][i]) + tend2->dxh_11[0][j][i]) + tend2->dyh_11[0][j][i]) + (4.0 * tend3->dxh_11[0][j][i])) + (4.0 * tend3->dyh_11[0][j][i])) * dt) / 6.0));
      state_new->u_11[0][j][i] = (state_old->u_11[0][j][i] + (((((((tend1->dxu_11[0][j][i] + tend1->dyu_11[0][j][i]) + tend2->dxu_11[0][j][i]) + tend2->dyu_11[0][j][i]) + (4.0 * tend3->dxu_11[0][j][i])) + (4.0 * tend3->dyu_11[0][j][i])) * dt) / 6.0));
      state_new->v_11[0][j][i] = (state_old->v_11[0][j][i] + (((((((tend1->dxv_11[0][j][i] + tend1->dyv_11[0][j][i]) + tend2->dxv_11[0][j][i]) + tend2->dyv_11[0][j][i]) + (4.0 * tend3->dxv_11[0][j][i])) + (4.0 * tend3->dyv_11[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_12[0][j][i] = (state_old->h_12[0][j][i] + (((((((tend1->dxh_12[0][j][i] + tend1->dyh_12[0][j][i]) + tend2->dxh_12[0][j][i]) + tend2->dyh_12[0][j][i]) + (4.0 * tend3->dxh_12[0][j][i])) + (4.0 * tend3->dyh_12[0][j][i])) * dt) / 6.0));
      state_new->u_12[0][j][i] = (state_old->u_12[0][j][i] + (((((((tend1->dxu_12[0][j][i] + tend1->dyu_12[0][j][i]) + tend2->dxu_12[0][j][i]) + tend2->dyu_12[0][j][i]) + (4.0 * tend3->dxu_12[0][j][i])) + (4.0 * tend3->dyu_12[0][j][i])) * dt) / 6.0));
      state_new->v_12[0][j][i] = (state_old->v_12[0][j][i] + (((((((tend1->dxv_12[0][j][i] + tend1->dyv_12[0][j][i]) + tend2->dxv_12[0][j][i]) + tend2->dyv_12[0][j][i]) + (4.0 * tend3->dxv_12[0][j][i])) + (4.0 * tend3->dyv_12[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_13[0][j][i] = (state_old->h_13[0][j][i] + (((((((tend1->dxh_13[0][j][i] + tend1->dyh_13[0][j][i]) + tend2->dxh_13[0][j][i]) + tend2->dyh_13[0][j][i]) + (4.0 * tend3->dxh_13[0][j][i])) + (4.0 * tend3->dyh_13[0][j][i])) * dt) / 6.0));
      state_new->u_13[0][j][i] = (state_old->u_13[0][j][i] + (((((((tend1->dxu_13[0][j][i] + tend1->dyu_13[0][j][i]) + tend2->dxu_13[0][j][i]) + tend2->dyu_13[0][j][i]) + (4.0 * tend3->dxu_13[0][j][i])) + (4.0 * tend3->dyu_13[0][j][i])) * dt) / 6.0));
      state_new->v_13[0][j][i] = (state_old->v_13[0][j][i] + (((((((tend1->dxv_13[0][j][i] + tend1->dyv_13[0][j][i]) + tend2->dxv_13[0][j][i]) + tend2->dyv_13[0][j][i]) + (4.0 * tend3->dxv_13[0][j][i])) + (4.0 * tend3->dyv_13[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_21[0][j][i] = (state_old->h_21[0][j][i] + (((((((tend1->dxh_21[0][j][i] + tend1->dyh_21[0][j][i]) + tend2->dxh_21[0][j][i]) + tend2->dyh_21[0][j][i]) + (4.0 * tend3->dxh_21[0][j][i])) + (4.0 * tend3->dyh_21[0][j][i])) * dt) / 6.0));
      state_new->u_21[0][j][i] = (state_old->u_21[0][j][i] + (((((((tend1->dxu_21[0][j][i] + tend1->dyu_21[0][j][i]) + tend2->dxu_21[0][j][i]) + tend2->dyu_21[0][j][i]) + (4.0 * tend3->dxu_21[0][j][i])) + (4.0 * tend3->dyu_21[0][j][i])) * dt) / 6.0));
      state_new->v_21[0][j][i] = (state_old->v_21[0][j][i] + (((((((tend1->dxv_21[0][j][i] + tend1->dyv_21[0][j][i]) + tend2->dxv_21[0][j][i]) + tend2->dyv_21[0][j][i]) + (4.0 * tend3->dxv_21[0][j][i])) + (4.0 * tend3->dyv_21[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_22[0][j][i] = (state_old->h_22[0][j][i] + (((((((tend1->dxh_22[0][j][i] + tend1->dyh_22[0][j][i]) + tend2->dxh_22[0][j][i]) + tend2->dyh_22[0][j][i]) + (4.0 * tend3->dxh_22[0][j][i])) + (4.0 * tend3->dyh_22[0][j][i])) * dt) / 6.0));
      state_new->u_22[0][j][i] = (state_old->u_22[0][j][i] + (((((((tend1->dxu_22[0][j][i] + tend1->dyu_22[0][j][i]) + tend2->dxu_22[0][j][i]) + tend2->dyu_22[0][j][i]) + (4.0 * tend3->dxu_22[0][j][i])) + (4.0 * tend3->dyu_22[0][j][i])) * dt) / 6.0));
      state_new->v_22[0][j][i] = (state_old->v_22[0][j][i] + (((((((tend1->dxv_22[0][j][i] + tend1->dyv_22[0][j][i]) + tend2->dxv_22[0][j][i]) + tend2->dyv_22[0][j][i]) + (4.0 * tend3->dxv_22[0][j][i])) + (4.0 * tend3->dyv_22[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_23[0][j][i] = (state_old->h_23[0][j][i] + (((((((tend1->dxh_23[0][j][i] + tend1->dyh_23[0][j][i]) + tend2->dxh_23[0][j][i]) + tend2->dyh_23[0][j][i]) + (4.0 * tend3->dxh_23[0][j][i])) + (4.0 * tend3->dyh_23[0][j][i])) * dt) / 6.0));
      state_new->u_23[0][j][i] = (state_old->u_23[0][j][i] + (((((((tend1->dxu_23[0][j][i] + tend1->dyu_23[0][j][i]) + tend2->dxu_23[0][j][i]) + tend2->dyu_23[0][j][i]) + (4.0 * tend3->dxu_23[0][j][i])) + (4.0 * tend3->dyu_23[0][j][i])) * dt) / 6.0));
      state_new->v_23[0][j][i] = (state_old->v_23[0][j][i] + (((((((tend1->dxv_23[0][j][i] + tend1->dyv_23[0][j][i]) + tend2->dxv_23[0][j][i]) + tend2->dyv_23[0][j][i]) + (4.0 * tend3->dxv_23[0][j][i])) + (4.0 * tend3->dyv_23[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_31[0][j][i] = (state_old->h_31[0][j][i] + (((((((tend1->dxh_31[0][j][i] + tend1->dyh_31[0][j][i]) + tend2->dxh_31[0][j][i]) + tend2->dyh_31[0][j][i]) + (4.0 * tend3->dxh_31[0][j][i])) + (4.0 * tend3->dyh_31[0][j][i])) * dt) / 6.0));
      state_new->u_31[0][j][i] = (state_old->u_31[0][j][i] + (((((((tend1->dxu_31[0][j][i] + tend1->dyu_31[0][j][i]) + tend2->dxu_31[0][j][i]) + tend2->dyu_31[0][j][i]) + (4.0 * tend3->dxu_31[0][j][i])) + (4.0 * tend3->dyu_31[0][j][i])) * dt) / 6.0));
      state_new->v_31[0][j][i] = (state_old->v_31[0][j][i] + (((((((tend1->dxv_31[0][j][i] + tend1->dyv_31[0][j][i]) + tend2->dxv_31[0][j][i]) + tend2->dyv_31[0][j][i]) + (4.0 * tend3->dxv_31[0][j][i])) + (4.0 * tend3->dyv_31[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_32[0][j][i] = (state_old->h_32[0][j][i] + (((((((tend1->dxh_32[0][j][i] + tend1->dyh_32[0][j][i]) + tend2->dxh_32[0][j][i]) + tend2->dyh_32[0][j][i]) + (4.0 * tend3->dxh_32[0][j][i])) + (4.0 * tend3->dyh_32[0][j][i])) * dt) / 6.0));
      state_new->u_32[0][j][i] = (state_old->u_32[0][j][i] + (((((((tend1->dxu_32[0][j][i] + tend1->dyu_32[0][j][i]) + tend2->dxu_32[0][j][i]) + tend2->dyu_32[0][j][i]) + (4.0 * tend3->dxu_32[0][j][i])) + (4.0 * tend3->dyu_32[0][j][i])) * dt) / 6.0));
      state_new->v_32[0][j][i] = (state_old->v_32[0][j][i] + (((((((tend1->dxv_32[0][j][i] + tend1->dyv_32[0][j][i]) + tend2->dxv_32[0][j][i]) + tend2->dyv_32[0][j][i]) + (4.0 * tend3->dxv_32[0][j][i])) + (4.0 * tend3->dyv_32[0][j][i])) * dt) / 6.0));
    }
  }
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
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 600)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 600)-Proc.lon_beg+Proc.lon_hw; i+=1){
      state_new->h_33[0][j][i] = (state_old->h_33[0][j][i] + (((((((tend1->dxh_33[0][j][i] + tend1->dyh_33[0][j][i]) + tend2->dxh_33[0][j][i]) + tend2->dyh_33[0][j][i]) + (4.0 * tend3->dxh_33[0][j][i])) + (4.0 * tend3->dyh_33[0][j][i])) * dt) / 6.0));
      state_new->u_33[0][j][i] = (state_old->u_33[0][j][i] + (((((((tend1->dxu_33[0][j][i] + tend1->dyu_33[0][j][i]) + tend2->dxu_33[0][j][i]) + tend2->dyu_33[0][j][i]) + (4.0 * tend3->dxu_33[0][j][i])) + (4.0 * tend3->dyu_33[0][j][i])) * dt) / 6.0));
      state_new->v_33[0][j][i] = (state_old->v_33[0][j][i] + (((((((tend1->dxv_33[0][j][i] + tend1->dyv_33[0][j][i]) + tend2->dxv_33[0][j][i]) + tend2->dyv_33[0][j][i]) + (4.0 * tend3->dxv_33[0][j][i])) + (4.0 * tend3->dyv_33[0][j][i])) * dt) / 6.0));
    }
  }
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
  mcvupdateXY_0(&state[2], &tend[Timeinfo.newt], &mesh[0]);
  update_rk2_0(&state[Timeinfo.oldt], &state[2], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], dt);
  mcvupdateXY_0(&state[2], &tend[2], &mesh[0]);
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