#include <stdio.h>
#include <stdlib.h>
#include <simd.h>
#include <math.h>
#include "slave.h"

#include "CPE.h"
#include "slave_struct.h"
#include "physical_variable.h"
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

__thread_local int fnx;
__thread_local int fny;
__thread_local int fnz;
__thread_local int hnx;
__thread_local int hny;
__thread_local int hnz;
__thread_local int ghx;
__thread_local int ghy;
__thread_local int ghz;

void global_prepare(void *_ptr){
  global_info *gi = (global_info*)(_ptr);
  fnx = gi->fnx;
  fny = gi->fny;
  fnz = gi->fnz;
  hnx = gi->hnx;
  hny = gi->hny;
  hnz = gi->hnz;
  ghx = gi->ghx;
  ghy = gi->ghy;
  ghz = gi->ghz;
}

static void CalcID(int *id, int *rid, int *cid, int *w, int *e, int *s, int *n, int mx, int my){
  *id = _MYID;
  *rid = ROW(_MYID);
  *cid = COL(_MYID);

  if ((*cid % mx) == 0) *w = -1;
  else *w = *id - 1;
  if ((*cid % mx) == (mx - 1)) *e = -1;
  else *e = *id + 1;
  if ((*rid % my) == 0) *s = -1;
  else *s = *id - 8;
  if ((*rid % my) == (my - 1)) *n = -1;
  else *n = *id + 8;
}

static void RoundRobin(int n, int p, int size, int *beg, int *end, int *len){
  int cur = 0;
  int i;
  int cursize;
  for (i = 0 ; i < p ; i++){
    if (size % (n - i) == 0) cursize = size / (n - i);
    else cursize = size / (n - i) + 1;
    size -= cursize;
    cur += cursize;
  }
  *beg = cur;
  if (size % (n - p) == 0) cursize = size / (n - p);
  else cursize = size / (n - p) + 1;
  *len = cursize;
  *end = *beg + cursize;
}

static void CalcRange(int id ,int rid, int cid, int nx,int ny, int nz, int mx, int my, int mz, \ 
                      int *xbeg, int *xend, int *xlen, int *ybeg, int *yend, int *ylen, int *zbeg, int *zend, int *zlen){
  int px,py,pz;
  pz = ((rid / my) * (8 / mx) + (cid / mx));
  RoundRobin(mz,pz,nz,zbeg,zend,zlen);
  py = rid % my;
  RoundRobin(my,py,ny,ybeg,yend,ylen);
  px = cid % mx;
  RoundRobin(mx,px,nx,xbeg,xend,xlen);
}

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












void updateState_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_0_info *data = (updateState_0_0_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_11 = data->state_oldh_11;
  __attribute__ ((aligned (32))) double state_oldh_11[6][130];
  
  double *_state_oldu_11 = data->state_oldu_11;
  __attribute__ ((aligned (32))) double state_oldu_11[6][130];
  
  double *_state_oldv_11 = data->state_oldv_11;
  __attribute__ ((aligned (32))) double state_oldv_11[6][130];
  
  double *_state_oldjab_11 = data->state_oldjab_11;
  __attribute__ ((aligned (32))) double state_oldjab_11[6][130];
  
  double *_state_oldjab5_11 = data->state_oldjab5_11;
  __attribute__ ((aligned (32))) double state_oldjab5_11[6][130];
  
  double *_state_oldjab6_11 = data->state_oldjab6_11;
  __attribute__ ((aligned (32))) double state_oldjab6_11[6][130];
  
  double *_state_oldjab7_11 = data->state_oldjab7_11;
  __attribute__ ((aligned (32))) double state_oldjab7_11[6][130];
  
  double *_state_olduc_11 = data->state_olduc_11;
  __attribute__ ((aligned (32))) double state_olduc_11[6][130];
  
  double *_state_oldvc_11 = data->state_oldvc_11;
  __attribute__ ((aligned (32))) double state_oldvc_11[6][130];
  
  double *_state_newh_11 = data->state_newh_11;
  __attribute__ ((aligned (32))) double state_newh_11[6][130];
  
  double *_state_newu_11 = data->state_newu_11;
  __attribute__ ((aligned (32))) double state_newu_11[6][130];
  
  double *_state_newv_11 = data->state_newv_11;
  __attribute__ ((aligned (32))) double state_newv_11[6][130];
  
  double *_state_newjab_11 = data->state_newjab_11;
  __attribute__ ((aligned (32))) double state_newjab_11[6][130];
  
  double *_state_newjab5_11 = data->state_newjab5_11;
  __attribute__ ((aligned (32))) double state_newjab5_11[6][130];
  
  double *_state_newjab6_11 = data->state_newjab6_11;
  __attribute__ ((aligned (32))) double state_newjab6_11[6][130];
  
  double *_state_newjab7_11 = data->state_newjab7_11;
  __attribute__ ((aligned (32))) double state_newjab7_11[6][130];
  
  double *_state_newuc_11 = data->state_newuc_11;
  __attribute__ ((aligned (32))) double state_newuc_11[6][130];
  
  double *_state_newvc_11 = data->state_newvc_11;
  __attribute__ ((aligned (32))) double state_newvc_11[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_11+(j*fnumx+ib), &state_oldjab_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_11+(j*fnumx+ib), &state_oldjab5_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_11+(j*fnumx+ib), &state_oldjab6_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_11+(j*fnumx+ib), &state_oldjab7_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_11+(j*fnumx+ib), &state_olduc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_11+(j*fnumx+ib), &state_oldvc_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_11+(j*fnumx+ib), &state_oldjab_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_11+(j*fnumx+ib), &state_oldjab5_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_11+(j*fnumx+ib), &state_oldjab6_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_11+(j*fnumx+ib), &state_oldjab7_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_11+(j*fnumx+ib), &state_olduc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_11+(j*fnumx+ib), &state_oldvc_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_11[j][i] = state_oldh_11[swindex[j]][i];
            state_newu_11[j][i] = state_oldu_11[swindex[j]][i];
            state_newv_11[j][i] = state_oldv_11[swindex[j]][i];
            state_newjab_11[j][i] = state_oldjab_11[swindex[j]][i];
            state_newjab5_11[j][i] = state_oldjab5_11[swindex[j]][i];
            state_newjab6_11[j][i] = state_oldjab6_11[swindex[j]][i];
            state_newjab7_11[j][i] = state_oldjab7_11[swindex[j]][i];
            state_newuc_11[j][i] = state_olduc_11[swindex[j]][i];
            state_newvc_11[j][i] = state_oldvc_11[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_1(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_1_info *data = (updateState_0_1_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_12 = data->state_oldh_12;
  __attribute__ ((aligned (32))) double state_oldh_12[6][130];
  
  double *_state_oldu_12 = data->state_oldu_12;
  __attribute__ ((aligned (32))) double state_oldu_12[6][130];
  
  double *_state_oldv_12 = data->state_oldv_12;
  __attribute__ ((aligned (32))) double state_oldv_12[6][130];
  
  double *_state_oldjab_12 = data->state_oldjab_12;
  __attribute__ ((aligned (32))) double state_oldjab_12[6][130];
  
  double *_state_oldjab5_12 = data->state_oldjab5_12;
  __attribute__ ((aligned (32))) double state_oldjab5_12[6][130];
  
  double *_state_oldjab6_12 = data->state_oldjab6_12;
  __attribute__ ((aligned (32))) double state_oldjab6_12[6][130];
  
  double *_state_oldjab7_12 = data->state_oldjab7_12;
  __attribute__ ((aligned (32))) double state_oldjab7_12[6][130];
  
  double *_state_olduc_12 = data->state_olduc_12;
  __attribute__ ((aligned (32))) double state_olduc_12[6][130];
  
  double *_state_oldvc_12 = data->state_oldvc_12;
  __attribute__ ((aligned (32))) double state_oldvc_12[6][130];
  
  double *_state_newh_12 = data->state_newh_12;
  __attribute__ ((aligned (32))) double state_newh_12[6][130];
  
  double *_state_newu_12 = data->state_newu_12;
  __attribute__ ((aligned (32))) double state_newu_12[6][130];
  
  double *_state_newv_12 = data->state_newv_12;
  __attribute__ ((aligned (32))) double state_newv_12[6][130];
  
  double *_state_newjab_12 = data->state_newjab_12;
  __attribute__ ((aligned (32))) double state_newjab_12[6][130];
  
  double *_state_newjab5_12 = data->state_newjab5_12;
  __attribute__ ((aligned (32))) double state_newjab5_12[6][130];
  
  double *_state_newjab6_12 = data->state_newjab6_12;
  __attribute__ ((aligned (32))) double state_newjab6_12[6][130];
  
  double *_state_newjab7_12 = data->state_newjab7_12;
  __attribute__ ((aligned (32))) double state_newjab7_12[6][130];
  
  double *_state_newuc_12 = data->state_newuc_12;
  __attribute__ ((aligned (32))) double state_newuc_12[6][130];
  
  double *_state_newvc_12 = data->state_newvc_12;
  __attribute__ ((aligned (32))) double state_newvc_12[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_12+(j*fnumx+ib), &state_oldjab_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_12+(j*fnumx+ib), &state_oldjab5_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_12+(j*fnumx+ib), &state_oldjab6_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_12+(j*fnumx+ib), &state_oldjab7_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_12+(j*fnumx+ib), &state_olduc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_12+(j*fnumx+ib), &state_oldvc_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_12+(j*fnumx+ib), &state_oldjab_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_12+(j*fnumx+ib), &state_oldjab5_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_12+(j*fnumx+ib), &state_oldjab6_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_12+(j*fnumx+ib), &state_oldjab7_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_12+(j*fnumx+ib), &state_olduc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_12+(j*fnumx+ib), &state_oldvc_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_12[j][i] = state_oldh_12[swindex[j]][i];
            state_newu_12[j][i] = state_oldu_12[swindex[j]][i];
            state_newv_12[j][i] = state_oldv_12[swindex[j]][i];
            state_newjab_12[j][i] = state_oldjab_12[swindex[j]][i];
            state_newjab5_12[j][i] = state_oldjab5_12[swindex[j]][i];
            state_newjab6_12[j][i] = state_oldjab6_12[swindex[j]][i];
            state_newjab7_12[j][i] = state_oldjab7_12[swindex[j]][i];
            state_newuc_12[j][i] = state_olduc_12[swindex[j]][i];
            state_newvc_12[j][i] = state_oldvc_12[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_2(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_2_info *data = (updateState_0_2_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_13 = data->state_oldh_13;
  __attribute__ ((aligned (32))) double state_oldh_13[6][130];
  
  double *_state_oldu_13 = data->state_oldu_13;
  __attribute__ ((aligned (32))) double state_oldu_13[6][130];
  
  double *_state_oldv_13 = data->state_oldv_13;
  __attribute__ ((aligned (32))) double state_oldv_13[6][130];
  
  double *_state_oldjab_13 = data->state_oldjab_13;
  __attribute__ ((aligned (32))) double state_oldjab_13[6][130];
  
  double *_state_oldjab5_13 = data->state_oldjab5_13;
  __attribute__ ((aligned (32))) double state_oldjab5_13[6][130];
  
  double *_state_oldjab6_13 = data->state_oldjab6_13;
  __attribute__ ((aligned (32))) double state_oldjab6_13[6][130];
  
  double *_state_oldjab7_13 = data->state_oldjab7_13;
  __attribute__ ((aligned (32))) double state_oldjab7_13[6][130];
  
  double *_state_olduc_13 = data->state_olduc_13;
  __attribute__ ((aligned (32))) double state_olduc_13[6][130];
  
  double *_state_oldvc_13 = data->state_oldvc_13;
  __attribute__ ((aligned (32))) double state_oldvc_13[6][130];
  
  double *_state_newh_13 = data->state_newh_13;
  __attribute__ ((aligned (32))) double state_newh_13[6][130];
  
  double *_state_newu_13 = data->state_newu_13;
  __attribute__ ((aligned (32))) double state_newu_13[6][130];
  
  double *_state_newv_13 = data->state_newv_13;
  __attribute__ ((aligned (32))) double state_newv_13[6][130];
  
  double *_state_newjab_13 = data->state_newjab_13;
  __attribute__ ((aligned (32))) double state_newjab_13[6][130];
  
  double *_state_newjab5_13 = data->state_newjab5_13;
  __attribute__ ((aligned (32))) double state_newjab5_13[6][130];
  
  double *_state_newjab6_13 = data->state_newjab6_13;
  __attribute__ ((aligned (32))) double state_newjab6_13[6][130];
  
  double *_state_newjab7_13 = data->state_newjab7_13;
  __attribute__ ((aligned (32))) double state_newjab7_13[6][130];
  
  double *_state_newuc_13 = data->state_newuc_13;
  __attribute__ ((aligned (32))) double state_newuc_13[6][130];
  
  double *_state_newvc_13 = data->state_newvc_13;
  __attribute__ ((aligned (32))) double state_newvc_13[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_13+(j*fnumx+ib), &state_oldjab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_13+(j*fnumx+ib), &state_oldjab5_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_13+(j*fnumx+ib), &state_oldjab6_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_13+(j*fnumx+ib), &state_oldjab7_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_13+(j*fnumx+ib), &state_olduc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_13+(j*fnumx+ib), &state_oldvc_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_13+(j*fnumx+ib), &state_oldjab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_13+(j*fnumx+ib), &state_oldjab5_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_13+(j*fnumx+ib), &state_oldjab6_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_13+(j*fnumx+ib), &state_oldjab7_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_13+(j*fnumx+ib), &state_olduc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_13+(j*fnumx+ib), &state_oldvc_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_13[j][i] = state_oldh_13[swindex[j]][i];
            state_newu_13[j][i] = state_oldu_13[swindex[j]][i];
            state_newv_13[j][i] = state_oldv_13[swindex[j]][i];
            state_newjab_13[j][i] = state_oldjab_13[swindex[j]][i];
            state_newjab5_13[j][i] = state_oldjab5_13[swindex[j]][i];
            state_newjab6_13[j][i] = state_oldjab6_13[swindex[j]][i];
            state_newjab7_13[j][i] = state_oldjab7_13[swindex[j]][i];
            state_newuc_13[j][i] = state_olduc_13[swindex[j]][i];
            state_newvc_13[j][i] = state_oldvc_13[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_3(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_3_info *data = (updateState_0_3_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_21 = data->state_oldh_21;
  __attribute__ ((aligned (32))) double state_oldh_21[6][130];
  
  double *_state_oldu_21 = data->state_oldu_21;
  __attribute__ ((aligned (32))) double state_oldu_21[6][130];
  
  double *_state_oldv_21 = data->state_oldv_21;
  __attribute__ ((aligned (32))) double state_oldv_21[6][130];
  
  double *_state_oldjab_21 = data->state_oldjab_21;
  __attribute__ ((aligned (32))) double state_oldjab_21[6][130];
  
  double *_state_oldjab5_21 = data->state_oldjab5_21;
  __attribute__ ((aligned (32))) double state_oldjab5_21[6][130];
  
  double *_state_oldjab6_21 = data->state_oldjab6_21;
  __attribute__ ((aligned (32))) double state_oldjab6_21[6][130];
  
  double *_state_oldjab7_21 = data->state_oldjab7_21;
  __attribute__ ((aligned (32))) double state_oldjab7_21[6][130];
  
  double *_state_olduc_21 = data->state_olduc_21;
  __attribute__ ((aligned (32))) double state_olduc_21[6][130];
  
  double *_state_oldvc_21 = data->state_oldvc_21;
  __attribute__ ((aligned (32))) double state_oldvc_21[6][130];
  
  double *_state_newh_21 = data->state_newh_21;
  __attribute__ ((aligned (32))) double state_newh_21[6][130];
  
  double *_state_newu_21 = data->state_newu_21;
  __attribute__ ((aligned (32))) double state_newu_21[6][130];
  
  double *_state_newv_21 = data->state_newv_21;
  __attribute__ ((aligned (32))) double state_newv_21[6][130];
  
  double *_state_newjab_21 = data->state_newjab_21;
  __attribute__ ((aligned (32))) double state_newjab_21[6][130];
  
  double *_state_newjab5_21 = data->state_newjab5_21;
  __attribute__ ((aligned (32))) double state_newjab5_21[6][130];
  
  double *_state_newjab6_21 = data->state_newjab6_21;
  __attribute__ ((aligned (32))) double state_newjab6_21[6][130];
  
  double *_state_newjab7_21 = data->state_newjab7_21;
  __attribute__ ((aligned (32))) double state_newjab7_21[6][130];
  
  double *_state_newuc_21 = data->state_newuc_21;
  __attribute__ ((aligned (32))) double state_newuc_21[6][130];
  
  double *_state_newvc_21 = data->state_newvc_21;
  __attribute__ ((aligned (32))) double state_newvc_21[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_21+(j*fnumx+ib), &state_oldjab_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_21+(j*fnumx+ib), &state_oldjab5_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_21+(j*fnumx+ib), &state_oldjab6_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_21+(j*fnumx+ib), &state_oldjab7_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_21+(j*fnumx+ib), &state_olduc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_21+(j*fnumx+ib), &state_oldvc_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_21+(j*fnumx+ib), &state_oldjab_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_21+(j*fnumx+ib), &state_oldjab5_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_21+(j*fnumx+ib), &state_oldjab6_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_21+(j*fnumx+ib), &state_oldjab7_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_21+(j*fnumx+ib), &state_olduc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_21+(j*fnumx+ib), &state_oldvc_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_21[j][i] = state_oldh_21[swindex[j]][i];
            state_newu_21[j][i] = state_oldu_21[swindex[j]][i];
            state_newv_21[j][i] = state_oldv_21[swindex[j]][i];
            state_newjab_21[j][i] = state_oldjab_21[swindex[j]][i];
            state_newjab5_21[j][i] = state_oldjab5_21[swindex[j]][i];
            state_newjab6_21[j][i] = state_oldjab6_21[swindex[j]][i];
            state_newjab7_21[j][i] = state_oldjab7_21[swindex[j]][i];
            state_newuc_21[j][i] = state_olduc_21[swindex[j]][i];
            state_newvc_21[j][i] = state_oldvc_21[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_4(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_4_info *data = (updateState_0_4_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_22 = data->state_oldh_22;
  __attribute__ ((aligned (32))) double state_oldh_22[6][130];
  
  double *_state_oldu_22 = data->state_oldu_22;
  __attribute__ ((aligned (32))) double state_oldu_22[6][130];
  
  double *_state_oldv_22 = data->state_oldv_22;
  __attribute__ ((aligned (32))) double state_oldv_22[6][130];
  
  double *_state_oldjab_22 = data->state_oldjab_22;
  __attribute__ ((aligned (32))) double state_oldjab_22[6][130];
  
  double *_state_oldjab5_22 = data->state_oldjab5_22;
  __attribute__ ((aligned (32))) double state_oldjab5_22[6][130];
  
  double *_state_oldjab6_22 = data->state_oldjab6_22;
  __attribute__ ((aligned (32))) double state_oldjab6_22[6][130];
  
  double *_state_oldjab7_22 = data->state_oldjab7_22;
  __attribute__ ((aligned (32))) double state_oldjab7_22[6][130];
  
  double *_state_olduc_22 = data->state_olduc_22;
  __attribute__ ((aligned (32))) double state_olduc_22[6][130];
  
  double *_state_oldvc_22 = data->state_oldvc_22;
  __attribute__ ((aligned (32))) double state_oldvc_22[6][130];
  
  double *_state_newh_22 = data->state_newh_22;
  __attribute__ ((aligned (32))) double state_newh_22[6][130];
  
  double *_state_newu_22 = data->state_newu_22;
  __attribute__ ((aligned (32))) double state_newu_22[6][130];
  
  double *_state_newv_22 = data->state_newv_22;
  __attribute__ ((aligned (32))) double state_newv_22[6][130];
  
  double *_state_newjab_22 = data->state_newjab_22;
  __attribute__ ((aligned (32))) double state_newjab_22[6][130];
  
  double *_state_newjab5_22 = data->state_newjab5_22;
  __attribute__ ((aligned (32))) double state_newjab5_22[6][130];
  
  double *_state_newjab6_22 = data->state_newjab6_22;
  __attribute__ ((aligned (32))) double state_newjab6_22[6][130];
  
  double *_state_newjab7_22 = data->state_newjab7_22;
  __attribute__ ((aligned (32))) double state_newjab7_22[6][130];
  
  double *_state_newuc_22 = data->state_newuc_22;
  __attribute__ ((aligned (32))) double state_newuc_22[6][130];
  
  double *_state_newvc_22 = data->state_newvc_22;
  __attribute__ ((aligned (32))) double state_newvc_22[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_22+(j*fnumx+ib), &state_oldjab_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_22+(j*fnumx+ib), &state_oldjab5_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_22+(j*fnumx+ib), &state_oldjab6_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_22+(j*fnumx+ib), &state_oldjab7_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_22+(j*fnumx+ib), &state_olduc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_22+(j*fnumx+ib), &state_oldvc_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_22+(j*fnumx+ib), &state_oldjab_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_22+(j*fnumx+ib), &state_oldjab5_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_22+(j*fnumx+ib), &state_oldjab6_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_22+(j*fnumx+ib), &state_oldjab7_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_22+(j*fnumx+ib), &state_olduc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_22+(j*fnumx+ib), &state_oldvc_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_22[j][i] = state_oldh_22[swindex[j]][i];
            state_newu_22[j][i] = state_oldu_22[swindex[j]][i];
            state_newv_22[j][i] = state_oldv_22[swindex[j]][i];
            state_newjab_22[j][i] = state_oldjab_22[swindex[j]][i];
            state_newjab5_22[j][i] = state_oldjab5_22[swindex[j]][i];
            state_newjab6_22[j][i] = state_oldjab6_22[swindex[j]][i];
            state_newjab7_22[j][i] = state_oldjab7_22[swindex[j]][i];
            state_newuc_22[j][i] = state_olduc_22[swindex[j]][i];
            state_newvc_22[j][i] = state_oldvc_22[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_5(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_5_info *data = (updateState_0_5_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_23 = data->state_oldh_23;
  __attribute__ ((aligned (32))) double state_oldh_23[6][130];
  
  double *_state_oldu_23 = data->state_oldu_23;
  __attribute__ ((aligned (32))) double state_oldu_23[6][130];
  
  double *_state_oldv_23 = data->state_oldv_23;
  __attribute__ ((aligned (32))) double state_oldv_23[6][130];
  
  double *_state_oldjab_23 = data->state_oldjab_23;
  __attribute__ ((aligned (32))) double state_oldjab_23[6][130];
  
  double *_state_oldjab5_23 = data->state_oldjab5_23;
  __attribute__ ((aligned (32))) double state_oldjab5_23[6][130];
  
  double *_state_oldjab6_23 = data->state_oldjab6_23;
  __attribute__ ((aligned (32))) double state_oldjab6_23[6][130];
  
  double *_state_oldjab7_23 = data->state_oldjab7_23;
  __attribute__ ((aligned (32))) double state_oldjab7_23[6][130];
  
  double *_state_olduc_23 = data->state_olduc_23;
  __attribute__ ((aligned (32))) double state_olduc_23[6][130];
  
  double *_state_oldvc_23 = data->state_oldvc_23;
  __attribute__ ((aligned (32))) double state_oldvc_23[6][130];
  
  double *_state_newh_23 = data->state_newh_23;
  __attribute__ ((aligned (32))) double state_newh_23[6][130];
  
  double *_state_newu_23 = data->state_newu_23;
  __attribute__ ((aligned (32))) double state_newu_23[6][130];
  
  double *_state_newv_23 = data->state_newv_23;
  __attribute__ ((aligned (32))) double state_newv_23[6][130];
  
  double *_state_newjab_23 = data->state_newjab_23;
  __attribute__ ((aligned (32))) double state_newjab_23[6][130];
  
  double *_state_newjab5_23 = data->state_newjab5_23;
  __attribute__ ((aligned (32))) double state_newjab5_23[6][130];
  
  double *_state_newjab6_23 = data->state_newjab6_23;
  __attribute__ ((aligned (32))) double state_newjab6_23[6][130];
  
  double *_state_newjab7_23 = data->state_newjab7_23;
  __attribute__ ((aligned (32))) double state_newjab7_23[6][130];
  
  double *_state_newuc_23 = data->state_newuc_23;
  __attribute__ ((aligned (32))) double state_newuc_23[6][130];
  
  double *_state_newvc_23 = data->state_newvc_23;
  __attribute__ ((aligned (32))) double state_newvc_23[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_23+(j*fnumx+ib), &state_oldjab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_23+(j*fnumx+ib), &state_oldjab5_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_23+(j*fnumx+ib), &state_oldjab6_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_23+(j*fnumx+ib), &state_oldjab7_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_23+(j*fnumx+ib), &state_olduc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_23+(j*fnumx+ib), &state_oldvc_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_23+(j*fnumx+ib), &state_oldjab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_23+(j*fnumx+ib), &state_oldjab5_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_23+(j*fnumx+ib), &state_oldjab6_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_23+(j*fnumx+ib), &state_oldjab7_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_23+(j*fnumx+ib), &state_olduc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_23+(j*fnumx+ib), &state_oldvc_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_23[j][i] = state_oldh_23[swindex[j]][i];
            state_newu_23[j][i] = state_oldu_23[swindex[j]][i];
            state_newv_23[j][i] = state_oldv_23[swindex[j]][i];
            state_newjab_23[j][i] = state_oldjab_23[swindex[j]][i];
            state_newjab5_23[j][i] = state_oldjab5_23[swindex[j]][i];
            state_newjab6_23[j][i] = state_oldjab6_23[swindex[j]][i];
            state_newjab7_23[j][i] = state_oldjab7_23[swindex[j]][i];
            state_newuc_23[j][i] = state_olduc_23[swindex[j]][i];
            state_newvc_23[j][i] = state_oldvc_23[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_6(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_6_info *data = (updateState_0_6_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_31 = data->state_oldh_31;
  __attribute__ ((aligned (32))) double state_oldh_31[6][130];
  
  double *_state_oldu_31 = data->state_oldu_31;
  __attribute__ ((aligned (32))) double state_oldu_31[6][130];
  
  double *_state_oldv_31 = data->state_oldv_31;
  __attribute__ ((aligned (32))) double state_oldv_31[6][130];
  
  double *_state_oldjab_31 = data->state_oldjab_31;
  __attribute__ ((aligned (32))) double state_oldjab_31[6][130];
  
  double *_state_oldjab5_31 = data->state_oldjab5_31;
  __attribute__ ((aligned (32))) double state_oldjab5_31[6][130];
  
  double *_state_oldjab6_31 = data->state_oldjab6_31;
  __attribute__ ((aligned (32))) double state_oldjab6_31[6][130];
  
  double *_state_oldjab7_31 = data->state_oldjab7_31;
  __attribute__ ((aligned (32))) double state_oldjab7_31[6][130];
  
  double *_state_olduc_31 = data->state_olduc_31;
  __attribute__ ((aligned (32))) double state_olduc_31[6][130];
  
  double *_state_oldvc_31 = data->state_oldvc_31;
  __attribute__ ((aligned (32))) double state_oldvc_31[6][130];
  
  double *_state_newh_31 = data->state_newh_31;
  __attribute__ ((aligned (32))) double state_newh_31[6][130];
  
  double *_state_newu_31 = data->state_newu_31;
  __attribute__ ((aligned (32))) double state_newu_31[6][130];
  
  double *_state_newv_31 = data->state_newv_31;
  __attribute__ ((aligned (32))) double state_newv_31[6][130];
  
  double *_state_newjab_31 = data->state_newjab_31;
  __attribute__ ((aligned (32))) double state_newjab_31[6][130];
  
  double *_state_newjab5_31 = data->state_newjab5_31;
  __attribute__ ((aligned (32))) double state_newjab5_31[6][130];
  
  double *_state_newjab6_31 = data->state_newjab6_31;
  __attribute__ ((aligned (32))) double state_newjab6_31[6][130];
  
  double *_state_newjab7_31 = data->state_newjab7_31;
  __attribute__ ((aligned (32))) double state_newjab7_31[6][130];
  
  double *_state_newuc_31 = data->state_newuc_31;
  __attribute__ ((aligned (32))) double state_newuc_31[6][130];
  
  double *_state_newvc_31 = data->state_newvc_31;
  __attribute__ ((aligned (32))) double state_newvc_31[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_31+(j*fnumx+ib), &state_oldjab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_31+(j*fnumx+ib), &state_oldjab5_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_31+(j*fnumx+ib), &state_oldjab6_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_31+(j*fnumx+ib), &state_oldjab7_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_31+(j*fnumx+ib), &state_olduc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_31+(j*fnumx+ib), &state_oldvc_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_31+(j*fnumx+ib), &state_oldjab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_31+(j*fnumx+ib), &state_oldjab5_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_31+(j*fnumx+ib), &state_oldjab6_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_31+(j*fnumx+ib), &state_oldjab7_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_31+(j*fnumx+ib), &state_olduc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_31+(j*fnumx+ib), &state_oldvc_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_31[j][i] = state_oldh_31[swindex[j]][i];
            state_newu_31[j][i] = state_oldu_31[swindex[j]][i];
            state_newv_31[j][i] = state_oldv_31[swindex[j]][i];
            state_newjab_31[j][i] = state_oldjab_31[swindex[j]][i];
            state_newjab5_31[j][i] = state_oldjab5_31[swindex[j]][i];
            state_newjab6_31[j][i] = state_oldjab6_31[swindex[j]][i];
            state_newjab7_31[j][i] = state_oldjab7_31[swindex[j]][i];
            state_newuc_31[j][i] = state_olduc_31[swindex[j]][i];
            state_newvc_31[j][i] = state_oldvc_31[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_7(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_7_info *data = (updateState_0_7_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_32 = data->state_oldh_32;
  __attribute__ ((aligned (32))) double state_oldh_32[6][130];
  
  double *_state_oldu_32 = data->state_oldu_32;
  __attribute__ ((aligned (32))) double state_oldu_32[6][130];
  
  double *_state_oldv_32 = data->state_oldv_32;
  __attribute__ ((aligned (32))) double state_oldv_32[6][130];
  
  double *_state_oldjab_32 = data->state_oldjab_32;
  __attribute__ ((aligned (32))) double state_oldjab_32[6][130];
  
  double *_state_oldjab5_32 = data->state_oldjab5_32;
  __attribute__ ((aligned (32))) double state_oldjab5_32[6][130];
  
  double *_state_oldjab6_32 = data->state_oldjab6_32;
  __attribute__ ((aligned (32))) double state_oldjab6_32[6][130];
  
  double *_state_oldjab7_32 = data->state_oldjab7_32;
  __attribute__ ((aligned (32))) double state_oldjab7_32[6][130];
  
  double *_state_olduc_32 = data->state_olduc_32;
  __attribute__ ((aligned (32))) double state_olduc_32[6][130];
  
  double *_state_oldvc_32 = data->state_oldvc_32;
  __attribute__ ((aligned (32))) double state_oldvc_32[6][130];
  
  double *_state_newh_32 = data->state_newh_32;
  __attribute__ ((aligned (32))) double state_newh_32[6][130];
  
  double *_state_newu_32 = data->state_newu_32;
  __attribute__ ((aligned (32))) double state_newu_32[6][130];
  
  double *_state_newv_32 = data->state_newv_32;
  __attribute__ ((aligned (32))) double state_newv_32[6][130];
  
  double *_state_newjab_32 = data->state_newjab_32;
  __attribute__ ((aligned (32))) double state_newjab_32[6][130];
  
  double *_state_newjab5_32 = data->state_newjab5_32;
  __attribute__ ((aligned (32))) double state_newjab5_32[6][130];
  
  double *_state_newjab6_32 = data->state_newjab6_32;
  __attribute__ ((aligned (32))) double state_newjab6_32[6][130];
  
  double *_state_newjab7_32 = data->state_newjab7_32;
  __attribute__ ((aligned (32))) double state_newjab7_32[6][130];
  
  double *_state_newuc_32 = data->state_newuc_32;
  __attribute__ ((aligned (32))) double state_newuc_32[6][130];
  
  double *_state_newvc_32 = data->state_newvc_32;
  __attribute__ ((aligned (32))) double state_newvc_32[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_32+(j*fnumx+ib), &state_oldjab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_32+(j*fnumx+ib), &state_oldjab5_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_32+(j*fnumx+ib), &state_oldjab6_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_32+(j*fnumx+ib), &state_oldjab7_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_32+(j*fnumx+ib), &state_olduc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_32+(j*fnumx+ib), &state_oldvc_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_32+(j*fnumx+ib), &state_oldjab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_32+(j*fnumx+ib), &state_oldjab5_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_32+(j*fnumx+ib), &state_oldjab6_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_32+(j*fnumx+ib), &state_oldjab7_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_32+(j*fnumx+ib), &state_olduc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_32+(j*fnumx+ib), &state_oldvc_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_32[j][i] = state_oldh_32[swindex[j]][i];
            state_newu_32[j][i] = state_oldu_32[swindex[j]][i];
            state_newv_32[j][i] = state_oldv_32[swindex[j]][i];
            state_newjab_32[j][i] = state_oldjab_32[swindex[j]][i];
            state_newjab5_32[j][i] = state_oldjab5_32[swindex[j]][i];
            state_newjab6_32[j][i] = state_oldjab6_32[swindex[j]][i];
            state_newjab7_32[j][i] = state_oldjab7_32[swindex[j]][i];
            state_newuc_32[j][i] = state_olduc_32[swindex[j]][i];
            state_newvc_32[j][i] = state_oldvc_32[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateState_0_8(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateState_0_8_info *data = (updateState_0_8_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_33 = data->state_oldh_33;
  __attribute__ ((aligned (32))) double state_oldh_33[6][130];
  
  double *_state_oldu_33 = data->state_oldu_33;
  __attribute__ ((aligned (32))) double state_oldu_33[6][130];
  
  double *_state_oldv_33 = data->state_oldv_33;
  __attribute__ ((aligned (32))) double state_oldv_33[6][130];
  
  double *_state_oldjab_33 = data->state_oldjab_33;
  __attribute__ ((aligned (32))) double state_oldjab_33[6][130];
  
  double *_state_oldjab5_33 = data->state_oldjab5_33;
  __attribute__ ((aligned (32))) double state_oldjab5_33[6][130];
  
  double *_state_oldjab6_33 = data->state_oldjab6_33;
  __attribute__ ((aligned (32))) double state_oldjab6_33[6][130];
  
  double *_state_oldjab7_33 = data->state_oldjab7_33;
  __attribute__ ((aligned (32))) double state_oldjab7_33[6][130];
  
  double *_state_olduc_33 = data->state_olduc_33;
  __attribute__ ((aligned (32))) double state_olduc_33[6][130];
  
  double *_state_oldvc_33 = data->state_oldvc_33;
  __attribute__ ((aligned (32))) double state_oldvc_33[6][130];
  
  double *_state_newh_33 = data->state_newh_33;
  __attribute__ ((aligned (32))) double state_newh_33[6][130];
  
  double *_state_newu_33 = data->state_newu_33;
  __attribute__ ((aligned (32))) double state_newu_33[6][130];
  
  double *_state_newv_33 = data->state_newv_33;
  __attribute__ ((aligned (32))) double state_newv_33[6][130];
  
  double *_state_newjab_33 = data->state_newjab_33;
  __attribute__ ((aligned (32))) double state_newjab_33[6][130];
  
  double *_state_newjab5_33 = data->state_newjab5_33;
  __attribute__ ((aligned (32))) double state_newjab5_33[6][130];
  
  double *_state_newjab6_33 = data->state_newjab6_33;
  __attribute__ ((aligned (32))) double state_newjab6_33[6][130];
  
  double *_state_newjab7_33 = data->state_newjab7_33;
  __attribute__ ((aligned (32))) double state_newjab7_33[6][130];
  
  double *_state_newuc_33 = data->state_newuc_33;
  __attribute__ ((aligned (32))) double state_newuc_33[6][130];
  
  double *_state_newvc_33 = data->state_newvc_33;
  __attribute__ ((aligned (32))) double state_newvc_33[6][130];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_33+(j*fnumx+ib), &state_oldjab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_33+(j*fnumx+ib), &state_oldjab5_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_33+(j*fnumx+ib), &state_oldjab6_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_33+(j*fnumx+ib), &state_oldjab7_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_33+(j*fnumx+ib), &state_olduc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_33+(j*fnumx+ib), &state_oldvc_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab_33+(j*fnumx+ib), &state_oldjab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab5_33+(j*fnumx+ib), &state_oldjab5_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab6_33+(j*fnumx+ib), &state_oldjab6_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldjab7_33+(j*fnumx+ib), &state_oldjab7_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_olduc_33+(j*fnumx+ib), &state_olduc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldvc_33+(j*fnumx+ib), &state_oldvc_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_33[j][i] = state_oldh_33[swindex[j]][i];
            state_newu_33[j][i] = state_oldu_33[swindex[j]][i];
            state_newv_33[j][i] = state_oldv_33[swindex[j]][i];
            state_newjab_33[j][i] = state_oldjab_33[swindex[j]][i];
            state_newjab5_33[j][i] = state_oldjab5_33[swindex[j]][i];
            state_newjab6_33[j][i] = state_oldjab6_33[swindex[j]][i];
            state_newjab7_33[j][i] = state_oldjab7_33[swindex[j]][i];
            state_newuc_33[j][i] = state_olduc_33[swindex[j]][i];
            state_newvc_33[j][i] = state_oldvc_33[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab5_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab5_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab6_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab6_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newjab7_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newjab7_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newuc_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newuc_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newvc_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newvc_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

void updateX_0_30(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_30_info *data = (updateX_0_30_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_statejab5_31 = data->statejab5_31;
  __attribute__ ((aligned (32))) double statejab5_31[6][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[6][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[6][130];
  
  double *_stateu_31 = data->stateu_31;
  __attribute__ ((aligned (32))) double stateu_31[6][130];
  
  double *_statev_31 = data->statev_31;
  __attribute__ ((aligned (32))) double statev_31[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_statejab_11 = data->statejab_11;
  __attribute__ ((aligned (32))) double statejab_11[6][130];
  
  double *_stateu_11 = data->stateu_11;
  __attribute__ ((aligned (32))) double stateu_11[6][130];
  
  double *_statev_11 = data->statev_11;
  __attribute__ ((aligned (32))) double statev_11[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double ulocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_31[swindex[j]][(i - 1)]) + sqrt((((statejab5_31[swindex[j]][(i - 1)] * 9.80616) * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)])));
            tendfh_1[j][i] = ((((stateuc_31[swindex[j]][(i - 1)] * stateh_31[swindex[j]][(i - 1)]) + (stateuc_11[swindex[j]][i] * stateh_11[swindex[j]][i])) + (ulocal * (stateh_31[swindex[j]][(i - 1)] - stateh_11[swindex[j]][i]))) / 2);
            tendfu_1[j][i] = ((((((9.80616 * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)]) + (0.5 * ((stateu_31[swindex[j]][(i - 1)] * stateuc_31[swindex[j]][(i - 1)]) + (statev_31[swindex[j]][(i - 1)] * statevc_31[swindex[j]][(i - 1)])))) + (((9.80616 * stateh_11[swindex[j]][i]) / statejab_11[swindex[j]][i]) + (0.5 * ((stateu_11[swindex[j]][i] * stateuc_11[swindex[j]][i]) + (statev_11[swindex[j]][i] * statevc_11[swindex[j]][i]))))) + (ulocal * (stateu_31[swindex[j]][(i - 1)] - stateu_11[swindex[j]][i]))) / 2);
            tendfv_1[j][i] = (((0.0 + 0.0) + (ulocal * (statev_31[swindex[j]][(i - 1)] - statev_11[swindex[j]][i]))) / 2);
            tendqh_1[j][i] = ((stateh_31[swindex[j]][(i - 1)] + stateh_11[swindex[j]][i]) / 2);
            tendqu_1[j][i] = ((stateu_31[swindex[j]][(i - 1)] + stateu_11[swindex[j]][i]) / 2);
            tendqv_1[j][i] = ((statev_31[swindex[j]][(i - 1)] + statev_11[swindex[j]][i]) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_31(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_31_info *data = (updateX_0_31_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[10][130];
  
  double *_statejab5_31 = data->statejab5_31;
  __attribute__ ((aligned (32))) double statejab5_31[10][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[10][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[10][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[10][130];
  
  double *_stateh_21 = data->stateh_21;
  __attribute__ ((aligned (32))) double stateh_21[10][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[10][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_31[swindex[j]][(i - 1)]) + sqrt((((statejab5_31[swindex[j]][(i - 1)] * 9.80616) * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)])));
            dql = ((stateh_11[swindex[j]][(i - 1)] - (4 * stateh_21[swindex[j]][(i - 1)])) + (3 * stateh_31[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateh_11[swindex[j]][i]) - (4 * stateh_21[swindex[j]][i])) + stateh_31[swindex[j]][i]);
            dfl = (((stateuc_11[swindex[j]][(i - 1)] * stateh_11[swindex[j]][(i - 1)]) - (4 * (stateuc_21[swindex[j]][(i - 1)] * stateh_21[swindex[j]][(i - 1)]))) + (3 * (stateuc_31[swindex[j]][(i - 1)] * stateh_31[swindex[j]][(i - 1)])));
            dfr = -(((3 * (stateuc_11[swindex[j]][i] * stateh_11[swindex[j]][i])) - (4 * (stateuc_21[swindex[j]][i] * stateh_21[swindex[j]][i]))) + (stateuc_31[swindex[j]][i] * stateh_31[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_32(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_32_info *data = (updateX_0_32_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_statejab5_31 = data->statejab5_31;
  __attribute__ ((aligned (32))) double statejab5_31[6][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[6][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[6][130];
  
  double *_stateu_11 = data->stateu_11;
  __attribute__ ((aligned (32))) double stateu_11[6][130];
  
  double *_stateu_21 = data->stateu_21;
  __attribute__ ((aligned (32))) double stateu_21[6][130];
  
  double *_stateu_31 = data->stateu_31;
  __attribute__ ((aligned (32))) double stateu_31[6][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[6][130];
  
  double *_statejab_11 = data->statejab_11;
  __attribute__ ((aligned (32))) double statejab_11[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_statev_11 = data->statev_11;
  __attribute__ ((aligned (32))) double statev_11[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_stateh_21 = data->stateh_21;
  __attribute__ ((aligned (32))) double stateh_21[6][130];
  
  double *_statejab_21 = data->statejab_21;
  __attribute__ ((aligned (32))) double statejab_21[6][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[6][130];
  
  double *_statev_21 = data->statev_21;
  __attribute__ ((aligned (32))) double statev_21[6][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[6][130];
  
  double *_statev_31 = data->statev_31;
  __attribute__ ((aligned (32))) double statev_31[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_31[swindex[j]][(i - 1)]) + sqrt((((statejab5_31[swindex[j]][(i - 1)] * 9.80616) * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)])));
            dql = ((stateu_11[swindex[j]][(i - 1)] - (4 * stateu_21[swindex[j]][(i - 1)])) + (3 * stateu_31[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateu_11[swindex[j]][i]) - (4 * stateu_21[swindex[j]][i])) + stateu_31[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_11[swindex[j]][(i - 1)]) / statejab_11[swindex[j]][(i - 1)]) + (0.5 * ((stateu_11[swindex[j]][(i - 1)] * stateuc_11[swindex[j]][(i - 1)]) + (statev_11[swindex[j]][(i - 1)] * statevc_11[swindex[j]][(i - 1)])))) - (4 * (((9.80616 * stateh_21[swindex[j]][(i - 1)]) / statejab_21[swindex[j]][(i - 1)]) + (0.5 * ((stateu_21[swindex[j]][(i - 1)] * stateuc_21[swindex[j]][(i - 1)]) + (statev_21[swindex[j]][(i - 1)] * statevc_21[swindex[j]][(i - 1)])))))) + (3 * (((9.80616 * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)]) + (0.5 * ((stateu_31[swindex[j]][(i - 1)] * stateuc_31[swindex[j]][(i - 1)]) + (statev_31[swindex[j]][(i - 1)] * statevc_31[swindex[j]][(i - 1)]))))));
            dfr = -(((3 * (((9.80616 * stateh_11[swindex[j]][i]) / statejab_11[swindex[j]][i]) + (0.5 * ((stateu_11[swindex[j]][i] * stateuc_11[swindex[j]][i]) + (statev_11[swindex[j]][i] * statevc_11[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_21[swindex[j]][i]) / statejab_21[swindex[j]][i]) + (0.5 * ((stateu_21[swindex[j]][i] * stateuc_21[swindex[j]][i]) + (statev_21[swindex[j]][i] * statevc_21[swindex[j]][i])))))) + (((9.80616 * stateh_31[swindex[j]][i]) / statejab_31[swindex[j]][i]) + (0.5 * ((stateu_31[swindex[j]][i] * stateuc_31[swindex[j]][i]) + (statev_31[swindex[j]][i] * statevc_31[swindex[j]][i])))));
            tendfu_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_33(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_33_info *data = (updateX_0_33_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[10][130];
  
  double *_statejab5_31 = data->statejab5_31;
  __attribute__ ((aligned (32))) double statejab5_31[10][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[10][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[10][130];
  
  double *_statev_11 = data->statev_11;
  __attribute__ ((aligned (32))) double statev_11[10][130];
  
  double *_statev_21 = data->statev_21;
  __attribute__ ((aligned (32))) double statev_21[10][130];
  
  double *_statev_31 = data->statev_31;
  __attribute__ ((aligned (32))) double statev_31[10][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[10][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_31+(j*fnumx+ib), &statejab5_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_31[swindex[j]][(i - 1)]) + sqrt((((statejab5_31[swindex[j]][(i - 1)] * 9.80616) * stateh_31[swindex[j]][(i - 1)]) / statejab_31[swindex[j]][(i - 1)])));
            dql = ((statev_11[swindex[j]][(i - 1)] - (4 * statev_21[swindex[j]][(i - 1)])) + (3 * statev_31[swindex[j]][(i - 1)]));
            dqr = -(((3 * statev_11[swindex[j]][i]) - (4 * statev_21[swindex[j]][i])) + statev_31[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfv_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_34(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_34_info *data = (updateX_0_34_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[6][130];
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_tenddxh_11 = data->tenddxh_11;
  __attribute__ ((aligned (32))) double tenddxh_11[6][130];
  
  double *_tenddxh_21 = data->tenddxh_21;
  __attribute__ ((aligned (32))) double tenddxh_21[6][130];
  
  double *_tenddxh_31 = data->tenddxh_31;
  __attribute__ ((aligned (32))) double tenddxh_31[6][130];
  
  double *_tenddxu_11 = data->tenddxu_11;
  __attribute__ ((aligned (32))) double tenddxu_11[6][130];
  
  double *_tenddxu_21 = data->tenddxu_21;
  __attribute__ ((aligned (32))) double tenddxu_21[6][130];
  
  double *_tenddxu_31 = data->tenddxu_31;
  __attribute__ ((aligned (32))) double tenddxu_31[6][130];
  
  double *_tenddxv_11 = data->tenddxv_11;
  __attribute__ ((aligned (32))) double tenddxv_11[6][130];
  
  double *_tenddxv_21 = data->tenddxv_21;
  __attribute__ ((aligned (32))) double tenddxv_21[6][130];
  
  double *_tenddxv_31 = data->tenddxv_31;
  __attribute__ ((aligned (32))) double tenddxv_31[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqv_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966);
            sc = (((((3.0 * (tendqv_1[swindex[j]][(i + 1)] - tendqv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_11[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddxh_21[j][i] = (((((-3.0 * (tendfh_1[swindex[j]][(i + 1)] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_31[j][i] = (-tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966);
            tenddxu_11[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) + (statevc_11[swindex[j]][i] * sl));
            tenddxu_21[j][i] = ((((((3.0 * (tendfu_1[swindex[j]][(i + 1)] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) + (statevc_21[swindex[j]][i] * sc));
            tenddxu_31[j][i] = ((-tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) + (statevc_31[swindex[j]][i] * sr));
            tenddxv_11[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) - (stateuc_11[swindex[j]][i] * sl));
            tenddxv_21[j][i] = ((((((-3.0 * (tendfv_1[swindex[j]][(i + 1)] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) - (stateuc_21[swindex[j]][i] * sc));
            tenddxv_31[j][i] = ((-tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) - (stateuc_31[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddxh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_35(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_35_info *data = (updateX_0_35_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[6][130];
  
  double *_statejab5_32 = data->statejab5_32;
  __attribute__ ((aligned (32))) double statejab5_32[6][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[6][130];
  
  double *_statejab_32 = data->statejab_32;
  __attribute__ ((aligned (32))) double statejab_32[6][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[6][130];
  
  double *_stateh_12 = data->stateh_12;
  __attribute__ ((aligned (32))) double stateh_12[6][130];
  
  double *_stateu_32 = data->stateu_32;
  __attribute__ ((aligned (32))) double stateu_32[6][130];
  
  double *_statev_32 = data->statev_32;
  __attribute__ ((aligned (32))) double statev_32[6][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[6][130];
  
  double *_statejab_12 = data->statejab_12;
  __attribute__ ((aligned (32))) double statejab_12[6][130];
  
  double *_stateu_12 = data->stateu_12;
  __attribute__ ((aligned (32))) double stateu_12[6][130];
  
  double *_statev_12 = data->statev_12;
  __attribute__ ((aligned (32))) double statev_12[6][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double ulocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_32[swindex[j]][(i - 1)]) + sqrt((((statejab5_32[swindex[j]][(i - 1)] * 9.80616) * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)])));
            tendfh_1[j][i] = ((((stateuc_32[swindex[j]][(i - 1)] * stateh_32[swindex[j]][(i - 1)]) + (stateuc_12[swindex[j]][i] * stateh_12[swindex[j]][i])) + (ulocal * (stateh_32[swindex[j]][(i - 1)] - stateh_12[swindex[j]][i]))) / 2);
            tendfu_1[j][i] = ((((((9.80616 * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)]) + (0.5 * ((stateu_32[swindex[j]][(i - 1)] * stateuc_32[swindex[j]][(i - 1)]) + (statev_32[swindex[j]][(i - 1)] * statevc_32[swindex[j]][(i - 1)])))) + (((9.80616 * stateh_12[swindex[j]][i]) / statejab_12[swindex[j]][i]) + (0.5 * ((stateu_12[swindex[j]][i] * stateuc_12[swindex[j]][i]) + (statev_12[swindex[j]][i] * statevc_12[swindex[j]][i]))))) + (ulocal * (stateu_32[swindex[j]][(i - 1)] - stateu_12[swindex[j]][i]))) / 2);
            tendfv_1[j][i] = (((0.0 + 0.0) + (ulocal * (statev_32[swindex[j]][(i - 1)] - statev_12[swindex[j]][i]))) / 2);
            tendqh_1[j][i] = ((stateh_32[swindex[j]][(i - 1)] + stateh_12[swindex[j]][i]) / 2);
            tendqu_1[j][i] = ((stateu_32[swindex[j]][(i - 1)] + stateu_12[swindex[j]][i]) / 2);
            tendqv_1[j][i] = ((statev_32[swindex[j]][(i - 1)] + statev_12[swindex[j]][i]) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_36(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_36_info *data = (updateX_0_36_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[10][130];
  
  double *_statejab5_32 = data->statejab5_32;
  __attribute__ ((aligned (32))) double statejab5_32[10][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[10][130];
  
  double *_statejab_32 = data->statejab_32;
  __attribute__ ((aligned (32))) double statejab_32[10][130];
  
  double *_stateh_12 = data->stateh_12;
  __attribute__ ((aligned (32))) double stateh_12[10][130];
  
  double *_stateh_22 = data->stateh_22;
  __attribute__ ((aligned (32))) double stateh_22[10][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[10][130];
  
  double *_stateuc_22 = data->stateuc_22;
  __attribute__ ((aligned (32))) double stateuc_22[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_32[swindex[j]][(i - 1)]) + sqrt((((statejab5_32[swindex[j]][(i - 1)] * 9.80616) * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)])));
            dql = ((stateh_12[swindex[j]][(i - 1)] - (4 * stateh_22[swindex[j]][(i - 1)])) + (3 * stateh_32[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateh_12[swindex[j]][i]) - (4 * stateh_22[swindex[j]][i])) + stateh_32[swindex[j]][i]);
            dfl = (((stateuc_12[swindex[j]][(i - 1)] * stateh_12[swindex[j]][(i - 1)]) - (4 * (stateuc_22[swindex[j]][(i - 1)] * stateh_22[swindex[j]][(i - 1)]))) + (3 * (stateuc_32[swindex[j]][(i - 1)] * stateh_32[swindex[j]][(i - 1)])));
            dfr = -(((3 * (stateuc_12[swindex[j]][i] * stateh_12[swindex[j]][i])) - (4 * (stateuc_22[swindex[j]][i] * stateh_22[swindex[j]][i]))) + (stateuc_32[swindex[j]][i] * stateh_32[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_37(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_37_info *data = (updateX_0_37_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[6][130];
  
  double *_statejab5_32 = data->statejab5_32;
  __attribute__ ((aligned (32))) double statejab5_32[6][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[6][130];
  
  double *_statejab_32 = data->statejab_32;
  __attribute__ ((aligned (32))) double statejab_32[6][130];
  
  double *_stateu_12 = data->stateu_12;
  __attribute__ ((aligned (32))) double stateu_12[6][130];
  
  double *_stateu_22 = data->stateu_22;
  __attribute__ ((aligned (32))) double stateu_22[6][130];
  
  double *_stateu_32 = data->stateu_32;
  __attribute__ ((aligned (32))) double stateu_32[6][130];
  
  double *_stateh_12 = data->stateh_12;
  __attribute__ ((aligned (32))) double stateh_12[6][130];
  
  double *_statejab_12 = data->statejab_12;
  __attribute__ ((aligned (32))) double statejab_12[6][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[6][130];
  
  double *_statev_12 = data->statev_12;
  __attribute__ ((aligned (32))) double statev_12[6][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[6][130];
  
  double *_stateh_22 = data->stateh_22;
  __attribute__ ((aligned (32))) double stateh_22[6][130];
  
  double *_statejab_22 = data->statejab_22;
  __attribute__ ((aligned (32))) double statejab_22[6][130];
  
  double *_stateuc_22 = data->stateuc_22;
  __attribute__ ((aligned (32))) double stateuc_22[6][130];
  
  double *_statev_22 = data->statev_22;
  __attribute__ ((aligned (32))) double statev_22[6][130];
  
  double *_statevc_22 = data->statevc_22;
  __attribute__ ((aligned (32))) double statevc_22[6][130];
  
  double *_statev_32 = data->statev_32;
  __attribute__ ((aligned (32))) double statev_32[6][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_22+(j*fnumx+ib), &statejab_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_22+(j*fnumx+ib), &statejab_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_32[swindex[j]][(i - 1)]) + sqrt((((statejab5_32[swindex[j]][(i - 1)] * 9.80616) * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)])));
            dql = ((stateu_12[swindex[j]][(i - 1)] - (4 * stateu_22[swindex[j]][(i - 1)])) + (3 * stateu_32[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateu_12[swindex[j]][i]) - (4 * stateu_22[swindex[j]][i])) + stateu_32[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_12[swindex[j]][(i - 1)]) / statejab_12[swindex[j]][(i - 1)]) + (0.5 * ((stateu_12[swindex[j]][(i - 1)] * stateuc_12[swindex[j]][(i - 1)]) + (statev_12[swindex[j]][(i - 1)] * statevc_12[swindex[j]][(i - 1)])))) - (4 * (((9.80616 * stateh_22[swindex[j]][(i - 1)]) / statejab_22[swindex[j]][(i - 1)]) + (0.5 * ((stateu_22[swindex[j]][(i - 1)] * stateuc_22[swindex[j]][(i - 1)]) + (statev_22[swindex[j]][(i - 1)] * statevc_22[swindex[j]][(i - 1)])))))) + (3 * (((9.80616 * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)]) + (0.5 * ((stateu_32[swindex[j]][(i - 1)] * stateuc_32[swindex[j]][(i - 1)]) + (statev_32[swindex[j]][(i - 1)] * statevc_32[swindex[j]][(i - 1)]))))));
            dfr = -(((3 * (((9.80616 * stateh_12[swindex[j]][i]) / statejab_12[swindex[j]][i]) + (0.5 * ((stateu_12[swindex[j]][i] * stateuc_12[swindex[j]][i]) + (statev_12[swindex[j]][i] * statevc_12[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_22[swindex[j]][i]) / statejab_22[swindex[j]][i]) + (0.5 * ((stateu_22[swindex[j]][i] * stateuc_22[swindex[j]][i]) + (statev_22[swindex[j]][i] * statevc_22[swindex[j]][i])))))) + (((9.80616 * stateh_32[swindex[j]][i]) / statejab_32[swindex[j]][i]) + (0.5 * ((stateu_32[swindex[j]][i] * stateuc_32[swindex[j]][i]) + (statev_32[swindex[j]][i] * statevc_32[swindex[j]][i])))));
            tendfu_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_38(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_38_info *data = (updateX_0_38_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[10][130];
  
  double *_statejab5_32 = data->statejab5_32;
  __attribute__ ((aligned (32))) double statejab5_32[10][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[10][130];
  
  double *_statejab_32 = data->statejab_32;
  __attribute__ ((aligned (32))) double statejab_32[10][130];
  
  double *_statev_12 = data->statev_12;
  __attribute__ ((aligned (32))) double statev_12[10][130];
  
  double *_statev_22 = data->statev_22;
  __attribute__ ((aligned (32))) double statev_22[10][130];
  
  double *_statev_32 = data->statev_32;
  __attribute__ ((aligned (32))) double statev_32[10][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[10][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_32+(j*fnumx+ib), &statejab5_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_32[swindex[j]][(i - 1)]) + sqrt((((statejab5_32[swindex[j]][(i - 1)] * 9.80616) * stateh_32[swindex[j]][(i - 1)]) / statejab_32[swindex[j]][(i - 1)])));
            dql = ((statev_12[swindex[j]][(i - 1)] - (4 * statev_22[swindex[j]][(i - 1)])) + (3 * statev_32[swindex[j]][(i - 1)]));
            dqr = -(((3 * statev_12[swindex[j]][i]) - (4 * statev_22[swindex[j]][i])) + statev_32[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfv_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_39(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_39_info *data = (updateX_0_39_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_22 = data->statevc_22;
  __attribute__ ((aligned (32))) double statevc_22[6][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_22 = data->stateuc_22;
  __attribute__ ((aligned (32))) double stateuc_22[6][130];
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[6][130];
  
  double *_tenddxh_12 = data->tenddxh_12;
  __attribute__ ((aligned (32))) double tenddxh_12[6][130];
  
  double *_tenddxh_22 = data->tenddxh_22;
  __attribute__ ((aligned (32))) double tenddxh_22[6][130];
  
  double *_tenddxh_32 = data->tenddxh_32;
  __attribute__ ((aligned (32))) double tenddxh_32[6][130];
  
  double *_tenddxu_12 = data->tenddxu_12;
  __attribute__ ((aligned (32))) double tenddxu_12[6][130];
  
  double *_tenddxu_22 = data->tenddxu_22;
  __attribute__ ((aligned (32))) double tenddxu_22[6][130];
  
  double *_tenddxu_32 = data->tenddxu_32;
  __attribute__ ((aligned (32))) double tenddxu_32[6][130];
  
  double *_tenddxv_12 = data->tenddxv_12;
  __attribute__ ((aligned (32))) double tenddxv_12[6][130];
  
  double *_tenddxv_22 = data->tenddxv_22;
  __attribute__ ((aligned (32))) double tenddxv_22[6][130];
  
  double *_tenddxv_32 = data->tenddxv_32;
  __attribute__ ((aligned (32))) double tenddxv_32[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqv_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966);
            sc = (((((3.0 * (tendqv_1[swindex[j]][(i + 1)] - tendqv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_12[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddxh_22[j][i] = (((((-3.0 * (tendfh_1[swindex[j]][(i + 1)] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_32[j][i] = (-tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966);
            tenddxu_12[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) + (statevc_12[swindex[j]][i] * sl));
            tenddxu_22[j][i] = ((((((3.0 * (tendfu_1[swindex[j]][(i + 1)] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) + (statevc_22[swindex[j]][i] * sc));
            tenddxu_32[j][i] = ((-tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) + (statevc_32[swindex[j]][i] * sr));
            tenddxv_12[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) - (stateuc_12[swindex[j]][i] * sl));
            tenddxv_22[j][i] = ((((((-3.0 * (tendfv_1[swindex[j]][(i + 1)] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) - (stateuc_22[swindex[j]][i] * sc));
            tenddxv_32[j][i] = ((-tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) - (stateuc_32[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddxh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_40(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_40_info *data = (updateX_0_40_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_statejab5_33 = data->statejab5_33;
  __attribute__ ((aligned (32))) double statejab5_33[6][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[6][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[6][130];
  
  double *_stateu_33 = data->stateu_33;
  __attribute__ ((aligned (32))) double stateu_33[6][130];
  
  double *_statev_33 = data->statev_33;
  __attribute__ ((aligned (32))) double statev_33[6][130];
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[6][130];
  
  double *_stateu_13 = data->stateu_13;
  __attribute__ ((aligned (32))) double stateu_13[6][130];
  
  double *_statev_13 = data->statev_13;
  __attribute__ ((aligned (32))) double statev_13[6][130];
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double ulocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_33[swindex[j]][(i - 1)]) + sqrt((((statejab5_33[swindex[j]][(i - 1)] * 9.80616) * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)])));
            tendfh_1[j][i] = ((((stateuc_33[swindex[j]][(i - 1)] * stateh_33[swindex[j]][(i - 1)]) + (stateuc_13[swindex[j]][i] * stateh_13[swindex[j]][i])) + (ulocal * (stateh_33[swindex[j]][(i - 1)] - stateh_13[swindex[j]][i]))) / 2);
            tendfu_1[j][i] = ((((((9.80616 * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)]) + (0.5 * ((stateu_33[swindex[j]][(i - 1)] * stateuc_33[swindex[j]][(i - 1)]) + (statev_33[swindex[j]][(i - 1)] * statevc_33[swindex[j]][(i - 1)])))) + (((9.80616 * stateh_13[swindex[j]][i]) / statejab_13[swindex[j]][i]) + (0.5 * ((stateu_13[swindex[j]][i] * stateuc_13[swindex[j]][i]) + (statev_13[swindex[j]][i] * statevc_13[swindex[j]][i]))))) + (ulocal * (stateu_33[swindex[j]][(i - 1)] - stateu_13[swindex[j]][i]))) / 2);
            tendfv_1[j][i] = (((0.0 + 0.0) + (ulocal * (statev_33[swindex[j]][(i - 1)] - statev_13[swindex[j]][i]))) / 2);
            tendqh_1[j][i] = ((stateh_33[swindex[j]][(i - 1)] + stateh_13[swindex[j]][i]) / 2);
            tendqu_1[j][i] = ((stateu_33[swindex[j]][(i - 1)] + stateu_13[swindex[j]][i]) / 2);
            tendqv_1[j][i] = ((statev_33[swindex[j]][(i - 1)] + statev_13[swindex[j]][i]) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_41(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_41_info *data = (updateX_0_41_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[10][130];
  
  double *_statejab5_33 = data->statejab5_33;
  __attribute__ ((aligned (32))) double statejab5_33[10][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[10][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[10][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[10][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[10][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[10][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_33[swindex[j]][(i - 1)]) + sqrt((((statejab5_33[swindex[j]][(i - 1)] * 9.80616) * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)])));
            dql = ((stateh_13[swindex[j]][(i - 1)] - (4 * stateh_23[swindex[j]][(i - 1)])) + (3 * stateh_33[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateh_13[swindex[j]][i]) - (4 * stateh_23[swindex[j]][i])) + stateh_33[swindex[j]][i]);
            dfl = (((stateuc_13[swindex[j]][(i - 1)] * stateh_13[swindex[j]][(i - 1)]) - (4 * (stateuc_23[swindex[j]][(i - 1)] * stateh_23[swindex[j]][(i - 1)]))) + (3 * (stateuc_33[swindex[j]][(i - 1)] * stateh_33[swindex[j]][(i - 1)])));
            dfr = -(((3 * (stateuc_13[swindex[j]][i] * stateh_13[swindex[j]][i])) - (4 * (stateuc_23[swindex[j]][i] * stateh_23[swindex[j]][i]))) + (stateuc_33[swindex[j]][i] * stateh_33[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_42(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_42_info *data = (updateX_0_42_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_statejab5_33 = data->statejab5_33;
  __attribute__ ((aligned (32))) double statejab5_33[6][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[6][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[6][130];
  
  double *_stateu_13 = data->stateu_13;
  __attribute__ ((aligned (32))) double stateu_13[6][130];
  
  double *_stateu_23 = data->stateu_23;
  __attribute__ ((aligned (32))) double stateu_23[6][130];
  
  double *_stateu_33 = data->stateu_33;
  __attribute__ ((aligned (32))) double stateu_33[6][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[6][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_statev_13 = data->statev_13;
  __attribute__ ((aligned (32))) double statev_13[6][130];
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[6][130];
  
  double *_statejab_23 = data->statejab_23;
  __attribute__ ((aligned (32))) double statejab_23[6][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[6][130];
  
  double *_statev_23 = data->statev_23;
  __attribute__ ((aligned (32))) double statev_23[6][130];
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[6][130];
  
  double *_statev_33 = data->statev_33;
  __attribute__ ((aligned (32))) double statev_33[6][130];
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_33[swindex[j]][(i - 1)]) + sqrt((((statejab5_33[swindex[j]][(i - 1)] * 9.80616) * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)])));
            dql = ((stateu_13[swindex[j]][(i - 1)] - (4 * stateu_23[swindex[j]][(i - 1)])) + (3 * stateu_33[swindex[j]][(i - 1)]));
            dqr = -(((3 * stateu_13[swindex[j]][i]) - (4 * stateu_23[swindex[j]][i])) + stateu_33[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_13[swindex[j]][(i - 1)]) / statejab_13[swindex[j]][(i - 1)]) + (0.5 * ((stateu_13[swindex[j]][(i - 1)] * stateuc_13[swindex[j]][(i - 1)]) + (statev_13[swindex[j]][(i - 1)] * statevc_13[swindex[j]][(i - 1)])))) - (4 * (((9.80616 * stateh_23[swindex[j]][(i - 1)]) / statejab_23[swindex[j]][(i - 1)]) + (0.5 * ((stateu_23[swindex[j]][(i - 1)] * stateuc_23[swindex[j]][(i - 1)]) + (statev_23[swindex[j]][(i - 1)] * statevc_23[swindex[j]][(i - 1)])))))) + (3 * (((9.80616 * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)]) + (0.5 * ((stateu_33[swindex[j]][(i - 1)] * stateuc_33[swindex[j]][(i - 1)]) + (statev_33[swindex[j]][(i - 1)] * statevc_33[swindex[j]][(i - 1)]))))));
            dfr = -(((3 * (((9.80616 * stateh_13[swindex[j]][i]) / statejab_13[swindex[j]][i]) + (0.5 * ((stateu_13[swindex[j]][i] * stateuc_13[swindex[j]][i]) + (statev_13[swindex[j]][i] * statevc_13[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_23[swindex[j]][i]) / statejab_23[swindex[j]][i]) + (0.5 * ((stateu_23[swindex[j]][i] * stateuc_23[swindex[j]][i]) + (statev_23[swindex[j]][i] * statevc_23[swindex[j]][i])))))) + (((9.80616 * stateh_33[swindex[j]][i]) / statejab_33[swindex[j]][i]) + (0.5 * ((stateu_33[swindex[j]][i] * stateuc_33[swindex[j]][i]) + (statev_33[swindex[j]][i] * statevc_33[swindex[j]][i])))));
            tendfu_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_43(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_43_info *data = (updateX_0_43_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[10][130];
  
  double *_statejab5_33 = data->statejab5_33;
  __attribute__ ((aligned (32))) double statejab5_33[10][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[10][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[10][130];
  
  double *_statev_13 = data->statev_13;
  __attribute__ ((aligned (32))) double statev_13[10][130];
  
  double *_statev_23 = data->statev_23;
  __attribute__ ((aligned (32))) double statev_23[10][130];
  
  double *_statev_33 = data->statev_33;
  __attribute__ ((aligned (32))) double statev_33[10][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[10][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[10][130];
  
  double ulocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab5_33+(j*fnumx+ib), &statejab5_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            ulocal = (fabs(stateuc_33[swindex[j]][(i - 1)]) + sqrt((((statejab5_33[swindex[j]][(i - 1)] * 9.80616) * stateh_33[swindex[j]][(i - 1)]) / statejab_33[swindex[j]][(i - 1)])));
            dql = ((statev_13[swindex[j]][(i - 1)] - (4 * statev_23[swindex[j]][(i - 1)])) + (3 * statev_33[swindex[j]][(i - 1)]));
            dqr = -(((3 * statev_13[swindex[j]][i]) - (4 * statev_23[swindex[j]][i])) + statev_33[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfv_2[j][i] = ((dfl + dfr) + ((ulocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateX_0_44(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateX_0_44_info *data = (updateX_0_44_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[6][130];
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[6][130];
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_tenddxh_13 = data->tenddxh_13;
  __attribute__ ((aligned (32))) double tenddxh_13[6][130];
  
  double *_tenddxh_23 = data->tenddxh_23;
  __attribute__ ((aligned (32))) double tenddxh_23[6][130];
  
  double *_tenddxh_33 = data->tenddxh_33;
  __attribute__ ((aligned (32))) double tenddxh_33[6][130];
  
  double *_tenddxu_13 = data->tenddxu_13;
  __attribute__ ((aligned (32))) double tenddxu_13[6][130];
  
  double *_tenddxu_23 = data->tenddxu_23;
  __attribute__ ((aligned (32))) double tenddxu_23[6][130];
  
  double *_tenddxu_33 = data->tenddxu_33;
  __attribute__ ((aligned (32))) double tenddxu_33[6][130];
  
  double *_tenddxv_13 = data->tenddxv_13;
  __attribute__ ((aligned (32))) double tenddxv_13[6][130];
  
  double *_tenddxv_23 = data->tenddxv_23;
  __attribute__ ((aligned (32))) double tenddxv_23[6][130];
  
  double *_tenddxv_33 = data->tenddxv_33;
  __attribute__ ((aligned (32))) double tenddxv_33[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqv_2+(j*fnumx+ib), &tendqv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqv_1+(j*fnumx+ib), &tendqv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqv_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966);
            sc = (((((3.0 * (tendqv_1[swindex[j]][(i + 1)] - tendqv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_13[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddxh_23[j][i] = (((((-3.0 * (tendfh_1[swindex[j]][(i + 1)] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0));
            tenddxh_33[j][i] = (-tendfh_2[swindex[j]][(i + 1)] / 16679.814955336966);
            tenddxu_13[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) + (statevc_13[swindex[j]][i] * sl));
            tenddxu_23[j][i] = ((((((3.0 * (tendfu_1[swindex[j]][(i + 1)] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) + (statevc_23[swindex[j]][i] * sc));
            tenddxu_33[j][i] = ((-tendfu_2[swindex[j]][(i + 1)] / 16679.814955336966) + (statevc_33[swindex[j]][i] * sr));
            tenddxv_13[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) - (stateuc_13[swindex[j]][i] * sl));
            tenddxv_23[j][i] = ((((((-3.0 * (tendfv_1[swindex[j]][(i + 1)] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) / 4.0)) - (stateuc_23[swindex[j]][i] * sc));
            tenddxv_33[j][i] = ((-tendfv_2[swindex[j]][(i + 1)] / 16679.814955336966) - (stateuc_33[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddxh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxh_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxh_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxu_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxu_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddxv_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddxv_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_30(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_30_info *data = (updateY_0_30_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_statejab7_13 = data->statejab7_13;
  __attribute__ ((aligned (32))) double statejab7_13[6][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[6][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[6][130];
  
  double *_stateu_13 = data->stateu_13;
  __attribute__ ((aligned (32))) double stateu_13[6][130];
  
  double *_stateu_11 = data->stateu_11;
  __attribute__ ((aligned (32))) double stateu_11[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_statev_13 = data->statev_13;
  __attribute__ ((aligned (32))) double statev_13[6][130];
  
  double *_statejab_11 = data->statejab_11;
  __attribute__ ((aligned (32))) double statejab_11[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_statev_11 = data->statev_11;
  __attribute__ ((aligned (32))) double statev_11[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double vlocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_13[swindex[(j - 1)]][i]) + sqrt((((statejab7_13[swindex[(j - 1)]][i] * 9.80616) * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i])));
            tendfh_1[j][i] = (((statevc_13[swindex[(j - 1)]][i] * stateh_13[swindex[(j - 1)]][i]) + (statevc_11[swindex[j]][i] * stateh_11[swindex[j]][i])) + ((vlocal * (stateh_13[swindex[(j - 1)]][i] - stateh_11[swindex[j]][i])) / 2.0));
            tendfu_1[j][i] = ((0.0 + 0.0) + ((vlocal * (stateu_13[swindex[(j - 1)]][i] - stateu_11[swindex[j]][i])) / 2.0));
            tendfv_1[j][i] = (((((9.80616 * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i]) + (0.5 * ((stateu_13[swindex[(j - 1)]][i] * stateuc_13[swindex[(j - 1)]][i]) + (statev_13[swindex[(j - 1)]][i] * statevc_13[swindex[(j - 1)]][i])))) + (((9.80616 * stateh_11[swindex[j]][i]) / statejab_11[swindex[j]][i]) + (0.5 * ((stateu_11[swindex[j]][i] * stateuc_11[swindex[j]][i]) + (statev_11[swindex[j]][i] * statevc_11[swindex[j]][i]))))) + ((vlocal * (statev_13[swindex[(j - 1)]][i] - statev_11[swindex[j]][i])) / 2.0));
            tendqh_1[j][i] = ((stateh_13[swindex[(j - 1)]][i] + stateh_11[swindex[j]][i]) / 2.0);
            tendqu_1[j][i] = ((stateu_13[swindex[(j - 1)]][i] + stateu_11[swindex[j]][i]) / 2.0);
            tendqv_1[j][i] = ((statev_13[swindex[(j - 1)]][i] + statev_11[swindex[j]][i]) / 2.0);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_31(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_31_info *data = (updateY_0_31_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[10][130];
  
  double *_statejab7_13 = data->statejab7_13;
  __attribute__ ((aligned (32))) double statejab7_13[10][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[10][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[10][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[10][130];
  
  double *_stateh_12 = data->stateh_12;
  __attribute__ ((aligned (32))) double stateh_12[10][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[10][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_13[swindex[(j - 1)]][i]) + sqrt((((statejab7_13[swindex[(j - 1)]][i] * 9.80616) * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i])));
            dql = ((stateh_11[swindex[(j - 1)]][i] - (4 * stateh_12[swindex[(j - 1)]][i])) + (3 * stateh_13[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateh_11[swindex[j]][i]) - (4 * stateh_12[swindex[j]][i])) + stateh_13[swindex[j]][i]);
            dfl = (((statevc_11[swindex[(j - 1)]][i] * stateh_11[swindex[(j - 1)]][i]) - (4 * (statevc_12[swindex[(j - 1)]][i] * stateh_12[swindex[(j - 1)]][i]))) + (3 * (statevc_13[swindex[(j - 1)]][i] * stateh_13[swindex[(j - 1)]][i])));
            dfr = -(((3 * (statevc_11[swindex[j]][i] * stateh_11[swindex[j]][i])) - (4 * (statevc_12[swindex[j]][i] * stateh_12[swindex[j]][i]))) + (statevc_13[swindex[j]][i] * stateh_13[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_32(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_32_info *data = (updateY_0_32_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[10][130];
  
  double *_statejab7_13 = data->statejab7_13;
  __attribute__ ((aligned (32))) double statejab7_13[10][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[10][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[10][130];
  
  double *_stateu_11 = data->stateu_11;
  __attribute__ ((aligned (32))) double stateu_11[10][130];
  
  double *_stateu_12 = data->stateu_12;
  __attribute__ ((aligned (32))) double stateu_12[10][130];
  
  double *_stateu_13 = data->stateu_13;
  __attribute__ ((aligned (32))) double stateu_13[10][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[10][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_13[swindex[(j - 1)]][i]) + sqrt((((statejab7_13[swindex[(j - 1)]][i] * 9.80616) * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i])));
            dql = ((stateu_11[swindex[(j - 1)]][i] - (4 * stateu_12[swindex[(j - 1)]][i])) + (3 * stateu_13[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateu_11[swindex[j]][i]) - (4 * stateu_12[swindex[j]][i])) + stateu_13[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfu_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_33(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_33_info *data = (updateY_0_33_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_statejab7_13 = data->statejab7_13;
  __attribute__ ((aligned (32))) double statejab7_13[6][130];
  
  double *_stateh_13 = data->stateh_13;
  __attribute__ ((aligned (32))) double stateh_13[6][130];
  
  double *_statejab_13 = data->statejab_13;
  __attribute__ ((aligned (32))) double statejab_13[6][130];
  
  double *_statev_11 = data->statev_11;
  __attribute__ ((aligned (32))) double statev_11[6][130];
  
  double *_statev_12 = data->statev_12;
  __attribute__ ((aligned (32))) double statev_12[6][130];
  
  double *_statev_13 = data->statev_13;
  __attribute__ ((aligned (32))) double statev_13[6][130];
  
  double *_stateh_11 = data->stateh_11;
  __attribute__ ((aligned (32))) double stateh_11[6][130];
  
  double *_statejab_11 = data->statejab_11;
  __attribute__ ((aligned (32))) double statejab_11[6][130];
  
  double *_stateu_11 = data->stateu_11;
  __attribute__ ((aligned (32))) double stateu_11[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_stateh_12 = data->stateh_12;
  __attribute__ ((aligned (32))) double stateh_12[6][130];
  
  double *_statejab_12 = data->statejab_12;
  __attribute__ ((aligned (32))) double statejab_12[6][130];
  
  double *_stateu_12 = data->stateu_12;
  __attribute__ ((aligned (32))) double stateu_12[6][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[6][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[6][130];
  
  double *_stateu_13 = data->stateu_13;
  __attribute__ ((aligned (32))) double stateu_13[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_13+(j*fnumx+ib), &statejab7_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_13+(j*fnumx+ib), &stateh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_13+(j*fnumx+ib), &statejab_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_11+(j*fnumx+ib), &statev_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_12+(j*fnumx+ib), &statev_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_13+(j*fnumx+ib), &statev_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_11+(j*fnumx+ib), &stateh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_11+(j*fnumx+ib), &statejab_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_11+(j*fnumx+ib), &stateu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_12+(j*fnumx+ib), &stateh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_12+(j*fnumx+ib), &statejab_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_12+(j*fnumx+ib), &stateu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_13+(j*fnumx+ib), &stateu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_13[swindex[(j - 1)]][i]) + sqrt((((statejab7_13[swindex[(j - 1)]][i] * 9.80616) * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i])));
            dql = ((statev_11[swindex[(j - 1)]][i] - (4 * statev_12[swindex[(j - 1)]][i])) + (3 * statev_13[swindex[(j - 1)]][i]));
            dqr = -(((3 * statev_11[swindex[j]][i]) - (4 * statev_12[swindex[j]][i])) + statev_13[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_11[swindex[(j - 1)]][i]) / statejab_11[swindex[(j - 1)]][i]) + (0.5 * ((stateu_11[swindex[(j - 1)]][i] * stateuc_11[swindex[(j - 1)]][i]) + (statev_11[swindex[(j - 1)]][i] * statevc_11[swindex[(j - 1)]][i])))) - (4 * (((9.80616 * stateh_12[swindex[(j - 1)]][i]) / statejab_12[swindex[(j - 1)]][i]) + (0.5 * ((stateu_12[swindex[(j - 1)]][i] * stateuc_12[swindex[(j - 1)]][i]) + (statev_12[swindex[(j - 1)]][i] * statevc_12[swindex[(j - 1)]][i])))))) + (3 * (((9.80616 * stateh_13[swindex[(j - 1)]][i]) / statejab_13[swindex[(j - 1)]][i]) + (0.5 * ((stateu_13[swindex[(j - 1)]][i] * stateuc_13[swindex[(j - 1)]][i]) + (statev_13[swindex[(j - 1)]][i] * statevc_13[swindex[(j - 1)]][i]))))));
            dfr = -(((3 * (((9.80616 * stateh_11[swindex[j]][i]) / statejab_11[swindex[j]][i]) + (0.5 * ((stateu_11[swindex[j]][i] * stateuc_11[swindex[j]][i]) + (statev_11[swindex[j]][i] * statevc_11[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_12[swindex[j]][i]) / statejab_12[swindex[j]][i]) + (0.5 * ((stateu_12[swindex[j]][i] * stateuc_12[swindex[j]][i]) + (statev_12[swindex[j]][i] * statevc_12[swindex[j]][i])))))) + (((9.80616 * stateh_13[swindex[j]][i]) / statejab_13[swindex[j]][i]) + (0.5 * ((stateu_13[swindex[j]][i] * stateuc_13[swindex[j]][i]) + (statev_13[swindex[j]][i] * statevc_13[swindex[j]][i])))));
            tendfv_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_34(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_34_info *data = (updateY_0_34_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_11 = data->statevc_11;
  __attribute__ ((aligned (32))) double statevc_11[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_12 = data->statevc_12;
  __attribute__ ((aligned (32))) double statevc_12[6][130];
  
  double *_statevc_13 = data->statevc_13;
  __attribute__ ((aligned (32))) double statevc_13[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_11 = data->stateuc_11;
  __attribute__ ((aligned (32))) double stateuc_11[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_12 = data->stateuc_12;
  __attribute__ ((aligned (32))) double stateuc_12[6][130];
  
  double *_stateuc_13 = data->stateuc_13;
  __attribute__ ((aligned (32))) double stateuc_13[6][130];
  
  double *_tenddyh_11 = data->tenddyh_11;
  __attribute__ ((aligned (32))) double tenddyh_11[6][130];
  
  double *_tenddyh_12 = data->tenddyh_12;
  __attribute__ ((aligned (32))) double tenddyh_12[6][130];
  
  double *_tenddyh_13 = data->tenddyh_13;
  __attribute__ ((aligned (32))) double tenddyh_13[6][130];
  
  double *_tenddyu_11 = data->tenddyu_11;
  __attribute__ ((aligned (32))) double tenddyu_11[6][130];
  
  double *_tenddyu_12 = data->tenddyu_12;
  __attribute__ ((aligned (32))) double tenddyu_12[6][130];
  
  double *_tenddyu_13 = data->tenddyu_13;
  __attribute__ ((aligned (32))) double tenddyu_13[6][130];
  
  double *_tenddyv_11 = data->tenddyv_11;
  __attribute__ ((aligned (32))) double tenddyv_11[6][130];
  
  double *_tenddyv_12 = data->tenddyv_12;
  __attribute__ ((aligned (32))) double tenddyv_12[6][130];
  
  double *_tenddyv_13 = data->tenddyv_13;
  __attribute__ ((aligned (32))) double tenddyv_13[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_11+(j*fnumx+ib), &statevc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_12+(j*fnumx+ib), &statevc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_13+(j*fnumx+ib), &statevc_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_11+(j*fnumx+ib), &stateuc_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_12+(j*fnumx+ib), &stateuc_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_13+(j*fnumx+ib), &stateuc_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqu_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966);
            sc = (((((3.0 * (tendqu_1[swindex[(j + 1)]][i] - tendqu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_11[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddyh_12[j][i] = (((((-3 * (tendfh_1[swindex[(j + 1)]][i] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_13[j][i] = (-tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966);
            tenddyu_11[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) - (statevc_11[swindex[j]][i] * sl));
            tenddyu_12[j][i] = ((((((-3.0 * (tendfu_1[swindex[(j + 1)]][i] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (statevc_12[swindex[j]][i] * sc));
            tenddyu_13[j][i] = ((-tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) - (statevc_13[swindex[j]][i] * sr));
            tenddyv_11[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) + (stateuc_11[swindex[j]][i] * sl));
            tenddyv_12[j][i] = ((((((-3.0 * (tendfv_1[swindex[(j + 1)]][i] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (stateuc_12[swindex[j]][i] * sc));
            tenddyv_13[j][i] = ((-tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) - (stateuc_13[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddyh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_35(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_35_info *data = (updateY_0_35_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[6][130];
  
  double *_statejab7_23 = data->statejab7_23;
  __attribute__ ((aligned (32))) double statejab7_23[6][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[6][130];
  
  double *_statejab_23 = data->statejab_23;
  __attribute__ ((aligned (32))) double statejab_23[6][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[6][130];
  
  double *_stateh_21 = data->stateh_21;
  __attribute__ ((aligned (32))) double stateh_21[6][130];
  
  double *_stateu_23 = data->stateu_23;
  __attribute__ ((aligned (32))) double stateu_23[6][130];
  
  double *_stateu_21 = data->stateu_21;
  __attribute__ ((aligned (32))) double stateu_21[6][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[6][130];
  
  double *_statev_23 = data->statev_23;
  __attribute__ ((aligned (32))) double statev_23[6][130];
  
  double *_statejab_21 = data->statejab_21;
  __attribute__ ((aligned (32))) double statejab_21[6][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[6][130];
  
  double *_statev_21 = data->statev_21;
  __attribute__ ((aligned (32))) double statev_21[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double vlocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_23[swindex[(j - 1)]][i]) + sqrt((((statejab7_23[swindex[(j - 1)]][i] * 9.80616) * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i])));
            tendfh_1[j][i] = (((statevc_23[swindex[(j - 1)]][i] * stateh_23[swindex[(j - 1)]][i]) + (statevc_21[swindex[j]][i] * stateh_21[swindex[j]][i])) + ((vlocal * (stateh_23[swindex[(j - 1)]][i] - stateh_21[swindex[j]][i])) / 2.0));
            tendfu_1[j][i] = ((0.0 + 0.0) + ((vlocal * (stateu_23[swindex[(j - 1)]][i] - stateu_21[swindex[j]][i])) / 2.0));
            tendfv_1[j][i] = (((((9.80616 * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i]) + (0.5 * ((stateu_23[swindex[(j - 1)]][i] * stateuc_23[swindex[(j - 1)]][i]) + (statev_23[swindex[(j - 1)]][i] * statevc_23[swindex[(j - 1)]][i])))) + (((9.80616 * stateh_21[swindex[j]][i]) / statejab_21[swindex[j]][i]) + (0.5 * ((stateu_21[swindex[j]][i] * stateuc_21[swindex[j]][i]) + (statev_21[swindex[j]][i] * statevc_21[swindex[j]][i]))))) + ((vlocal * (statev_23[swindex[(j - 1)]][i] - statev_21[swindex[j]][i])) / 2.0));
            tendqh_1[j][i] = ((stateh_23[swindex[(j - 1)]][i] + stateh_21[swindex[j]][i]) / 2.0);
            tendqu_1[j][i] = ((stateu_23[swindex[(j - 1)]][i] + stateu_21[swindex[j]][i]) / 2.0);
            tendqv_1[j][i] = ((statev_23[swindex[(j - 1)]][i] + statev_21[swindex[j]][i]) / 2.0);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_36(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_36_info *data = (updateY_0_36_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[10][130];
  
  double *_statejab7_23 = data->statejab7_23;
  __attribute__ ((aligned (32))) double statejab7_23[10][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[10][130];
  
  double *_statejab_23 = data->statejab_23;
  __attribute__ ((aligned (32))) double statejab_23[10][130];
  
  double *_stateh_21 = data->stateh_21;
  __attribute__ ((aligned (32))) double stateh_21[10][130];
  
  double *_stateh_22 = data->stateh_22;
  __attribute__ ((aligned (32))) double stateh_22[10][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[10][130];
  
  double *_statevc_22 = data->statevc_22;
  __attribute__ ((aligned (32))) double statevc_22[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_23[swindex[(j - 1)]][i]) + sqrt((((statejab7_23[swindex[(j - 1)]][i] * 9.80616) * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i])));
            dql = ((stateh_21[swindex[(j - 1)]][i] - (4 * stateh_22[swindex[(j - 1)]][i])) + (3 * stateh_23[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateh_21[swindex[j]][i]) - (4 * stateh_22[swindex[j]][i])) + stateh_23[swindex[j]][i]);
            dfl = (((statevc_21[swindex[(j - 1)]][i] * stateh_21[swindex[(j - 1)]][i]) - (4 * (statevc_22[swindex[(j - 1)]][i] * stateh_22[swindex[(j - 1)]][i]))) + (3 * (statevc_23[swindex[(j - 1)]][i] * stateh_23[swindex[(j - 1)]][i])));
            dfr = -(((3 * (statevc_21[swindex[j]][i] * stateh_21[swindex[j]][i])) - (4 * (statevc_22[swindex[j]][i] * stateh_22[swindex[j]][i]))) + (statevc_23[swindex[j]][i] * stateh_23[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_37(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_37_info *data = (updateY_0_37_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[10][130];
  
  double *_statejab7_23 = data->statejab7_23;
  __attribute__ ((aligned (32))) double statejab7_23[10][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[10][130];
  
  double *_statejab_23 = data->statejab_23;
  __attribute__ ((aligned (32))) double statejab_23[10][130];
  
  double *_stateu_21 = data->stateu_21;
  __attribute__ ((aligned (32))) double stateu_21[10][130];
  
  double *_stateu_22 = data->stateu_22;
  __attribute__ ((aligned (32))) double stateu_22[10][130];
  
  double *_stateu_23 = data->stateu_23;
  __attribute__ ((aligned (32))) double stateu_23[10][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[10][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_23[swindex[(j - 1)]][i]) + sqrt((((statejab7_23[swindex[(j - 1)]][i] * 9.80616) * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i])));
            dql = ((stateu_21[swindex[(j - 1)]][i] - (4 * stateu_22[swindex[(j - 1)]][i])) + (3 * stateu_23[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateu_21[swindex[j]][i]) - (4 * stateu_22[swindex[j]][i])) + stateu_23[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfu_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_38(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_38_info *data = (updateY_0_38_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[6][130];
  
  double *_statejab7_23 = data->statejab7_23;
  __attribute__ ((aligned (32))) double statejab7_23[6][130];
  
  double *_stateh_23 = data->stateh_23;
  __attribute__ ((aligned (32))) double stateh_23[6][130];
  
  double *_statejab_23 = data->statejab_23;
  __attribute__ ((aligned (32))) double statejab_23[6][130];
  
  double *_statev_21 = data->statev_21;
  __attribute__ ((aligned (32))) double statev_21[6][130];
  
  double *_statev_22 = data->statev_22;
  __attribute__ ((aligned (32))) double statev_22[6][130];
  
  double *_statev_23 = data->statev_23;
  __attribute__ ((aligned (32))) double statev_23[6][130];
  
  double *_stateh_21 = data->stateh_21;
  __attribute__ ((aligned (32))) double stateh_21[6][130];
  
  double *_statejab_21 = data->statejab_21;
  __attribute__ ((aligned (32))) double statejab_21[6][130];
  
  double *_stateu_21 = data->stateu_21;
  __attribute__ ((aligned (32))) double stateu_21[6][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[6][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[6][130];
  
  double *_stateh_22 = data->stateh_22;
  __attribute__ ((aligned (32))) double stateh_22[6][130];
  
  double *_statejab_22 = data->statejab_22;
  __attribute__ ((aligned (32))) double statejab_22[6][130];
  
  double *_stateu_22 = data->stateu_22;
  __attribute__ ((aligned (32))) double stateu_22[6][130];
  
  double *_stateuc_22 = data->stateuc_22;
  __attribute__ ((aligned (32))) double stateuc_22[6][130];
  
  double *_statevc_22 = data->statevc_22;
  __attribute__ ((aligned (32))) double statevc_22[6][130];
  
  double *_stateu_23 = data->stateu_23;
  __attribute__ ((aligned (32))) double stateu_23[6][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_22+(j*fnumx+ib), &statejab_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_23+(j*fnumx+ib), &statejab7_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_23+(j*fnumx+ib), &stateh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_23+(j*fnumx+ib), &statejab_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_21+(j*fnumx+ib), &statev_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_22+(j*fnumx+ib), &statev_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_23+(j*fnumx+ib), &statev_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_21+(j*fnumx+ib), &stateh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_21+(j*fnumx+ib), &statejab_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_21+(j*fnumx+ib), &stateu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_22+(j*fnumx+ib), &stateh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_22+(j*fnumx+ib), &statejab_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_22+(j*fnumx+ib), &stateu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_23+(j*fnumx+ib), &stateu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_23[swindex[(j - 1)]][i]) + sqrt((((statejab7_23[swindex[(j - 1)]][i] * 9.80616) * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i])));
            dql = ((statev_21[swindex[(j - 1)]][i] - (4 * statev_22[swindex[(j - 1)]][i])) + (3 * statev_23[swindex[(j - 1)]][i]));
            dqr = -(((3 * statev_21[swindex[j]][i]) - (4 * statev_22[swindex[j]][i])) + statev_23[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_21[swindex[(j - 1)]][i]) / statejab_21[swindex[(j - 1)]][i]) + (0.5 * ((stateu_21[swindex[(j - 1)]][i] * stateuc_21[swindex[(j - 1)]][i]) + (statev_21[swindex[(j - 1)]][i] * statevc_21[swindex[(j - 1)]][i])))) - (4 * (((9.80616 * stateh_22[swindex[(j - 1)]][i]) / statejab_22[swindex[(j - 1)]][i]) + (0.5 * ((stateu_22[swindex[(j - 1)]][i] * stateuc_22[swindex[(j - 1)]][i]) + (statev_22[swindex[(j - 1)]][i] * statevc_22[swindex[(j - 1)]][i])))))) + (3 * (((9.80616 * stateh_23[swindex[(j - 1)]][i]) / statejab_23[swindex[(j - 1)]][i]) + (0.5 * ((stateu_23[swindex[(j - 1)]][i] * stateuc_23[swindex[(j - 1)]][i]) + (statev_23[swindex[(j - 1)]][i] * statevc_23[swindex[(j - 1)]][i]))))));
            dfr = -(((3 * (((9.80616 * stateh_21[swindex[j]][i]) / statejab_21[swindex[j]][i]) + (0.5 * ((stateu_21[swindex[j]][i] * stateuc_21[swindex[j]][i]) + (statev_21[swindex[j]][i] * statevc_21[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_22[swindex[j]][i]) / statejab_22[swindex[j]][i]) + (0.5 * ((stateu_22[swindex[j]][i] * stateuc_22[swindex[j]][i]) + (statev_22[swindex[j]][i] * statevc_22[swindex[j]][i])))))) + (((9.80616 * stateh_23[swindex[j]][i]) / statejab_23[swindex[j]][i]) + (0.5 * ((stateu_23[swindex[j]][i] * stateuc_23[swindex[j]][i]) + (statev_23[swindex[j]][i] * statevc_23[swindex[j]][i])))));
            tendfv_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_39(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_39_info *data = (updateY_0_39_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_21 = data->statevc_21;
  __attribute__ ((aligned (32))) double statevc_21[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_22 = data->statevc_22;
  __attribute__ ((aligned (32))) double statevc_22[6][130];
  
  double *_statevc_23 = data->statevc_23;
  __attribute__ ((aligned (32))) double statevc_23[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_21 = data->stateuc_21;
  __attribute__ ((aligned (32))) double stateuc_21[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_22 = data->stateuc_22;
  __attribute__ ((aligned (32))) double stateuc_22[6][130];
  
  double *_stateuc_23 = data->stateuc_23;
  __attribute__ ((aligned (32))) double stateuc_23[6][130];
  
  double *_tenddyh_11 = data->tenddyh_11;
  __attribute__ ((aligned (32))) double tenddyh_11[6][130];
  
  double *_tenddyh_12 = data->tenddyh_12;
  __attribute__ ((aligned (32))) double tenddyh_12[6][130];
  
  double *_tenddyh_13 = data->tenddyh_13;
  __attribute__ ((aligned (32))) double tenddyh_13[6][130];
  
  double *_tenddyu_11 = data->tenddyu_11;
  __attribute__ ((aligned (32))) double tenddyu_11[6][130];
  
  double *_tenddyu_12 = data->tenddyu_12;
  __attribute__ ((aligned (32))) double tenddyu_12[6][130];
  
  double *_tenddyu_13 = data->tenddyu_13;
  __attribute__ ((aligned (32))) double tenddyu_13[6][130];
  
  double *_tenddyv_11 = data->tenddyv_11;
  __attribute__ ((aligned (32))) double tenddyv_11[6][130];
  
  double *_tenddyv_12 = data->tenddyv_12;
  __attribute__ ((aligned (32))) double tenddyv_12[6][130];
  
  double *_tenddyv_13 = data->tenddyv_13;
  __attribute__ ((aligned (32))) double tenddyv_13[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_21+(j*fnumx+ib), &statevc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_22+(j*fnumx+ib), &statevc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_23+(j*fnumx+ib), &statevc_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_21+(j*fnumx+ib), &stateuc_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_22+(j*fnumx+ib), &stateuc_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_23+(j*fnumx+ib), &stateuc_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqu_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966);
            sc = (((((3.0 * (tendqu_1[swindex[(j + 1)]][i] - tendqu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_11[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddyh_12[j][i] = (((((-3 * (tendfh_1[swindex[(j + 1)]][i] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_13[j][i] = (-tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966);
            tenddyu_11[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) - (statevc_21[swindex[j]][i] * sl));
            tenddyu_12[j][i] = ((((((-3.0 * (tendfu_1[swindex[(j + 1)]][i] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (statevc_22[swindex[j]][i] * sc));
            tenddyu_13[j][i] = ((-tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) - (statevc_23[swindex[j]][i] * sr));
            tenddyv_11[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) + (stateuc_21[swindex[j]][i] * sl));
            tenddyv_12[j][i] = ((((((-3.0 * (tendfv_1[swindex[(j + 1)]][i] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (stateuc_22[swindex[j]][i] * sc));
            tenddyv_13[j][i] = ((-tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) - (stateuc_23[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddyh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_40(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_40_info *data = (updateY_0_40_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_statejab7_33 = data->statejab7_33;
  __attribute__ ((aligned (32))) double statejab7_33[6][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[6][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[6][130];
  
  double *_stateu_33 = data->stateu_33;
  __attribute__ ((aligned (32))) double stateu_33[6][130];
  
  double *_stateu_31 = data->stateu_31;
  __attribute__ ((aligned (32))) double stateu_31[6][130];
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_statev_33 = data->statev_33;
  __attribute__ ((aligned (32))) double statev_33[6][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[6][130];
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_statev_31 = data->statev_31;
  __attribute__ ((aligned (32))) double statev_31[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_tendqh_1 = data->tendqh_1;
  __attribute__ ((aligned (32))) double tendqh_1[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendqv_1 = data->tendqv_1;
  __attribute__ ((aligned (32))) double tendqv_1[6][130];
  
  double vlocal;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_33[swindex[(j - 1)]][i]) + sqrt((((statejab7_33[swindex[(j - 1)]][i] * 9.80616) * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i])));
            tendfh_1[j][i] = (((statevc_33[swindex[(j - 1)]][i] * stateh_33[swindex[(j - 1)]][i]) + (statevc_31[swindex[j]][i] * stateh_31[swindex[j]][i])) + ((vlocal * (stateh_33[swindex[(j - 1)]][i] - stateh_31[swindex[j]][i])) / 2.0));
            tendfu_1[j][i] = ((0.0 + 0.0) + ((vlocal * (stateu_33[swindex[(j - 1)]][i] - stateu_31[swindex[j]][i])) / 2.0));
            tendfv_1[j][i] = (((((9.80616 * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i]) + (0.5 * ((stateu_33[swindex[(j - 1)]][i] * stateuc_33[swindex[(j - 1)]][i]) + (statev_33[swindex[(j - 1)]][i] * statevc_33[swindex[(j - 1)]][i])))) + (((9.80616 * stateh_31[swindex[j]][i]) / statejab_31[swindex[j]][i]) + (0.5 * ((stateu_31[swindex[j]][i] * stateuc_31[swindex[j]][i]) + (statev_31[swindex[j]][i] * statevc_31[swindex[j]][i]))))) + ((vlocal * (statev_33[swindex[(j - 1)]][i] - statev_31[swindex[j]][i])) / 2.0));
            tendqh_1[j][i] = ((stateh_33[swindex[(j - 1)]][i] + stateh_31[swindex[j]][i]) / 2.0);
            tendqu_1[j][i] = ((stateu_33[swindex[(j - 1)]][i] + stateu_31[swindex[j]][i]) / 2.0);
            tendqv_1[j][i] = ((statev_33[swindex[(j - 1)]][i] + statev_31[swindex[j]][i]) / 2.0);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendfv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_1+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_1[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_41(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_41_info *data = (updateY_0_41_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[10][130];
  
  double *_statejab7_33 = data->statejab7_33;
  __attribute__ ((aligned (32))) double statejab7_33[10][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[10][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[10][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[10][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[10][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[10][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[10][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[10][130];
  
  double *_tendqh_2 = data->tendqh_2;
  __attribute__ ((aligned (32))) double tendqh_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_33[swindex[(j - 1)]][i]) + sqrt((((statejab7_33[swindex[(j - 1)]][i] * 9.80616) * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i])));
            dql = ((stateh_31[swindex[(j - 1)]][i] - (4 * stateh_32[swindex[(j - 1)]][i])) + (3 * stateh_33[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateh_31[swindex[j]][i]) - (4 * stateh_32[swindex[j]][i])) + stateh_33[swindex[j]][i]);
            dfl = (((statevc_31[swindex[(j - 1)]][i] * stateh_31[swindex[(j - 1)]][i]) - (4 * (statevc_32[swindex[(j - 1)]][i] * stateh_32[swindex[(j - 1)]][i]))) + (3 * (statevc_33[swindex[(j - 1)]][i] * stateh_33[swindex[(j - 1)]][i])));
            dfr = -(((3 * (statevc_31[swindex[j]][i] * stateh_31[swindex[j]][i])) - (4 * (statevc_32[swindex[j]][i] * stateh_32[swindex[j]][i]))) + (statevc_33[swindex[j]][i] * stateh_33[swindex[j]][i]));
            tendfh_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqh_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqh_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqh_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_42(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_42_info *data = (updateY_0_42_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[10][130];
  
  double *_statejab7_33 = data->statejab7_33;
  __attribute__ ((aligned (32))) double statejab7_33[10][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[10][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[10][130];
  
  double *_stateu_31 = data->stateu_31;
  __attribute__ ((aligned (32))) double stateu_31[10][130];
  
  double *_stateu_32 = data->stateu_32;
  __attribute__ ((aligned (32))) double stateu_32[10][130];
  
  double *_stateu_33 = data->stateu_33;
  __attribute__ ((aligned (32))) double stateu_33[10][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[10][130];
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[10][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_33[swindex[(j - 1)]][i]) + sqrt((((statejab7_33[swindex[(j - 1)]][i] * 9.80616) * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i])));
            dql = ((stateu_31[swindex[(j - 1)]][i] - (4 * stateu_32[swindex[(j - 1)]][i])) + (3 * stateu_33[swindex[(j - 1)]][i]));
            dqr = -(((3 * stateu_31[swindex[j]][i]) - (4 * stateu_32[swindex[j]][i])) + stateu_33[swindex[j]][i]);
            dfl = 0.0;
            dfr = 0.0;
            tendfu_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqu_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqu_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqu_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_43(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_43_info *data = (updateY_0_43_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_statejab7_33 = data->statejab7_33;
  __attribute__ ((aligned (32))) double statejab7_33[6][130];
  
  double *_stateh_33 = data->stateh_33;
  __attribute__ ((aligned (32))) double stateh_33[6][130];
  
  double *_statejab_33 = data->statejab_33;
  __attribute__ ((aligned (32))) double statejab_33[6][130];
  
  double *_statev_31 = data->statev_31;
  __attribute__ ((aligned (32))) double statev_31[6][130];
  
  double *_statev_32 = data->statev_32;
  __attribute__ ((aligned (32))) double statev_32[6][130];
  
  double *_statev_33 = data->statev_33;
  __attribute__ ((aligned (32))) double statev_33[6][130];
  
  double *_stateh_31 = data->stateh_31;
  __attribute__ ((aligned (32))) double stateh_31[6][130];
  
  double *_statejab_31 = data->statejab_31;
  __attribute__ ((aligned (32))) double statejab_31[6][130];
  
  double *_stateu_31 = data->stateu_31;
  __attribute__ ((aligned (32))) double stateu_31[6][130];
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_stateh_32 = data->stateh_32;
  __attribute__ ((aligned (32))) double stateh_32[6][130];
  
  double *_statejab_32 = data->statejab_32;
  __attribute__ ((aligned (32))) double statejab_32[6][130];
  
  double *_stateu_32 = data->stateu_32;
  __attribute__ ((aligned (32))) double stateu_32[6][130];
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[6][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[6][130];
  
  double *_stateu_33 = data->stateu_33;
  __attribute__ ((aligned (32))) double stateu_33[6][130];
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_tendqv_2 = data->tendqv_2;
  __attribute__ ((aligned (32))) double tendqv_2[6][130];
  
  double vlocal;
  double dql;
  double dqr;
  double dfl;
  double dfr;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab7_33+(j*fnumx+ib), &statejab7_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_33+(j*fnumx+ib), &stateh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_33+(j*fnumx+ib), &statejab_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_31+(j*fnumx+ib), &statev_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_32+(j*fnumx+ib), &statev_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statev_33+(j*fnumx+ib), &statev_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_31+(j*fnumx+ib), &stateh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_31+(j*fnumx+ib), &statejab_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_31+(j*fnumx+ib), &stateu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateh_32+(j*fnumx+ib), &stateh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statejab_32+(j*fnumx+ib), &statejab_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_32+(j*fnumx+ib), &stateu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu_33+(j*fnumx+ib), &stateu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vlocal = (fabs(statevc_33[swindex[(j - 1)]][i]) + sqrt((((statejab7_33[swindex[(j - 1)]][i] * 9.80616) * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i])));
            dql = ((statev_31[swindex[(j - 1)]][i] - (4 * statev_32[swindex[(j - 1)]][i])) + (3 * statev_33[swindex[(j - 1)]][i]));
            dqr = -(((3 * statev_31[swindex[j]][i]) - (4 * statev_32[swindex[j]][i])) + statev_33[swindex[j]][i]);
            dfl = (((((9.80616 * stateh_31[swindex[(j - 1)]][i]) / statejab_31[swindex[(j - 1)]][i]) + (0.5 * ((stateu_31[swindex[(j - 1)]][i] * stateuc_31[swindex[(j - 1)]][i]) + (statev_31[swindex[(j - 1)]][i] * statevc_31[swindex[(j - 1)]][i])))) - (4 * (((9.80616 * stateh_32[swindex[(j - 1)]][i]) / statejab_32[swindex[(j - 1)]][i]) + (0.5 * ((stateu_32[swindex[(j - 1)]][i] * stateuc_32[swindex[(j - 1)]][i]) + (statev_32[swindex[(j - 1)]][i] * statevc_32[swindex[(j - 1)]][i])))))) + (3 * (((9.80616 * stateh_33[swindex[(j - 1)]][i]) / statejab_33[swindex[(j - 1)]][i]) + (0.5 * ((stateu_33[swindex[(j - 1)]][i] * stateuc_33[swindex[(j - 1)]][i]) + (statev_33[swindex[(j - 1)]][i] * statevc_33[swindex[(j - 1)]][i]))))));
            dfr = -(((3 * (((9.80616 * stateh_31[swindex[j]][i]) / statejab_31[swindex[j]][i]) + (0.5 * ((stateu_31[swindex[j]][i] * stateuc_31[swindex[j]][i]) + (statev_31[swindex[j]][i] * statevc_31[swindex[j]][i]))))) - (4 * (((9.80616 * stateh_32[swindex[j]][i]) / statejab_32[swindex[j]][i]) + (0.5 * ((stateu_32[swindex[j]][i] * stateuc_32[swindex[j]][i]) + (statev_32[swindex[j]][i] * statevc_32[swindex[j]][i])))))) + (((9.80616 * stateh_33[swindex[j]][i]) / statejab_33[swindex[j]][i]) + (0.5 * ((stateu_33[swindex[j]][i] * stateuc_33[swindex[j]][i]) + (statev_33[swindex[j]][i] * statevc_33[swindex[j]][i])))));
            tendfv_2[j][i] = ((dfl + dfr) + ((vlocal * (dql - dqr)) / 2));
            tendqv_2[j][i] = ((dql + dqr) / 2);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tendfv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendfv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tendqv_2+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tendqv_2[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void updateY_0_44(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  updateY_0_44_info *data = (updateY_0_44_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_tendqu_2 = data->tendqu_2;
  __attribute__ ((aligned (32))) double tendqu_2[6][130];
  
  double *_tendqu_1 = data->tendqu_1;
  __attribute__ ((aligned (32))) double tendqu_1[6][130];
  
  double *_tendfh_2 = data->tendfh_2;
  __attribute__ ((aligned (32))) double tendfh_2[6][130];
  
  double *_tendfh_1 = data->tendfh_1;
  __attribute__ ((aligned (32))) double tendfh_1[6][130];
  
  double *_tendfu_2 = data->tendfu_2;
  __attribute__ ((aligned (32))) double tendfu_2[6][130];
  
  double *_statevc_31 = data->statevc_31;
  __attribute__ ((aligned (32))) double statevc_31[6][130];
  
  double *_tendfu_1 = data->tendfu_1;
  __attribute__ ((aligned (32))) double tendfu_1[6][130];
  
  double *_statevc_32 = data->statevc_32;
  __attribute__ ((aligned (32))) double statevc_32[6][130];
  
  double *_statevc_33 = data->statevc_33;
  __attribute__ ((aligned (32))) double statevc_33[6][130];
  
  double *_tendfv_2 = data->tendfv_2;
  __attribute__ ((aligned (32))) double tendfv_2[6][130];
  
  double *_stateuc_31 = data->stateuc_31;
  __attribute__ ((aligned (32))) double stateuc_31[6][130];
  
  double *_tendfv_1 = data->tendfv_1;
  __attribute__ ((aligned (32))) double tendfv_1[6][130];
  
  double *_stateuc_32 = data->stateuc_32;
  __attribute__ ((aligned (32))) double stateuc_32[6][130];
  
  double *_stateuc_33 = data->stateuc_33;
  __attribute__ ((aligned (32))) double stateuc_33[6][130];
  
  double *_tenddyh_11 = data->tenddyh_11;
  __attribute__ ((aligned (32))) double tenddyh_11[6][130];
  
  double *_tenddyh_12 = data->tenddyh_12;
  __attribute__ ((aligned (32))) double tenddyh_12[6][130];
  
  double *_tenddyh_13 = data->tenddyh_13;
  __attribute__ ((aligned (32))) double tenddyh_13[6][130];
  
  double *_tenddyu_11 = data->tenddyu_11;
  __attribute__ ((aligned (32))) double tenddyu_11[6][130];
  
  double *_tenddyu_12 = data->tenddyu_12;
  __attribute__ ((aligned (32))) double tenddyu_12[6][130];
  
  double *_tenddyu_13 = data->tenddyu_13;
  __attribute__ ((aligned (32))) double tenddyu_13[6][130];
  
  double *_tenddyv_11 = data->tenddyv_11;
  __attribute__ ((aligned (32))) double tenddyv_11[6][130];
  
  double *_tenddyv_12 = data->tenddyv_12;
  __attribute__ ((aligned (32))) double tenddyv_12[6][130];
  
  double *_tenddyv_13 = data->tenddyv_13;
  __attribute__ ((aligned (32))) double tenddyv_13[6][130];
  
  double sl;
  double sr;
  double sc;
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_tendqu_2+(j*fnumx+ib), &tendqu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendqu_1+(j*fnumx+ib), &tendqu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_2+(j*fnumx+ib), &tendfh_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfh_1+(j*fnumx+ib), &tendfh_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_2+(j*fnumx+ib), &tendfu_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_31+(j*fnumx+ib), &statevc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfu_1+(j*fnumx+ib), &tendfu_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_32+(j*fnumx+ib), &statevc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statevc_33+(j*fnumx+ib), &statevc_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_2+(j*fnumx+ib), &tendfv_2[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_31+(j*fnumx+ib), &stateuc_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tendfv_1+(j*fnumx+ib), &tendfv_1[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_32+(j*fnumx+ib), &stateuc_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateuc_33+(j*fnumx+ib), &stateuc_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            sl = (tendqu_2[swindex[j]][i] / 16679.814955336966);
            sr = (tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966);
            sc = (((((3.0 * (tendqu_1[swindex[(j + 1)]][i] - tendqu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) - ((tendqu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) - ((tendqu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_11[j][i] = (-tendfh_2[swindex[j]][i] / 16679.814955336966);
            tenddyh_12[j][i] = (((((-3 * (tendfh_1[swindex[(j + 1)]][i] - tendfh_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfh_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0));
            tenddyh_13[j][i] = (-tendfh_2[swindex[(j + 1)]][i] / 16679.814955336966);
            tenddyu_11[j][i] = ((-tendfu_2[swindex[j]][i] / 16679.814955336966) - (statevc_31[swindex[j]][i] * sl));
            tenddyu_12[j][i] = ((((((-3.0 * (tendfu_1[swindex[(j + 1)]][i] - tendfu_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfu_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (statevc_32[swindex[j]][i] * sc));
            tenddyu_13[j][i] = ((-tendfu_2[swindex[(j + 1)]][i] / 16679.814955336966) - (statevc_33[swindex[j]][i] * sr));
            tenddyv_11[j][i] = ((-tendfv_2[swindex[j]][i] / 16679.814955336966) + (stateuc_31[swindex[j]][i] * sl));
            tenddyv_12[j][i] = ((((((-3.0 * (tendfv_1[swindex[(j + 1)]][i] - tendfv_1[swindex[j]][i])) / 16679.814955336966) / 2.0) + ((tendfv_2[swindex[j]][i] / 16679.814955336966) / 4.0)) + ((tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) / 4.0)) - (stateuc_32[swindex[j]][i] * sc));
            tenddyv_13[j][i] = ((-tendfv_2[swindex[(j + 1)]][i] / 16679.814955336966) - (stateuc_33[swindex[j]][i] * sr));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddyh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddyv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddyv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

void update_rk1_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_0_info *data = (update_rk1_0_0_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_11 = data->state_oldh_11;
  __attribute__ ((aligned (32))) double state_oldh_11[10][130];
  
  double *_tend1dxh_11 = data->tend1dxh_11;
  __attribute__ ((aligned (32))) double tend1dxh_11[10][130];
  
  double *_tend1dyh_11 = data->tend1dyh_11;
  __attribute__ ((aligned (32))) double tend1dyh_11[10][130];
  
  double *_state_oldu_11 = data->state_oldu_11;
  __attribute__ ((aligned (32))) double state_oldu_11[10][130];
  
  double *_tend1dxu_11 = data->tend1dxu_11;
  __attribute__ ((aligned (32))) double tend1dxu_11[10][130];
  
  double *_tend1dyu_11 = data->tend1dyu_11;
  __attribute__ ((aligned (32))) double tend1dyu_11[10][130];
  
  double *_state_oldv_11 = data->state_oldv_11;
  __attribute__ ((aligned (32))) double state_oldv_11[10][130];
  
  double *_tend1dxv_11 = data->tend1dxv_11;
  __attribute__ ((aligned (32))) double tend1dxv_11[10][130];
  
  double *_tend1dyv_11 = data->tend1dyv_11;
  __attribute__ ((aligned (32))) double tend1dyv_11[10][130];
  
  double *_state_newh_11 = data->state_newh_11;
  __attribute__ ((aligned (32))) double state_newh_11[10][130];
  
  double *_state_newu_11 = data->state_newu_11;
  __attribute__ ((aligned (32))) double state_newu_11[10][130];
  
  double *_state_newv_11 = data->state_newv_11;
  __attribute__ ((aligned (32))) double state_newv_11[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_11[j][i] = (state_oldh_11[swindex[j]][i] + ((tend1dxh_11[swindex[j]][i] + tend1dyh_11[swindex[j]][i]) * dt));
            state_newu_11[j][i] = (state_oldu_11[swindex[j]][i] + ((tend1dxu_11[swindex[j]][i] + tend1dyu_11[swindex[j]][i]) * dt));
            state_newv_11[j][i] = (state_oldv_11[swindex[j]][i] + ((tend1dxv_11[swindex[j]][i] + tend1dyv_11[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_1(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_1_info *data = (update_rk1_0_1_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_12 = data->state_oldh_12;
  __attribute__ ((aligned (32))) double state_oldh_12[10][130];
  
  double *_tend1dxh_12 = data->tend1dxh_12;
  __attribute__ ((aligned (32))) double tend1dxh_12[10][130];
  
  double *_tend1dyh_12 = data->tend1dyh_12;
  __attribute__ ((aligned (32))) double tend1dyh_12[10][130];
  
  double *_state_oldu_12 = data->state_oldu_12;
  __attribute__ ((aligned (32))) double state_oldu_12[10][130];
  
  double *_tend1dxu_12 = data->tend1dxu_12;
  __attribute__ ((aligned (32))) double tend1dxu_12[10][130];
  
  double *_tend1dyu_12 = data->tend1dyu_12;
  __attribute__ ((aligned (32))) double tend1dyu_12[10][130];
  
  double *_state_oldv_12 = data->state_oldv_12;
  __attribute__ ((aligned (32))) double state_oldv_12[10][130];
  
  double *_tend1dxv_12 = data->tend1dxv_12;
  __attribute__ ((aligned (32))) double tend1dxv_12[10][130];
  
  double *_tend1dyv_12 = data->tend1dyv_12;
  __attribute__ ((aligned (32))) double tend1dyv_12[10][130];
  
  double *_state_newh_12 = data->state_newh_12;
  __attribute__ ((aligned (32))) double state_newh_12[10][130];
  
  double *_state_newu_12 = data->state_newu_12;
  __attribute__ ((aligned (32))) double state_newu_12[10][130];
  
  double *_state_newv_12 = data->state_newv_12;
  __attribute__ ((aligned (32))) double state_newv_12[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_12[j][i] = (state_oldh_12[swindex[j]][i] + ((tend1dxh_12[swindex[j]][i] + tend1dyh_12[swindex[j]][i]) * dt));
            state_newu_12[j][i] = (state_oldu_12[swindex[j]][i] + ((tend1dxu_12[swindex[j]][i] + tend1dyu_12[swindex[j]][i]) * dt));
            state_newv_12[j][i] = (state_oldv_12[swindex[j]][i] + ((tend1dxv_12[swindex[j]][i] + tend1dyv_12[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_2(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_2_info *data = (update_rk1_0_2_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_13 = data->state_oldh_13;
  __attribute__ ((aligned (32))) double state_oldh_13[10][130];
  
  double *_tend1dxh_13 = data->tend1dxh_13;
  __attribute__ ((aligned (32))) double tend1dxh_13[10][130];
  
  double *_tend1dyh_13 = data->tend1dyh_13;
  __attribute__ ((aligned (32))) double tend1dyh_13[10][130];
  
  double *_state_oldu_13 = data->state_oldu_13;
  __attribute__ ((aligned (32))) double state_oldu_13[10][130];
  
  double *_tend1dxu_13 = data->tend1dxu_13;
  __attribute__ ((aligned (32))) double tend1dxu_13[10][130];
  
  double *_tend1dyu_13 = data->tend1dyu_13;
  __attribute__ ((aligned (32))) double tend1dyu_13[10][130];
  
  double *_state_oldv_13 = data->state_oldv_13;
  __attribute__ ((aligned (32))) double state_oldv_13[10][130];
  
  double *_tend1dxv_13 = data->tend1dxv_13;
  __attribute__ ((aligned (32))) double tend1dxv_13[10][130];
  
  double *_tend1dyv_13 = data->tend1dyv_13;
  __attribute__ ((aligned (32))) double tend1dyv_13[10][130];
  
  double *_state_newh_13 = data->state_newh_13;
  __attribute__ ((aligned (32))) double state_newh_13[10][130];
  
  double *_state_newu_13 = data->state_newu_13;
  __attribute__ ((aligned (32))) double state_newu_13[10][130];
  
  double *_state_newv_13 = data->state_newv_13;
  __attribute__ ((aligned (32))) double state_newv_13[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_13[j][i] = (state_oldh_13[swindex[j]][i] + ((tend1dxh_13[swindex[j]][i] + tend1dyh_13[swindex[j]][i]) * dt));
            state_newu_13[j][i] = (state_oldu_13[swindex[j]][i] + ((tend1dxu_13[swindex[j]][i] + tend1dyu_13[swindex[j]][i]) * dt));
            state_newv_13[j][i] = (state_oldv_13[swindex[j]][i] + ((tend1dxv_13[swindex[j]][i] + tend1dyv_13[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_3(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_3_info *data = (update_rk1_0_3_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_21 = data->state_oldh_21;
  __attribute__ ((aligned (32))) double state_oldh_21[10][130];
  
  double *_tend1dxh_21 = data->tend1dxh_21;
  __attribute__ ((aligned (32))) double tend1dxh_21[10][130];
  
  double *_tend1dyh_21 = data->tend1dyh_21;
  __attribute__ ((aligned (32))) double tend1dyh_21[10][130];
  
  double *_state_oldu_21 = data->state_oldu_21;
  __attribute__ ((aligned (32))) double state_oldu_21[10][130];
  
  double *_tend1dxu_21 = data->tend1dxu_21;
  __attribute__ ((aligned (32))) double tend1dxu_21[10][130];
  
  double *_tend1dyu_21 = data->tend1dyu_21;
  __attribute__ ((aligned (32))) double tend1dyu_21[10][130];
  
  double *_state_oldv_21 = data->state_oldv_21;
  __attribute__ ((aligned (32))) double state_oldv_21[10][130];
  
  double *_tend1dxv_21 = data->tend1dxv_21;
  __attribute__ ((aligned (32))) double tend1dxv_21[10][130];
  
  double *_tend1dyv_21 = data->tend1dyv_21;
  __attribute__ ((aligned (32))) double tend1dyv_21[10][130];
  
  double *_state_newh_21 = data->state_newh_21;
  __attribute__ ((aligned (32))) double state_newh_21[10][130];
  
  double *_state_newu_21 = data->state_newu_21;
  __attribute__ ((aligned (32))) double state_newu_21[10][130];
  
  double *_state_newv_21 = data->state_newv_21;
  __attribute__ ((aligned (32))) double state_newv_21[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_21[j][i] = (state_oldh_21[swindex[j]][i] + ((tend1dxh_21[swindex[j]][i] + tend1dyh_21[swindex[j]][i]) * dt));
            state_newu_21[j][i] = (state_oldu_21[swindex[j]][i] + ((tend1dxu_21[swindex[j]][i] + tend1dyu_21[swindex[j]][i]) * dt));
            state_newv_21[j][i] = (state_oldv_21[swindex[j]][i] + ((tend1dxv_21[swindex[j]][i] + tend1dyv_21[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_4(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_4_info *data = (update_rk1_0_4_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_22 = data->state_oldh_22;
  __attribute__ ((aligned (32))) double state_oldh_22[10][130];
  
  double *_tend1dxh_22 = data->tend1dxh_22;
  __attribute__ ((aligned (32))) double tend1dxh_22[10][130];
  
  double *_tend1dyh_22 = data->tend1dyh_22;
  __attribute__ ((aligned (32))) double tend1dyh_22[10][130];
  
  double *_state_oldu_22 = data->state_oldu_22;
  __attribute__ ((aligned (32))) double state_oldu_22[10][130];
  
  double *_tend1dxu_22 = data->tend1dxu_22;
  __attribute__ ((aligned (32))) double tend1dxu_22[10][130];
  
  double *_tend1dyu_22 = data->tend1dyu_22;
  __attribute__ ((aligned (32))) double tend1dyu_22[10][130];
  
  double *_state_oldv_22 = data->state_oldv_22;
  __attribute__ ((aligned (32))) double state_oldv_22[10][130];
  
  double *_tend1dxv_22 = data->tend1dxv_22;
  __attribute__ ((aligned (32))) double tend1dxv_22[10][130];
  
  double *_tend1dyv_22 = data->tend1dyv_22;
  __attribute__ ((aligned (32))) double tend1dyv_22[10][130];
  
  double *_state_newh_22 = data->state_newh_22;
  __attribute__ ((aligned (32))) double state_newh_22[10][130];
  
  double *_state_newu_22 = data->state_newu_22;
  __attribute__ ((aligned (32))) double state_newu_22[10][130];
  
  double *_state_newv_22 = data->state_newv_22;
  __attribute__ ((aligned (32))) double state_newv_22[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_22[j][i] = (state_oldh_22[swindex[j]][i] + ((tend1dxh_22[swindex[j]][i] + tend1dyh_22[swindex[j]][i]) * dt));
            state_newu_22[j][i] = (state_oldu_22[swindex[j]][i] + ((tend1dxu_22[swindex[j]][i] + tend1dyu_22[swindex[j]][i]) * dt));
            state_newv_22[j][i] = (state_oldv_22[swindex[j]][i] + ((tend1dxv_22[swindex[j]][i] + tend1dyv_22[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_5(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_5_info *data = (update_rk1_0_5_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_23 = data->state_oldh_23;
  __attribute__ ((aligned (32))) double state_oldh_23[10][130];
  
  double *_tend1dxh_23 = data->tend1dxh_23;
  __attribute__ ((aligned (32))) double tend1dxh_23[10][130];
  
  double *_tend1dyh_23 = data->tend1dyh_23;
  __attribute__ ((aligned (32))) double tend1dyh_23[10][130];
  
  double *_state_oldu_23 = data->state_oldu_23;
  __attribute__ ((aligned (32))) double state_oldu_23[10][130];
  
  double *_tend1dxu_23 = data->tend1dxu_23;
  __attribute__ ((aligned (32))) double tend1dxu_23[10][130];
  
  double *_tend1dyu_23 = data->tend1dyu_23;
  __attribute__ ((aligned (32))) double tend1dyu_23[10][130];
  
  double *_state_oldv_23 = data->state_oldv_23;
  __attribute__ ((aligned (32))) double state_oldv_23[10][130];
  
  double *_tend1dxv_23 = data->tend1dxv_23;
  __attribute__ ((aligned (32))) double tend1dxv_23[10][130];
  
  double *_tend1dyv_23 = data->tend1dyv_23;
  __attribute__ ((aligned (32))) double tend1dyv_23[10][130];
  
  double *_state_newh_23 = data->state_newh_23;
  __attribute__ ((aligned (32))) double state_newh_23[10][130];
  
  double *_state_newu_23 = data->state_newu_23;
  __attribute__ ((aligned (32))) double state_newu_23[10][130];
  
  double *_state_newv_23 = data->state_newv_23;
  __attribute__ ((aligned (32))) double state_newv_23[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_23[j][i] = (state_oldh_23[swindex[j]][i] + ((tend1dxh_23[swindex[j]][i] + tend1dyh_23[swindex[j]][i]) * dt));
            state_newu_23[j][i] = (state_oldu_23[swindex[j]][i] + ((tend1dxu_23[swindex[j]][i] + tend1dyu_23[swindex[j]][i]) * dt));
            state_newv_23[j][i] = (state_oldv_23[swindex[j]][i] + ((tend1dxv_23[swindex[j]][i] + tend1dyv_23[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_6(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_6_info *data = (update_rk1_0_6_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_31 = data->state_oldh_31;
  __attribute__ ((aligned (32))) double state_oldh_31[10][130];
  
  double *_tend1dxh_31 = data->tend1dxh_31;
  __attribute__ ((aligned (32))) double tend1dxh_31[10][130];
  
  double *_tend1dyh_31 = data->tend1dyh_31;
  __attribute__ ((aligned (32))) double tend1dyh_31[10][130];
  
  double *_state_oldu_31 = data->state_oldu_31;
  __attribute__ ((aligned (32))) double state_oldu_31[10][130];
  
  double *_tend1dxu_31 = data->tend1dxu_31;
  __attribute__ ((aligned (32))) double tend1dxu_31[10][130];
  
  double *_tend1dyu_31 = data->tend1dyu_31;
  __attribute__ ((aligned (32))) double tend1dyu_31[10][130];
  
  double *_state_oldv_31 = data->state_oldv_31;
  __attribute__ ((aligned (32))) double state_oldv_31[10][130];
  
  double *_tend1dxv_31 = data->tend1dxv_31;
  __attribute__ ((aligned (32))) double tend1dxv_31[10][130];
  
  double *_tend1dyv_31 = data->tend1dyv_31;
  __attribute__ ((aligned (32))) double tend1dyv_31[10][130];
  
  double *_state_newh_31 = data->state_newh_31;
  __attribute__ ((aligned (32))) double state_newh_31[10][130];
  
  double *_state_newu_31 = data->state_newu_31;
  __attribute__ ((aligned (32))) double state_newu_31[10][130];
  
  double *_state_newv_31 = data->state_newv_31;
  __attribute__ ((aligned (32))) double state_newv_31[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_31[j][i] = (state_oldh_31[swindex[j]][i] + ((tend1dxh_31[swindex[j]][i] + tend1dyh_31[swindex[j]][i]) * dt));
            state_newu_31[j][i] = (state_oldu_31[swindex[j]][i] + ((tend1dxu_31[swindex[j]][i] + tend1dyu_31[swindex[j]][i]) * dt));
            state_newv_31[j][i] = (state_oldv_31[swindex[j]][i] + ((tend1dxv_31[swindex[j]][i] + tend1dyv_31[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_7(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_7_info *data = (update_rk1_0_7_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_32 = data->state_oldh_32;
  __attribute__ ((aligned (32))) double state_oldh_32[10][130];
  
  double *_tend1dxh_32 = data->tend1dxh_32;
  __attribute__ ((aligned (32))) double tend1dxh_32[10][130];
  
  double *_tend1dyh_32 = data->tend1dyh_32;
  __attribute__ ((aligned (32))) double tend1dyh_32[10][130];
  
  double *_state_oldu_32 = data->state_oldu_32;
  __attribute__ ((aligned (32))) double state_oldu_32[10][130];
  
  double *_tend1dxu_32 = data->tend1dxu_32;
  __attribute__ ((aligned (32))) double tend1dxu_32[10][130];
  
  double *_tend1dyu_32 = data->tend1dyu_32;
  __attribute__ ((aligned (32))) double tend1dyu_32[10][130];
  
  double *_state_oldv_32 = data->state_oldv_32;
  __attribute__ ((aligned (32))) double state_oldv_32[10][130];
  
  double *_tend1dxv_32 = data->tend1dxv_32;
  __attribute__ ((aligned (32))) double tend1dxv_32[10][130];
  
  double *_tend1dyv_32 = data->tend1dyv_32;
  __attribute__ ((aligned (32))) double tend1dyv_32[10][130];
  
  double *_state_newh_32 = data->state_newh_32;
  __attribute__ ((aligned (32))) double state_newh_32[10][130];
  
  double *_state_newu_32 = data->state_newu_32;
  __attribute__ ((aligned (32))) double state_newu_32[10][130];
  
  double *_state_newv_32 = data->state_newv_32;
  __attribute__ ((aligned (32))) double state_newv_32[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_32[j][i] = (state_oldh_32[swindex[j]][i] + ((tend1dxh_32[swindex[j]][i] + tend1dyh_32[swindex[j]][i]) * dt));
            state_newu_32[j][i] = (state_oldu_32[swindex[j]][i] + ((tend1dxu_32[swindex[j]][i] + tend1dyu_32[swindex[j]][i]) * dt));
            state_newv_32[j][i] = (state_oldv_32[swindex[j]][i] + ((tend1dxv_32[swindex[j]][i] + tend1dyv_32[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk1_0_8(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk1_0_8_info *data = (update_rk1_0_8_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_33 = data->state_oldh_33;
  __attribute__ ((aligned (32))) double state_oldh_33[10][130];
  
  double *_tend1dxh_33 = data->tend1dxh_33;
  __attribute__ ((aligned (32))) double tend1dxh_33[10][130];
  
  double *_tend1dyh_33 = data->tend1dyh_33;
  __attribute__ ((aligned (32))) double tend1dyh_33[10][130];
  
  double *_state_oldu_33 = data->state_oldu_33;
  __attribute__ ((aligned (32))) double state_oldu_33[10][130];
  
  double *_tend1dxu_33 = data->tend1dxu_33;
  __attribute__ ((aligned (32))) double tend1dxu_33[10][130];
  
  double *_tend1dyu_33 = data->tend1dyu_33;
  __attribute__ ((aligned (32))) double tend1dyu_33[10][130];
  
  double *_state_oldv_33 = data->state_oldv_33;
  __attribute__ ((aligned (32))) double state_oldv_33[10][130];
  
  double *_tend1dxv_33 = data->tend1dxv_33;
  __attribute__ ((aligned (32))) double tend1dxv_33[10][130];
  
  double *_tend1dyv_33 = data->tend1dyv_33;
  __attribute__ ((aligned (32))) double tend1dyv_33[10][130];
  
  double *_state_newh_33 = data->state_newh_33;
  __attribute__ ((aligned (32))) double state_newh_33[10][130];
  
  double *_state_newu_33 = data->state_newu_33;
  __attribute__ ((aligned (32))) double state_newu_33[10][130];
  
  double *_state_newv_33 = data->state_newv_33;
  __attribute__ ((aligned (32))) double state_newv_33[10][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_33[j][i] = (state_oldh_33[swindex[j]][i] + ((tend1dxh_33[swindex[j]][i] + tend1dyh_33[swindex[j]][i]) * dt));
            state_newu_33[j][i] = (state_oldu_33[swindex[j]][i] + ((tend1dxu_33[swindex[j]][i] + tend1dyu_33[swindex[j]][i]) * dt));
            state_newv_33[j][i] = (state_oldv_33[swindex[j]][i] + ((tend1dxv_33[swindex[j]][i] + tend1dyv_33[swindex[j]][i]) * dt));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

void update_rk2_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_0_info *data = (update_rk2_0_0_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_11 = data->state_oldh_11;
  __attribute__ ((aligned (32))) double state_oldh_11[6][130];
  
  double *_tend1dxh_11 = data->tend1dxh_11;
  __attribute__ ((aligned (32))) double tend1dxh_11[6][130];
  
  double *_tend1dyh_11 = data->tend1dyh_11;
  __attribute__ ((aligned (32))) double tend1dyh_11[6][130];
  
  double *_state_oldu_11 = data->state_oldu_11;
  __attribute__ ((aligned (32))) double state_oldu_11[6][130];
  
  double *_tend1dxu_11 = data->tend1dxu_11;
  __attribute__ ((aligned (32))) double tend1dxu_11[6][130];
  
  double *_tend1dyu_11 = data->tend1dyu_11;
  __attribute__ ((aligned (32))) double tend1dyu_11[6][130];
  
  double *_state_oldv_11 = data->state_oldv_11;
  __attribute__ ((aligned (32))) double state_oldv_11[6][130];
  
  double *_tend1dxv_11 = data->tend1dxv_11;
  __attribute__ ((aligned (32))) double tend1dxv_11[6][130];
  
  double *_tend1dyv_11 = data->tend1dyv_11;
  __attribute__ ((aligned (32))) double tend1dyv_11[6][130];
  
  double *_tend2dxh_11 = data->tend2dxh_11;
  __attribute__ ((aligned (32))) double tend2dxh_11[6][130];
  
  double *_tend2dyh_11 = data->tend2dyh_11;
  __attribute__ ((aligned (32))) double tend2dyh_11[6][130];
  
  double *_tend2dxu_11 = data->tend2dxu_11;
  __attribute__ ((aligned (32))) double tend2dxu_11[6][130];
  
  double *_tend2dyu_11 = data->tend2dyu_11;
  __attribute__ ((aligned (32))) double tend2dyu_11[6][130];
  
  double *_tend2dxv_11 = data->tend2dxv_11;
  __attribute__ ((aligned (32))) double tend2dxv_11[6][130];
  
  double *_tend2dyv_11 = data->tend2dyv_11;
  __attribute__ ((aligned (32))) double tend2dyv_11[6][130];
  
  double *_state_newh_11 = data->state_newh_11;
  __attribute__ ((aligned (32))) double state_newh_11[6][130];
  
  double *_state_newu_11 = data->state_newu_11;
  __attribute__ ((aligned (32))) double state_newu_11[6][130];
  
  double *_state_newv_11 = data->state_newv_11;
  __attribute__ ((aligned (32))) double state_newv_11[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_11+(j*fnumx+ib), &tend2dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_11+(j*fnumx+ib), &tend2dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_11+(j*fnumx+ib), &tend2dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_11+(j*fnumx+ib), &tend2dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_11+(j*fnumx+ib), &tend2dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_11+(j*fnumx+ib), &tend2dyv_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_11+(j*fnumx+ib), &tend2dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_11+(j*fnumx+ib), &tend2dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_11+(j*fnumx+ib), &tend2dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_11+(j*fnumx+ib), &tend2dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_11+(j*fnumx+ib), &tend2dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_11+(j*fnumx+ib), &tend2dyv_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_11[j][i] = (state_oldh_11[swindex[j]][i] + (((((tend1dxh_11[swindex[j]][i] + tend1dyh_11[swindex[j]][i]) + tend2dxh_11[swindex[j]][i]) + tend2dyh_11[swindex[j]][i]) * dt) / 4.0));
            state_newu_11[j][i] = (state_oldu_11[swindex[j]][i] + (((((tend1dxu_11[swindex[j]][i] + tend1dyu_11[swindex[j]][i]) + tend2dxu_11[swindex[j]][i]) + tend2dyu_11[swindex[j]][i]) * dt) / 4.0));
            state_newv_11[j][i] = (state_oldv_11[swindex[j]][i] + (((((tend1dxv_11[swindex[j]][i] + tend1dyv_11[swindex[j]][i]) + tend2dxv_11[swindex[j]][i]) + tend2dyv_11[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_1(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_1_info *data = (update_rk2_0_1_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_12 = data->state_oldh_12;
  __attribute__ ((aligned (32))) double state_oldh_12[6][130];
  
  double *_tend1dxh_12 = data->tend1dxh_12;
  __attribute__ ((aligned (32))) double tend1dxh_12[6][130];
  
  double *_tend1dyh_12 = data->tend1dyh_12;
  __attribute__ ((aligned (32))) double tend1dyh_12[6][130];
  
  double *_state_oldu_12 = data->state_oldu_12;
  __attribute__ ((aligned (32))) double state_oldu_12[6][130];
  
  double *_tend1dxu_12 = data->tend1dxu_12;
  __attribute__ ((aligned (32))) double tend1dxu_12[6][130];
  
  double *_tend1dyu_12 = data->tend1dyu_12;
  __attribute__ ((aligned (32))) double tend1dyu_12[6][130];
  
  double *_state_oldv_12 = data->state_oldv_12;
  __attribute__ ((aligned (32))) double state_oldv_12[6][130];
  
  double *_tend1dxv_12 = data->tend1dxv_12;
  __attribute__ ((aligned (32))) double tend1dxv_12[6][130];
  
  double *_tend1dyv_12 = data->tend1dyv_12;
  __attribute__ ((aligned (32))) double tend1dyv_12[6][130];
  
  double *_tend2dxh_12 = data->tend2dxh_12;
  __attribute__ ((aligned (32))) double tend2dxh_12[6][130];
  
  double *_tend2dyh_12 = data->tend2dyh_12;
  __attribute__ ((aligned (32))) double tend2dyh_12[6][130];
  
  double *_tend2dxu_12 = data->tend2dxu_12;
  __attribute__ ((aligned (32))) double tend2dxu_12[6][130];
  
  double *_tend2dyu_12 = data->tend2dyu_12;
  __attribute__ ((aligned (32))) double tend2dyu_12[6][130];
  
  double *_tend2dxv_12 = data->tend2dxv_12;
  __attribute__ ((aligned (32))) double tend2dxv_12[6][130];
  
  double *_tend2dyv_12 = data->tend2dyv_12;
  __attribute__ ((aligned (32))) double tend2dyv_12[6][130];
  
  double *_state_newh_12 = data->state_newh_12;
  __attribute__ ((aligned (32))) double state_newh_12[6][130];
  
  double *_state_newu_12 = data->state_newu_12;
  __attribute__ ((aligned (32))) double state_newu_12[6][130];
  
  double *_state_newv_12 = data->state_newv_12;
  __attribute__ ((aligned (32))) double state_newv_12[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_12+(j*fnumx+ib), &tend2dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_12+(j*fnumx+ib), &tend2dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_12+(j*fnumx+ib), &tend2dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_12+(j*fnumx+ib), &tend2dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_12+(j*fnumx+ib), &tend2dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_12+(j*fnumx+ib), &tend2dyv_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_12+(j*fnumx+ib), &tend2dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_12+(j*fnumx+ib), &tend2dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_12+(j*fnumx+ib), &tend2dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_12+(j*fnumx+ib), &tend2dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_12+(j*fnumx+ib), &tend2dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_12+(j*fnumx+ib), &tend2dyv_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_12[j][i] = (state_oldh_12[swindex[j]][i] + (((((tend1dxh_12[swindex[j]][i] + tend1dyh_12[swindex[j]][i]) + tend2dxh_12[swindex[j]][i]) + tend2dyh_12[swindex[j]][i]) * dt) / 4.0));
            state_newu_12[j][i] = (state_oldu_12[swindex[j]][i] + (((((tend1dxu_12[swindex[j]][i] + tend1dyu_12[swindex[j]][i]) + tend2dxu_12[swindex[j]][i]) + tend2dyu_12[swindex[j]][i]) * dt) / 4.0));
            state_newv_12[j][i] = (state_oldv_12[swindex[j]][i] + (((((tend1dxv_12[swindex[j]][i] + tend1dyv_12[swindex[j]][i]) + tend2dxv_12[swindex[j]][i]) + tend2dyv_12[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_2(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_2_info *data = (update_rk2_0_2_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_13 = data->state_oldh_13;
  __attribute__ ((aligned (32))) double state_oldh_13[6][130];
  
  double *_tend1dxh_13 = data->tend1dxh_13;
  __attribute__ ((aligned (32))) double tend1dxh_13[6][130];
  
  double *_tend1dyh_13 = data->tend1dyh_13;
  __attribute__ ((aligned (32))) double tend1dyh_13[6][130];
  
  double *_state_oldu_13 = data->state_oldu_13;
  __attribute__ ((aligned (32))) double state_oldu_13[6][130];
  
  double *_tend1dxu_13 = data->tend1dxu_13;
  __attribute__ ((aligned (32))) double tend1dxu_13[6][130];
  
  double *_tend1dyu_13 = data->tend1dyu_13;
  __attribute__ ((aligned (32))) double tend1dyu_13[6][130];
  
  double *_state_oldv_13 = data->state_oldv_13;
  __attribute__ ((aligned (32))) double state_oldv_13[6][130];
  
  double *_tend1dxv_13 = data->tend1dxv_13;
  __attribute__ ((aligned (32))) double tend1dxv_13[6][130];
  
  double *_tend1dyv_13 = data->tend1dyv_13;
  __attribute__ ((aligned (32))) double tend1dyv_13[6][130];
  
  double *_tend2dxh_13 = data->tend2dxh_13;
  __attribute__ ((aligned (32))) double tend2dxh_13[6][130];
  
  double *_tend2dyh_13 = data->tend2dyh_13;
  __attribute__ ((aligned (32))) double tend2dyh_13[6][130];
  
  double *_tend2dxu_13 = data->tend2dxu_13;
  __attribute__ ((aligned (32))) double tend2dxu_13[6][130];
  
  double *_tend2dyu_13 = data->tend2dyu_13;
  __attribute__ ((aligned (32))) double tend2dyu_13[6][130];
  
  double *_tend2dxv_13 = data->tend2dxv_13;
  __attribute__ ((aligned (32))) double tend2dxv_13[6][130];
  
  double *_tend2dyv_13 = data->tend2dyv_13;
  __attribute__ ((aligned (32))) double tend2dyv_13[6][130];
  
  double *_state_newh_13 = data->state_newh_13;
  __attribute__ ((aligned (32))) double state_newh_13[6][130];
  
  double *_state_newu_13 = data->state_newu_13;
  __attribute__ ((aligned (32))) double state_newu_13[6][130];
  
  double *_state_newv_13 = data->state_newv_13;
  __attribute__ ((aligned (32))) double state_newv_13[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_13+(j*fnumx+ib), &tend2dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_13+(j*fnumx+ib), &tend2dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_13+(j*fnumx+ib), &tend2dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_13+(j*fnumx+ib), &tend2dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_13+(j*fnumx+ib), &tend2dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_13+(j*fnumx+ib), &tend2dyv_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_13+(j*fnumx+ib), &tend2dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_13+(j*fnumx+ib), &tend2dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_13+(j*fnumx+ib), &tend2dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_13+(j*fnumx+ib), &tend2dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_13+(j*fnumx+ib), &tend2dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_13+(j*fnumx+ib), &tend2dyv_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_13[j][i] = (state_oldh_13[swindex[j]][i] + (((((tend1dxh_13[swindex[j]][i] + tend1dyh_13[swindex[j]][i]) + tend2dxh_13[swindex[j]][i]) + tend2dyh_13[swindex[j]][i]) * dt) / 4.0));
            state_newu_13[j][i] = (state_oldu_13[swindex[j]][i] + (((((tend1dxu_13[swindex[j]][i] + tend1dyu_13[swindex[j]][i]) + tend2dxu_13[swindex[j]][i]) + tend2dyu_13[swindex[j]][i]) * dt) / 4.0));
            state_newv_13[j][i] = (state_oldv_13[swindex[j]][i] + (((((tend1dxv_13[swindex[j]][i] + tend1dyv_13[swindex[j]][i]) + tend2dxv_13[swindex[j]][i]) + tend2dyv_13[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_3(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_3_info *data = (update_rk2_0_3_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_21 = data->state_oldh_21;
  __attribute__ ((aligned (32))) double state_oldh_21[6][130];
  
  double *_tend1dxh_21 = data->tend1dxh_21;
  __attribute__ ((aligned (32))) double tend1dxh_21[6][130];
  
  double *_tend1dyh_21 = data->tend1dyh_21;
  __attribute__ ((aligned (32))) double tend1dyh_21[6][130];
  
  double *_state_oldu_21 = data->state_oldu_21;
  __attribute__ ((aligned (32))) double state_oldu_21[6][130];
  
  double *_tend1dxu_21 = data->tend1dxu_21;
  __attribute__ ((aligned (32))) double tend1dxu_21[6][130];
  
  double *_tend1dyu_21 = data->tend1dyu_21;
  __attribute__ ((aligned (32))) double tend1dyu_21[6][130];
  
  double *_state_oldv_21 = data->state_oldv_21;
  __attribute__ ((aligned (32))) double state_oldv_21[6][130];
  
  double *_tend1dxv_21 = data->tend1dxv_21;
  __attribute__ ((aligned (32))) double tend1dxv_21[6][130];
  
  double *_tend1dyv_21 = data->tend1dyv_21;
  __attribute__ ((aligned (32))) double tend1dyv_21[6][130];
  
  double *_tend2dxh_21 = data->tend2dxh_21;
  __attribute__ ((aligned (32))) double tend2dxh_21[6][130];
  
  double *_tend2dyh_21 = data->tend2dyh_21;
  __attribute__ ((aligned (32))) double tend2dyh_21[6][130];
  
  double *_tend2dxu_21 = data->tend2dxu_21;
  __attribute__ ((aligned (32))) double tend2dxu_21[6][130];
  
  double *_tend2dyu_21 = data->tend2dyu_21;
  __attribute__ ((aligned (32))) double tend2dyu_21[6][130];
  
  double *_tend2dxv_21 = data->tend2dxv_21;
  __attribute__ ((aligned (32))) double tend2dxv_21[6][130];
  
  double *_tend2dyv_21 = data->tend2dyv_21;
  __attribute__ ((aligned (32))) double tend2dyv_21[6][130];
  
  double *_state_newh_21 = data->state_newh_21;
  __attribute__ ((aligned (32))) double state_newh_21[6][130];
  
  double *_state_newu_21 = data->state_newu_21;
  __attribute__ ((aligned (32))) double state_newu_21[6][130];
  
  double *_state_newv_21 = data->state_newv_21;
  __attribute__ ((aligned (32))) double state_newv_21[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_21+(j*fnumx+ib), &tend2dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_21+(j*fnumx+ib), &tend2dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_21+(j*fnumx+ib), &tend2dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_21+(j*fnumx+ib), &tend2dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_21+(j*fnumx+ib), &tend2dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_21+(j*fnumx+ib), &tend2dyv_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_21+(j*fnumx+ib), &tend2dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_21+(j*fnumx+ib), &tend2dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_21+(j*fnumx+ib), &tend2dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_21+(j*fnumx+ib), &tend2dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_21+(j*fnumx+ib), &tend2dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_21+(j*fnumx+ib), &tend2dyv_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_21[j][i] = (state_oldh_21[swindex[j]][i] + (((((tend1dxh_21[swindex[j]][i] + tend1dyh_21[swindex[j]][i]) + tend2dxh_21[swindex[j]][i]) + tend2dyh_21[swindex[j]][i]) * dt) / 4.0));
            state_newu_21[j][i] = (state_oldu_21[swindex[j]][i] + (((((tend1dxu_21[swindex[j]][i] + tend1dyu_21[swindex[j]][i]) + tend2dxu_21[swindex[j]][i]) + tend2dyu_21[swindex[j]][i]) * dt) / 4.0));
            state_newv_21[j][i] = (state_oldv_21[swindex[j]][i] + (((((tend1dxv_21[swindex[j]][i] + tend1dyv_21[swindex[j]][i]) + tend2dxv_21[swindex[j]][i]) + tend2dyv_21[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_4(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_4_info *data = (update_rk2_0_4_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_22 = data->state_oldh_22;
  __attribute__ ((aligned (32))) double state_oldh_22[6][130];
  
  double *_tend1dxh_22 = data->tend1dxh_22;
  __attribute__ ((aligned (32))) double tend1dxh_22[6][130];
  
  double *_tend1dyh_22 = data->tend1dyh_22;
  __attribute__ ((aligned (32))) double tend1dyh_22[6][130];
  
  double *_state_oldu_22 = data->state_oldu_22;
  __attribute__ ((aligned (32))) double state_oldu_22[6][130];
  
  double *_tend1dxu_22 = data->tend1dxu_22;
  __attribute__ ((aligned (32))) double tend1dxu_22[6][130];
  
  double *_tend1dyu_22 = data->tend1dyu_22;
  __attribute__ ((aligned (32))) double tend1dyu_22[6][130];
  
  double *_state_oldv_22 = data->state_oldv_22;
  __attribute__ ((aligned (32))) double state_oldv_22[6][130];
  
  double *_tend1dxv_22 = data->tend1dxv_22;
  __attribute__ ((aligned (32))) double tend1dxv_22[6][130];
  
  double *_tend1dyv_22 = data->tend1dyv_22;
  __attribute__ ((aligned (32))) double tend1dyv_22[6][130];
  
  double *_tend2dxh_22 = data->tend2dxh_22;
  __attribute__ ((aligned (32))) double tend2dxh_22[6][130];
  
  double *_tend2dyh_22 = data->tend2dyh_22;
  __attribute__ ((aligned (32))) double tend2dyh_22[6][130];
  
  double *_tend2dxu_22 = data->tend2dxu_22;
  __attribute__ ((aligned (32))) double tend2dxu_22[6][130];
  
  double *_tend2dyu_22 = data->tend2dyu_22;
  __attribute__ ((aligned (32))) double tend2dyu_22[6][130];
  
  double *_tend2dxv_22 = data->tend2dxv_22;
  __attribute__ ((aligned (32))) double tend2dxv_22[6][130];
  
  double *_tend2dyv_22 = data->tend2dyv_22;
  __attribute__ ((aligned (32))) double tend2dyv_22[6][130];
  
  double *_state_newh_22 = data->state_newh_22;
  __attribute__ ((aligned (32))) double state_newh_22[6][130];
  
  double *_state_newu_22 = data->state_newu_22;
  __attribute__ ((aligned (32))) double state_newu_22[6][130];
  
  double *_state_newv_22 = data->state_newv_22;
  __attribute__ ((aligned (32))) double state_newv_22[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_22+(j*fnumx+ib), &tend2dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_22+(j*fnumx+ib), &tend2dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_22+(j*fnumx+ib), &tend2dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_22+(j*fnumx+ib), &tend2dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_22+(j*fnumx+ib), &tend2dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_22+(j*fnumx+ib), &tend2dyv_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_22+(j*fnumx+ib), &tend2dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_22+(j*fnumx+ib), &tend2dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_22+(j*fnumx+ib), &tend2dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_22+(j*fnumx+ib), &tend2dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_22+(j*fnumx+ib), &tend2dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_22+(j*fnumx+ib), &tend2dyv_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_22[j][i] = (state_oldh_22[swindex[j]][i] + (((((tend1dxh_22[swindex[j]][i] + tend1dyh_22[swindex[j]][i]) + tend2dxh_22[swindex[j]][i]) + tend2dyh_22[swindex[j]][i]) * dt) / 4.0));
            state_newu_22[j][i] = (state_oldu_22[swindex[j]][i] + (((((tend1dxu_22[swindex[j]][i] + tend1dyu_22[swindex[j]][i]) + tend2dxu_22[swindex[j]][i]) + tend2dyu_22[swindex[j]][i]) * dt) / 4.0));
            state_newv_22[j][i] = (state_oldv_22[swindex[j]][i] + (((((tend1dxv_22[swindex[j]][i] + tend1dyv_22[swindex[j]][i]) + tend2dxv_22[swindex[j]][i]) + tend2dyv_22[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_5(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_5_info *data = (update_rk2_0_5_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_23 = data->state_oldh_23;
  __attribute__ ((aligned (32))) double state_oldh_23[6][130];
  
  double *_tend1dxh_23 = data->tend1dxh_23;
  __attribute__ ((aligned (32))) double tend1dxh_23[6][130];
  
  double *_tend1dyh_23 = data->tend1dyh_23;
  __attribute__ ((aligned (32))) double tend1dyh_23[6][130];
  
  double *_state_oldu_23 = data->state_oldu_23;
  __attribute__ ((aligned (32))) double state_oldu_23[6][130];
  
  double *_tend1dxu_23 = data->tend1dxu_23;
  __attribute__ ((aligned (32))) double tend1dxu_23[6][130];
  
  double *_tend1dyu_23 = data->tend1dyu_23;
  __attribute__ ((aligned (32))) double tend1dyu_23[6][130];
  
  double *_state_oldv_23 = data->state_oldv_23;
  __attribute__ ((aligned (32))) double state_oldv_23[6][130];
  
  double *_tend1dxv_23 = data->tend1dxv_23;
  __attribute__ ((aligned (32))) double tend1dxv_23[6][130];
  
  double *_tend1dyv_23 = data->tend1dyv_23;
  __attribute__ ((aligned (32))) double tend1dyv_23[6][130];
  
  double *_tend2dxh_23 = data->tend2dxh_23;
  __attribute__ ((aligned (32))) double tend2dxh_23[6][130];
  
  double *_tend2dyh_23 = data->tend2dyh_23;
  __attribute__ ((aligned (32))) double tend2dyh_23[6][130];
  
  double *_tend2dxu_23 = data->tend2dxu_23;
  __attribute__ ((aligned (32))) double tend2dxu_23[6][130];
  
  double *_tend2dyu_23 = data->tend2dyu_23;
  __attribute__ ((aligned (32))) double tend2dyu_23[6][130];
  
  double *_tend2dxv_23 = data->tend2dxv_23;
  __attribute__ ((aligned (32))) double tend2dxv_23[6][130];
  
  double *_tend2dyv_23 = data->tend2dyv_23;
  __attribute__ ((aligned (32))) double tend2dyv_23[6][130];
  
  double *_state_newh_23 = data->state_newh_23;
  __attribute__ ((aligned (32))) double state_newh_23[6][130];
  
  double *_state_newu_23 = data->state_newu_23;
  __attribute__ ((aligned (32))) double state_newu_23[6][130];
  
  double *_state_newv_23 = data->state_newv_23;
  __attribute__ ((aligned (32))) double state_newv_23[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_23+(j*fnumx+ib), &tend2dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_23+(j*fnumx+ib), &tend2dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_23+(j*fnumx+ib), &tend2dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_23+(j*fnumx+ib), &tend2dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_23+(j*fnumx+ib), &tend2dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_23+(j*fnumx+ib), &tend2dyv_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_23+(j*fnumx+ib), &tend2dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_23+(j*fnumx+ib), &tend2dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_23+(j*fnumx+ib), &tend2dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_23+(j*fnumx+ib), &tend2dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_23+(j*fnumx+ib), &tend2dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_23+(j*fnumx+ib), &tend2dyv_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_23[j][i] = (state_oldh_23[swindex[j]][i] + (((((tend1dxh_23[swindex[j]][i] + tend1dyh_23[swindex[j]][i]) + tend2dxh_23[swindex[j]][i]) + tend2dyh_23[swindex[j]][i]) * dt) / 4.0));
            state_newu_23[j][i] = (state_oldu_23[swindex[j]][i] + (((((tend1dxu_23[swindex[j]][i] + tend1dyu_23[swindex[j]][i]) + tend2dxu_23[swindex[j]][i]) + tend2dyu_23[swindex[j]][i]) * dt) / 4.0));
            state_newv_23[j][i] = (state_oldv_23[swindex[j]][i] + (((((tend1dxv_23[swindex[j]][i] + tend1dyv_23[swindex[j]][i]) + tend2dxv_23[swindex[j]][i]) + tend2dyv_23[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_6(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_6_info *data = (update_rk2_0_6_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_31 = data->state_oldh_31;
  __attribute__ ((aligned (32))) double state_oldh_31[6][130];
  
  double *_tend1dxh_31 = data->tend1dxh_31;
  __attribute__ ((aligned (32))) double tend1dxh_31[6][130];
  
  double *_tend1dyh_31 = data->tend1dyh_31;
  __attribute__ ((aligned (32))) double tend1dyh_31[6][130];
  
  double *_state_oldu_31 = data->state_oldu_31;
  __attribute__ ((aligned (32))) double state_oldu_31[6][130];
  
  double *_tend1dxu_31 = data->tend1dxu_31;
  __attribute__ ((aligned (32))) double tend1dxu_31[6][130];
  
  double *_tend1dyu_31 = data->tend1dyu_31;
  __attribute__ ((aligned (32))) double tend1dyu_31[6][130];
  
  double *_state_oldv_31 = data->state_oldv_31;
  __attribute__ ((aligned (32))) double state_oldv_31[6][130];
  
  double *_tend1dxv_31 = data->tend1dxv_31;
  __attribute__ ((aligned (32))) double tend1dxv_31[6][130];
  
  double *_tend1dyv_31 = data->tend1dyv_31;
  __attribute__ ((aligned (32))) double tend1dyv_31[6][130];
  
  double *_tend2dxh_31 = data->tend2dxh_31;
  __attribute__ ((aligned (32))) double tend2dxh_31[6][130];
  
  double *_tend2dyh_31 = data->tend2dyh_31;
  __attribute__ ((aligned (32))) double tend2dyh_31[6][130];
  
  double *_tend2dxu_31 = data->tend2dxu_31;
  __attribute__ ((aligned (32))) double tend2dxu_31[6][130];
  
  double *_tend2dyu_31 = data->tend2dyu_31;
  __attribute__ ((aligned (32))) double tend2dyu_31[6][130];
  
  double *_tend2dxv_31 = data->tend2dxv_31;
  __attribute__ ((aligned (32))) double tend2dxv_31[6][130];
  
  double *_tend2dyv_31 = data->tend2dyv_31;
  __attribute__ ((aligned (32))) double tend2dyv_31[6][130];
  
  double *_state_newh_31 = data->state_newh_31;
  __attribute__ ((aligned (32))) double state_newh_31[6][130];
  
  double *_state_newu_31 = data->state_newu_31;
  __attribute__ ((aligned (32))) double state_newu_31[6][130];
  
  double *_state_newv_31 = data->state_newv_31;
  __attribute__ ((aligned (32))) double state_newv_31[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_31+(j*fnumx+ib), &tend2dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_31+(j*fnumx+ib), &tend2dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_31+(j*fnumx+ib), &tend2dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_31+(j*fnumx+ib), &tend2dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_31+(j*fnumx+ib), &tend2dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_31+(j*fnumx+ib), &tend2dyv_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_31+(j*fnumx+ib), &tend2dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_31+(j*fnumx+ib), &tend2dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_31+(j*fnumx+ib), &tend2dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_31+(j*fnumx+ib), &tend2dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_31+(j*fnumx+ib), &tend2dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_31+(j*fnumx+ib), &tend2dyv_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_31[j][i] = (state_oldh_31[swindex[j]][i] + (((((tend1dxh_31[swindex[j]][i] + tend1dyh_31[swindex[j]][i]) + tend2dxh_31[swindex[j]][i]) + tend2dyh_31[swindex[j]][i]) * dt) / 4.0));
            state_newu_31[j][i] = (state_oldu_31[swindex[j]][i] + (((((tend1dxu_31[swindex[j]][i] + tend1dyu_31[swindex[j]][i]) + tend2dxu_31[swindex[j]][i]) + tend2dyu_31[swindex[j]][i]) * dt) / 4.0));
            state_newv_31[j][i] = (state_oldv_31[swindex[j]][i] + (((((tend1dxv_31[swindex[j]][i] + tend1dyv_31[swindex[j]][i]) + tend2dxv_31[swindex[j]][i]) + tend2dyv_31[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_7(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_7_info *data = (update_rk2_0_7_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_32 = data->state_oldh_32;
  __attribute__ ((aligned (32))) double state_oldh_32[6][130];
  
  double *_tend1dxh_32 = data->tend1dxh_32;
  __attribute__ ((aligned (32))) double tend1dxh_32[6][130];
  
  double *_tend1dyh_32 = data->tend1dyh_32;
  __attribute__ ((aligned (32))) double tend1dyh_32[6][130];
  
  double *_state_oldu_32 = data->state_oldu_32;
  __attribute__ ((aligned (32))) double state_oldu_32[6][130];
  
  double *_tend1dxu_32 = data->tend1dxu_32;
  __attribute__ ((aligned (32))) double tend1dxu_32[6][130];
  
  double *_tend1dyu_32 = data->tend1dyu_32;
  __attribute__ ((aligned (32))) double tend1dyu_32[6][130];
  
  double *_state_oldv_32 = data->state_oldv_32;
  __attribute__ ((aligned (32))) double state_oldv_32[6][130];
  
  double *_tend1dxv_32 = data->tend1dxv_32;
  __attribute__ ((aligned (32))) double tend1dxv_32[6][130];
  
  double *_tend1dyv_32 = data->tend1dyv_32;
  __attribute__ ((aligned (32))) double tend1dyv_32[6][130];
  
  double *_tend2dxh_32 = data->tend2dxh_32;
  __attribute__ ((aligned (32))) double tend2dxh_32[6][130];
  
  double *_tend2dyh_32 = data->tend2dyh_32;
  __attribute__ ((aligned (32))) double tend2dyh_32[6][130];
  
  double *_tend2dxu_32 = data->tend2dxu_32;
  __attribute__ ((aligned (32))) double tend2dxu_32[6][130];
  
  double *_tend2dyu_32 = data->tend2dyu_32;
  __attribute__ ((aligned (32))) double tend2dyu_32[6][130];
  
  double *_tend2dxv_32 = data->tend2dxv_32;
  __attribute__ ((aligned (32))) double tend2dxv_32[6][130];
  
  double *_tend2dyv_32 = data->tend2dyv_32;
  __attribute__ ((aligned (32))) double tend2dyv_32[6][130];
  
  double *_state_newh_32 = data->state_newh_32;
  __attribute__ ((aligned (32))) double state_newh_32[6][130];
  
  double *_state_newu_32 = data->state_newu_32;
  __attribute__ ((aligned (32))) double state_newu_32[6][130];
  
  double *_state_newv_32 = data->state_newv_32;
  __attribute__ ((aligned (32))) double state_newv_32[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_32+(j*fnumx+ib), &tend2dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_32+(j*fnumx+ib), &tend2dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_32+(j*fnumx+ib), &tend2dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_32+(j*fnumx+ib), &tend2dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_32+(j*fnumx+ib), &tend2dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_32+(j*fnumx+ib), &tend2dyv_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_32+(j*fnumx+ib), &tend2dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_32+(j*fnumx+ib), &tend2dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_32+(j*fnumx+ib), &tend2dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_32+(j*fnumx+ib), &tend2dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_32+(j*fnumx+ib), &tend2dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_32+(j*fnumx+ib), &tend2dyv_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_32[j][i] = (state_oldh_32[swindex[j]][i] + (((((tend1dxh_32[swindex[j]][i] + tend1dyh_32[swindex[j]][i]) + tend2dxh_32[swindex[j]][i]) + tend2dyh_32[swindex[j]][i]) * dt) / 4.0));
            state_newu_32[j][i] = (state_oldu_32[swindex[j]][i] + (((((tend1dxu_32[swindex[j]][i] + tend1dyu_32[swindex[j]][i]) + tend2dxu_32[swindex[j]][i]) + tend2dyu_32[swindex[j]][i]) * dt) / 4.0));
            state_newv_32[j][i] = (state_oldv_32[swindex[j]][i] + (((((tend1dxv_32[swindex[j]][i] + tend1dyv_32[swindex[j]][i]) + tend2dxv_32[swindex[j]][i]) + tend2dyv_32[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk2_0_8(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk2_0_8_info *data = (update_rk2_0_8_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_33 = data->state_oldh_33;
  __attribute__ ((aligned (32))) double state_oldh_33[6][130];
  
  double *_tend1dxh_33 = data->tend1dxh_33;
  __attribute__ ((aligned (32))) double tend1dxh_33[6][130];
  
  double *_tend1dyh_33 = data->tend1dyh_33;
  __attribute__ ((aligned (32))) double tend1dyh_33[6][130];
  
  double *_state_oldu_33 = data->state_oldu_33;
  __attribute__ ((aligned (32))) double state_oldu_33[6][130];
  
  double *_tend1dxu_33 = data->tend1dxu_33;
  __attribute__ ((aligned (32))) double tend1dxu_33[6][130];
  
  double *_tend1dyu_33 = data->tend1dyu_33;
  __attribute__ ((aligned (32))) double tend1dyu_33[6][130];
  
  double *_state_oldv_33 = data->state_oldv_33;
  __attribute__ ((aligned (32))) double state_oldv_33[6][130];
  
  double *_tend1dxv_33 = data->tend1dxv_33;
  __attribute__ ((aligned (32))) double tend1dxv_33[6][130];
  
  double *_tend1dyv_33 = data->tend1dyv_33;
  __attribute__ ((aligned (32))) double tend1dyv_33[6][130];
  
  double *_tend2dxh_33 = data->tend2dxh_33;
  __attribute__ ((aligned (32))) double tend2dxh_33[6][130];
  
  double *_tend2dyh_33 = data->tend2dyh_33;
  __attribute__ ((aligned (32))) double tend2dyh_33[6][130];
  
  double *_tend2dxu_33 = data->tend2dxu_33;
  __attribute__ ((aligned (32))) double tend2dxu_33[6][130];
  
  double *_tend2dyu_33 = data->tend2dyu_33;
  __attribute__ ((aligned (32))) double tend2dyu_33[6][130];
  
  double *_tend2dxv_33 = data->tend2dxv_33;
  __attribute__ ((aligned (32))) double tend2dxv_33[6][130];
  
  double *_tend2dyv_33 = data->tend2dyv_33;
  __attribute__ ((aligned (32))) double tend2dyv_33[6][130];
  
  double *_state_newh_33 = data->state_newh_33;
  __attribute__ ((aligned (32))) double state_newh_33[6][130];
  
  double *_state_newu_33 = data->state_newu_33;
  __attribute__ ((aligned (32))) double state_newu_33[6][130];
  
  double *_state_newv_33 = data->state_newv_33;
  __attribute__ ((aligned (32))) double state_newv_33[6][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_33+(j*fnumx+ib), &tend2dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_33+(j*fnumx+ib), &tend2dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_33+(j*fnumx+ib), &tend2dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_33+(j*fnumx+ib), &tend2dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_33+(j*fnumx+ib), &tend2dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_33+(j*fnumx+ib), &tend2dyv_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_33+(j*fnumx+ib), &tend2dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_33+(j*fnumx+ib), &tend2dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_33+(j*fnumx+ib), &tend2dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_33+(j*fnumx+ib), &tend2dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_33+(j*fnumx+ib), &tend2dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_33+(j*fnumx+ib), &tend2dyv_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_33[j][i] = (state_oldh_33[swindex[j]][i] + (((((tend1dxh_33[swindex[j]][i] + tend1dyh_33[swindex[j]][i]) + tend2dxh_33[swindex[j]][i]) + tend2dyh_33[swindex[j]][i]) * dt) / 4.0));
            state_newu_33[j][i] = (state_oldu_33[swindex[j]][i] + (((((tend1dxu_33[swindex[j]][i] + tend1dyu_33[swindex[j]][i]) + tend2dxu_33[swindex[j]][i]) + tend2dyu_33[swindex[j]][i]) * dt) / 4.0));
            state_newv_33[j][i] = (state_oldv_33[swindex[j]][i] + (((((tend1dxv_33[swindex[j]][i] + tend1dyv_33[swindex[j]][i]) + tend2dxv_33[swindex[j]][i]) + tend2dyv_33[swindex[j]][i]) * dt) / 4.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

void update_rk3_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_0_info *data = (update_rk3_0_0_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_11 = data->state_oldh_11;
  __attribute__ ((aligned (32))) double state_oldh_11[4][130];
  
  double *_tend1dxh_11 = data->tend1dxh_11;
  __attribute__ ((aligned (32))) double tend1dxh_11[4][130];
  
  double *_tend1dyh_11 = data->tend1dyh_11;
  __attribute__ ((aligned (32))) double tend1dyh_11[4][130];
  
  double *_state_oldu_11 = data->state_oldu_11;
  __attribute__ ((aligned (32))) double state_oldu_11[4][130];
  
  double *_tend1dxu_11 = data->tend1dxu_11;
  __attribute__ ((aligned (32))) double tend1dxu_11[4][130];
  
  double *_tend1dyu_11 = data->tend1dyu_11;
  __attribute__ ((aligned (32))) double tend1dyu_11[4][130];
  
  double *_state_oldv_11 = data->state_oldv_11;
  __attribute__ ((aligned (32))) double state_oldv_11[4][130];
  
  double *_tend1dxv_11 = data->tend1dxv_11;
  __attribute__ ((aligned (32))) double tend1dxv_11[4][130];
  
  double *_tend1dyv_11 = data->tend1dyv_11;
  __attribute__ ((aligned (32))) double tend1dyv_11[4][130];
  
  double *_tend2dxh_11 = data->tend2dxh_11;
  __attribute__ ((aligned (32))) double tend2dxh_11[4][130];
  
  double *_tend2dyh_11 = data->tend2dyh_11;
  __attribute__ ((aligned (32))) double tend2dyh_11[4][130];
  
  double *_tend3dxh_11 = data->tend3dxh_11;
  __attribute__ ((aligned (32))) double tend3dxh_11[4][130];
  
  double *_tend3dyh_11 = data->tend3dyh_11;
  __attribute__ ((aligned (32))) double tend3dyh_11[4][130];
  
  double *_tend2dxu_11 = data->tend2dxu_11;
  __attribute__ ((aligned (32))) double tend2dxu_11[4][130];
  
  double *_tend2dyu_11 = data->tend2dyu_11;
  __attribute__ ((aligned (32))) double tend2dyu_11[4][130];
  
  double *_tend3dxu_11 = data->tend3dxu_11;
  __attribute__ ((aligned (32))) double tend3dxu_11[4][130];
  
  double *_tend3dyu_11 = data->tend3dyu_11;
  __attribute__ ((aligned (32))) double tend3dyu_11[4][130];
  
  double *_tend2dxv_11 = data->tend2dxv_11;
  __attribute__ ((aligned (32))) double tend2dxv_11[4][130];
  
  double *_tend2dyv_11 = data->tend2dyv_11;
  __attribute__ ((aligned (32))) double tend2dyv_11[4][130];
  
  double *_tend3dxv_11 = data->tend3dxv_11;
  __attribute__ ((aligned (32))) double tend3dxv_11[4][130];
  
  double *_tend3dyv_11 = data->tend3dyv_11;
  __attribute__ ((aligned (32))) double tend3dyv_11[4][130];
  
  double *_state_newh_11 = data->state_newh_11;
  __attribute__ ((aligned (32))) double state_newh_11[4][130];
  
  double *_state_newu_11 = data->state_newu_11;
  __attribute__ ((aligned (32))) double state_newu_11[4][130];
  
  double *_state_newv_11 = data->state_newv_11;
  __attribute__ ((aligned (32))) double state_newv_11[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_11+(j*fnumx+ib), &tend2dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_11+(j*fnumx+ib), &tend2dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_11+(j*fnumx+ib), &tend3dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_11+(j*fnumx+ib), &tend3dxh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_11+(j*fnumx+ib), &tend3dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_11+(j*fnumx+ib), &tend3dyh_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_11+(j*fnumx+ib), &tend2dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_11+(j*fnumx+ib), &tend2dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_11+(j*fnumx+ib), &tend3dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_11+(j*fnumx+ib), &tend3dxu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_11+(j*fnumx+ib), &tend3dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_11+(j*fnumx+ib), &tend3dyu_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_11+(j*fnumx+ib), &tend2dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_11+(j*fnumx+ib), &tend2dyv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_11+(j*fnumx+ib), &tend3dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_11+(j*fnumx+ib), &tend3dxv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_11+(j*fnumx+ib), &tend3dyv_11[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_11+(j*fnumx+ib), &tend3dyv_11[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_11+(j*fnumx+ib), &state_oldh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_11+(j*fnumx+ib), &tend1dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_11+(j*fnumx+ib), &tend1dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_11+(j*fnumx+ib), &state_oldu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_11+(j*fnumx+ib), &tend1dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_11+(j*fnumx+ib), &tend1dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_11+(j*fnumx+ib), &state_oldv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_11+(j*fnumx+ib), &tend1dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_11+(j*fnumx+ib), &tend1dyv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_11+(j*fnumx+ib), &tend2dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_11+(j*fnumx+ib), &tend2dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_11+(j*fnumx+ib), &tend3dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_11+(j*fnumx+ib), &tend3dxh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_11+(j*fnumx+ib), &tend3dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_11+(j*fnumx+ib), &tend3dyh_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_11+(j*fnumx+ib), &tend2dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_11+(j*fnumx+ib), &tend2dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_11+(j*fnumx+ib), &tend3dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_11+(j*fnumx+ib), &tend3dxu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_11+(j*fnumx+ib), &tend3dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_11+(j*fnumx+ib), &tend3dyu_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_11+(j*fnumx+ib), &tend2dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_11+(j*fnumx+ib), &tend2dyv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_11+(j*fnumx+ib), &tend3dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_11+(j*fnumx+ib), &tend3dxv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_11+(j*fnumx+ib), &tend3dyv_11[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_11+(j*fnumx+ib), &tend3dyv_11[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_11[j][i] = (state_oldh_11[swindex[j]][i] + (((((((tend1dxh_11[swindex[j]][i] + tend1dyh_11[swindex[j]][i]) + tend2dxh_11[swindex[j]][i]) + tend2dyh_11[swindex[j]][i]) + (4.0 * tend3dxh_11[swindex[j]][i])) + (4.0 * tend3dyh_11[swindex[j]][i])) * dt) / 6.0));
            state_newu_11[j][i] = (state_oldu_11[swindex[j]][i] + (((((((tend1dxu_11[swindex[j]][i] + tend1dyu_11[swindex[j]][i]) + tend2dxu_11[swindex[j]][i]) + tend2dyu_11[swindex[j]][i]) + (4.0 * tend3dxu_11[swindex[j]][i])) + (4.0 * tend3dyu_11[swindex[j]][i])) * dt) / 6.0));
            state_newv_11[j][i] = (state_oldv_11[swindex[j]][i] + (((((((tend1dxv_11[swindex[j]][i] + tend1dyv_11[swindex[j]][i]) + tend2dxv_11[swindex[j]][i]) + tend2dyv_11[swindex[j]][i]) + (4.0 * tend3dxv_11[swindex[j]][i])) + (4.0 * tend3dyv_11[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_11+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_11[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_1(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_1_info *data = (update_rk3_0_1_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_12 = data->state_oldh_12;
  __attribute__ ((aligned (32))) double state_oldh_12[4][130];
  
  double *_tend1dxh_12 = data->tend1dxh_12;
  __attribute__ ((aligned (32))) double tend1dxh_12[4][130];
  
  double *_tend1dyh_12 = data->tend1dyh_12;
  __attribute__ ((aligned (32))) double tend1dyh_12[4][130];
  
  double *_state_oldu_12 = data->state_oldu_12;
  __attribute__ ((aligned (32))) double state_oldu_12[4][130];
  
  double *_tend1dxu_12 = data->tend1dxu_12;
  __attribute__ ((aligned (32))) double tend1dxu_12[4][130];
  
  double *_tend1dyu_12 = data->tend1dyu_12;
  __attribute__ ((aligned (32))) double tend1dyu_12[4][130];
  
  double *_state_oldv_12 = data->state_oldv_12;
  __attribute__ ((aligned (32))) double state_oldv_12[4][130];
  
  double *_tend1dxv_12 = data->tend1dxv_12;
  __attribute__ ((aligned (32))) double tend1dxv_12[4][130];
  
  double *_tend1dyv_12 = data->tend1dyv_12;
  __attribute__ ((aligned (32))) double tend1dyv_12[4][130];
  
  double *_tend2dxh_12 = data->tend2dxh_12;
  __attribute__ ((aligned (32))) double tend2dxh_12[4][130];
  
  double *_tend2dyh_12 = data->tend2dyh_12;
  __attribute__ ((aligned (32))) double tend2dyh_12[4][130];
  
  double *_tend3dxh_12 = data->tend3dxh_12;
  __attribute__ ((aligned (32))) double tend3dxh_12[4][130];
  
  double *_tend3dyh_12 = data->tend3dyh_12;
  __attribute__ ((aligned (32))) double tend3dyh_12[4][130];
  
  double *_tend2dxu_12 = data->tend2dxu_12;
  __attribute__ ((aligned (32))) double tend2dxu_12[4][130];
  
  double *_tend2dyu_12 = data->tend2dyu_12;
  __attribute__ ((aligned (32))) double tend2dyu_12[4][130];
  
  double *_tend3dxu_12 = data->tend3dxu_12;
  __attribute__ ((aligned (32))) double tend3dxu_12[4][130];
  
  double *_tend3dyu_12 = data->tend3dyu_12;
  __attribute__ ((aligned (32))) double tend3dyu_12[4][130];
  
  double *_tend2dxv_12 = data->tend2dxv_12;
  __attribute__ ((aligned (32))) double tend2dxv_12[4][130];
  
  double *_tend2dyv_12 = data->tend2dyv_12;
  __attribute__ ((aligned (32))) double tend2dyv_12[4][130];
  
  double *_tend3dxv_12 = data->tend3dxv_12;
  __attribute__ ((aligned (32))) double tend3dxv_12[4][130];
  
  double *_tend3dyv_12 = data->tend3dyv_12;
  __attribute__ ((aligned (32))) double tend3dyv_12[4][130];
  
  double *_state_newh_12 = data->state_newh_12;
  __attribute__ ((aligned (32))) double state_newh_12[4][130];
  
  double *_state_newu_12 = data->state_newu_12;
  __attribute__ ((aligned (32))) double state_newu_12[4][130];
  
  double *_state_newv_12 = data->state_newv_12;
  __attribute__ ((aligned (32))) double state_newv_12[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_12+(j*fnumx+ib), &tend2dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_12+(j*fnumx+ib), &tend2dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_12+(j*fnumx+ib), &tend3dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_12+(j*fnumx+ib), &tend3dxh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_12+(j*fnumx+ib), &tend3dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_12+(j*fnumx+ib), &tend3dyh_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_12+(j*fnumx+ib), &tend2dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_12+(j*fnumx+ib), &tend2dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_12+(j*fnumx+ib), &tend3dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_12+(j*fnumx+ib), &tend3dxu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_12+(j*fnumx+ib), &tend3dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_12+(j*fnumx+ib), &tend3dyu_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_12+(j*fnumx+ib), &tend2dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_12+(j*fnumx+ib), &tend2dyv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_12+(j*fnumx+ib), &tend3dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_12+(j*fnumx+ib), &tend3dxv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_12+(j*fnumx+ib), &tend3dyv_12[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_12+(j*fnumx+ib), &tend3dyv_12[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_12+(j*fnumx+ib), &state_oldh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_12+(j*fnumx+ib), &tend1dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_12+(j*fnumx+ib), &tend1dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_12+(j*fnumx+ib), &state_oldu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_12+(j*fnumx+ib), &tend1dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_12+(j*fnumx+ib), &tend1dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_12+(j*fnumx+ib), &state_oldv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_12+(j*fnumx+ib), &tend1dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_12+(j*fnumx+ib), &tend1dyv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_12+(j*fnumx+ib), &tend2dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_12+(j*fnumx+ib), &tend2dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_12+(j*fnumx+ib), &tend3dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_12+(j*fnumx+ib), &tend3dxh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_12+(j*fnumx+ib), &tend3dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_12+(j*fnumx+ib), &tend3dyh_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_12+(j*fnumx+ib), &tend2dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_12+(j*fnumx+ib), &tend2dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_12+(j*fnumx+ib), &tend3dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_12+(j*fnumx+ib), &tend3dxu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_12+(j*fnumx+ib), &tend3dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_12+(j*fnumx+ib), &tend3dyu_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_12+(j*fnumx+ib), &tend2dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_12+(j*fnumx+ib), &tend2dyv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_12+(j*fnumx+ib), &tend3dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_12+(j*fnumx+ib), &tend3dxv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_12+(j*fnumx+ib), &tend3dyv_12[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_12+(j*fnumx+ib), &tend3dyv_12[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_12[j][i] = (state_oldh_12[swindex[j]][i] + (((((((tend1dxh_12[swindex[j]][i] + tend1dyh_12[swindex[j]][i]) + tend2dxh_12[swindex[j]][i]) + tend2dyh_12[swindex[j]][i]) + (4.0 * tend3dxh_12[swindex[j]][i])) + (4.0 * tend3dyh_12[swindex[j]][i])) * dt) / 6.0));
            state_newu_12[j][i] = (state_oldu_12[swindex[j]][i] + (((((((tend1dxu_12[swindex[j]][i] + tend1dyu_12[swindex[j]][i]) + tend2dxu_12[swindex[j]][i]) + tend2dyu_12[swindex[j]][i]) + (4.0 * tend3dxu_12[swindex[j]][i])) + (4.0 * tend3dyu_12[swindex[j]][i])) * dt) / 6.0));
            state_newv_12[j][i] = (state_oldv_12[swindex[j]][i] + (((((((tend1dxv_12[swindex[j]][i] + tend1dyv_12[swindex[j]][i]) + tend2dxv_12[swindex[j]][i]) + tend2dyv_12[swindex[j]][i]) + (4.0 * tend3dxv_12[swindex[j]][i])) + (4.0 * tend3dyv_12[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_12+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_12[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_2(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_2_info *data = (update_rk3_0_2_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_13 = data->state_oldh_13;
  __attribute__ ((aligned (32))) double state_oldh_13[4][130];
  
  double *_tend1dxh_13 = data->tend1dxh_13;
  __attribute__ ((aligned (32))) double tend1dxh_13[4][130];
  
  double *_tend1dyh_13 = data->tend1dyh_13;
  __attribute__ ((aligned (32))) double tend1dyh_13[4][130];
  
  double *_state_oldu_13 = data->state_oldu_13;
  __attribute__ ((aligned (32))) double state_oldu_13[4][130];
  
  double *_tend1dxu_13 = data->tend1dxu_13;
  __attribute__ ((aligned (32))) double tend1dxu_13[4][130];
  
  double *_tend1dyu_13 = data->tend1dyu_13;
  __attribute__ ((aligned (32))) double tend1dyu_13[4][130];
  
  double *_state_oldv_13 = data->state_oldv_13;
  __attribute__ ((aligned (32))) double state_oldv_13[4][130];
  
  double *_tend1dxv_13 = data->tend1dxv_13;
  __attribute__ ((aligned (32))) double tend1dxv_13[4][130];
  
  double *_tend1dyv_13 = data->tend1dyv_13;
  __attribute__ ((aligned (32))) double tend1dyv_13[4][130];
  
  double *_tend2dxh_13 = data->tend2dxh_13;
  __attribute__ ((aligned (32))) double tend2dxh_13[4][130];
  
  double *_tend2dyh_13 = data->tend2dyh_13;
  __attribute__ ((aligned (32))) double tend2dyh_13[4][130];
  
  double *_tend3dxh_13 = data->tend3dxh_13;
  __attribute__ ((aligned (32))) double tend3dxh_13[4][130];
  
  double *_tend3dyh_13 = data->tend3dyh_13;
  __attribute__ ((aligned (32))) double tend3dyh_13[4][130];
  
  double *_tend2dxu_13 = data->tend2dxu_13;
  __attribute__ ((aligned (32))) double tend2dxu_13[4][130];
  
  double *_tend2dyu_13 = data->tend2dyu_13;
  __attribute__ ((aligned (32))) double tend2dyu_13[4][130];
  
  double *_tend3dxu_13 = data->tend3dxu_13;
  __attribute__ ((aligned (32))) double tend3dxu_13[4][130];
  
  double *_tend3dyu_13 = data->tend3dyu_13;
  __attribute__ ((aligned (32))) double tend3dyu_13[4][130];
  
  double *_tend2dxv_13 = data->tend2dxv_13;
  __attribute__ ((aligned (32))) double tend2dxv_13[4][130];
  
  double *_tend2dyv_13 = data->tend2dyv_13;
  __attribute__ ((aligned (32))) double tend2dyv_13[4][130];
  
  double *_tend3dxv_13 = data->tend3dxv_13;
  __attribute__ ((aligned (32))) double tend3dxv_13[4][130];
  
  double *_tend3dyv_13 = data->tend3dyv_13;
  __attribute__ ((aligned (32))) double tend3dyv_13[4][130];
  
  double *_state_newh_13 = data->state_newh_13;
  __attribute__ ((aligned (32))) double state_newh_13[4][130];
  
  double *_state_newu_13 = data->state_newu_13;
  __attribute__ ((aligned (32))) double state_newu_13[4][130];
  
  double *_state_newv_13 = data->state_newv_13;
  __attribute__ ((aligned (32))) double state_newv_13[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_13+(j*fnumx+ib), &tend2dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_13+(j*fnumx+ib), &tend2dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_13+(j*fnumx+ib), &tend3dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_13+(j*fnumx+ib), &tend3dxh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_13+(j*fnumx+ib), &tend3dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_13+(j*fnumx+ib), &tend3dyh_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_13+(j*fnumx+ib), &tend2dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_13+(j*fnumx+ib), &tend2dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_13+(j*fnumx+ib), &tend3dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_13+(j*fnumx+ib), &tend3dxu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_13+(j*fnumx+ib), &tend3dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_13+(j*fnumx+ib), &tend3dyu_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_13+(j*fnumx+ib), &tend2dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_13+(j*fnumx+ib), &tend2dyv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_13+(j*fnumx+ib), &tend3dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_13+(j*fnumx+ib), &tend3dxv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_13+(j*fnumx+ib), &tend3dyv_13[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_13+(j*fnumx+ib), &tend3dyv_13[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_13+(j*fnumx+ib), &state_oldh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_13+(j*fnumx+ib), &tend1dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_13+(j*fnumx+ib), &tend1dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_13+(j*fnumx+ib), &state_oldu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_13+(j*fnumx+ib), &tend1dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_13+(j*fnumx+ib), &tend1dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_13+(j*fnumx+ib), &state_oldv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_13+(j*fnumx+ib), &tend1dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_13+(j*fnumx+ib), &tend1dyv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_13+(j*fnumx+ib), &tend2dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_13+(j*fnumx+ib), &tend2dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_13+(j*fnumx+ib), &tend3dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_13+(j*fnumx+ib), &tend3dxh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_13+(j*fnumx+ib), &tend3dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_13+(j*fnumx+ib), &tend3dyh_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_13+(j*fnumx+ib), &tend2dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_13+(j*fnumx+ib), &tend2dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_13+(j*fnumx+ib), &tend3dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_13+(j*fnumx+ib), &tend3dxu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_13+(j*fnumx+ib), &tend3dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_13+(j*fnumx+ib), &tend3dyu_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_13+(j*fnumx+ib), &tend2dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_13+(j*fnumx+ib), &tend2dyv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_13+(j*fnumx+ib), &tend3dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_13+(j*fnumx+ib), &tend3dxv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_13+(j*fnumx+ib), &tend3dyv_13[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_13+(j*fnumx+ib), &tend3dyv_13[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_13[j][i] = (state_oldh_13[swindex[j]][i] + (((((((tend1dxh_13[swindex[j]][i] + tend1dyh_13[swindex[j]][i]) + tend2dxh_13[swindex[j]][i]) + tend2dyh_13[swindex[j]][i]) + (4.0 * tend3dxh_13[swindex[j]][i])) + (4.0 * tend3dyh_13[swindex[j]][i])) * dt) / 6.0));
            state_newu_13[j][i] = (state_oldu_13[swindex[j]][i] + (((((((tend1dxu_13[swindex[j]][i] + tend1dyu_13[swindex[j]][i]) + tend2dxu_13[swindex[j]][i]) + tend2dyu_13[swindex[j]][i]) + (4.0 * tend3dxu_13[swindex[j]][i])) + (4.0 * tend3dyu_13[swindex[j]][i])) * dt) / 6.0));
            state_newv_13[j][i] = (state_oldv_13[swindex[j]][i] + (((((((tend1dxv_13[swindex[j]][i] + tend1dyv_13[swindex[j]][i]) + tend2dxv_13[swindex[j]][i]) + tend2dyv_13[swindex[j]][i]) + (4.0 * tend3dxv_13[swindex[j]][i])) + (4.0 * tend3dyv_13[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_13+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_13[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_3(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_3_info *data = (update_rk3_0_3_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_21 = data->state_oldh_21;
  __attribute__ ((aligned (32))) double state_oldh_21[4][130];
  
  double *_tend1dxh_21 = data->tend1dxh_21;
  __attribute__ ((aligned (32))) double tend1dxh_21[4][130];
  
  double *_tend1dyh_21 = data->tend1dyh_21;
  __attribute__ ((aligned (32))) double tend1dyh_21[4][130];
  
  double *_state_oldu_21 = data->state_oldu_21;
  __attribute__ ((aligned (32))) double state_oldu_21[4][130];
  
  double *_tend1dxu_21 = data->tend1dxu_21;
  __attribute__ ((aligned (32))) double tend1dxu_21[4][130];
  
  double *_tend1dyu_21 = data->tend1dyu_21;
  __attribute__ ((aligned (32))) double tend1dyu_21[4][130];
  
  double *_state_oldv_21 = data->state_oldv_21;
  __attribute__ ((aligned (32))) double state_oldv_21[4][130];
  
  double *_tend1dxv_21 = data->tend1dxv_21;
  __attribute__ ((aligned (32))) double tend1dxv_21[4][130];
  
  double *_tend1dyv_21 = data->tend1dyv_21;
  __attribute__ ((aligned (32))) double tend1dyv_21[4][130];
  
  double *_tend2dxh_21 = data->tend2dxh_21;
  __attribute__ ((aligned (32))) double tend2dxh_21[4][130];
  
  double *_tend2dyh_21 = data->tend2dyh_21;
  __attribute__ ((aligned (32))) double tend2dyh_21[4][130];
  
  double *_tend3dxh_21 = data->tend3dxh_21;
  __attribute__ ((aligned (32))) double tend3dxh_21[4][130];
  
  double *_tend3dyh_21 = data->tend3dyh_21;
  __attribute__ ((aligned (32))) double tend3dyh_21[4][130];
  
  double *_tend2dxu_21 = data->tend2dxu_21;
  __attribute__ ((aligned (32))) double tend2dxu_21[4][130];
  
  double *_tend2dyu_21 = data->tend2dyu_21;
  __attribute__ ((aligned (32))) double tend2dyu_21[4][130];
  
  double *_tend3dxu_21 = data->tend3dxu_21;
  __attribute__ ((aligned (32))) double tend3dxu_21[4][130];
  
  double *_tend3dyu_21 = data->tend3dyu_21;
  __attribute__ ((aligned (32))) double tend3dyu_21[4][130];
  
  double *_tend2dxv_21 = data->tend2dxv_21;
  __attribute__ ((aligned (32))) double tend2dxv_21[4][130];
  
  double *_tend2dyv_21 = data->tend2dyv_21;
  __attribute__ ((aligned (32))) double tend2dyv_21[4][130];
  
  double *_tend3dxv_21 = data->tend3dxv_21;
  __attribute__ ((aligned (32))) double tend3dxv_21[4][130];
  
  double *_tend3dyv_21 = data->tend3dyv_21;
  __attribute__ ((aligned (32))) double tend3dyv_21[4][130];
  
  double *_state_newh_21 = data->state_newh_21;
  __attribute__ ((aligned (32))) double state_newh_21[4][130];
  
  double *_state_newu_21 = data->state_newu_21;
  __attribute__ ((aligned (32))) double state_newu_21[4][130];
  
  double *_state_newv_21 = data->state_newv_21;
  __attribute__ ((aligned (32))) double state_newv_21[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_21+(j*fnumx+ib), &tend2dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_21+(j*fnumx+ib), &tend2dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_21+(j*fnumx+ib), &tend3dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_21+(j*fnumx+ib), &tend3dxh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_21+(j*fnumx+ib), &tend3dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_21+(j*fnumx+ib), &tend3dyh_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_21+(j*fnumx+ib), &tend2dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_21+(j*fnumx+ib), &tend2dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_21+(j*fnumx+ib), &tend3dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_21+(j*fnumx+ib), &tend3dxu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_21+(j*fnumx+ib), &tend3dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_21+(j*fnumx+ib), &tend3dyu_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_21+(j*fnumx+ib), &tend2dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_21+(j*fnumx+ib), &tend2dyv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_21+(j*fnumx+ib), &tend3dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_21+(j*fnumx+ib), &tend3dxv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_21+(j*fnumx+ib), &tend3dyv_21[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_21+(j*fnumx+ib), &tend3dyv_21[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_21+(j*fnumx+ib), &state_oldh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_21+(j*fnumx+ib), &tend1dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_21+(j*fnumx+ib), &tend1dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_21+(j*fnumx+ib), &state_oldu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_21+(j*fnumx+ib), &tend1dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_21+(j*fnumx+ib), &tend1dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_21+(j*fnumx+ib), &state_oldv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_21+(j*fnumx+ib), &tend1dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_21+(j*fnumx+ib), &tend1dyv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_21+(j*fnumx+ib), &tend2dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_21+(j*fnumx+ib), &tend2dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_21+(j*fnumx+ib), &tend3dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_21+(j*fnumx+ib), &tend3dxh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_21+(j*fnumx+ib), &tend3dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_21+(j*fnumx+ib), &tend3dyh_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_21+(j*fnumx+ib), &tend2dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_21+(j*fnumx+ib), &tend2dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_21+(j*fnumx+ib), &tend3dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_21+(j*fnumx+ib), &tend3dxu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_21+(j*fnumx+ib), &tend3dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_21+(j*fnumx+ib), &tend3dyu_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_21+(j*fnumx+ib), &tend2dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_21+(j*fnumx+ib), &tend2dyv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_21+(j*fnumx+ib), &tend3dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_21+(j*fnumx+ib), &tend3dxv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_21+(j*fnumx+ib), &tend3dyv_21[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_21+(j*fnumx+ib), &tend3dyv_21[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_21[j][i] = (state_oldh_21[swindex[j]][i] + (((((((tend1dxh_21[swindex[j]][i] + tend1dyh_21[swindex[j]][i]) + tend2dxh_21[swindex[j]][i]) + tend2dyh_21[swindex[j]][i]) + (4.0 * tend3dxh_21[swindex[j]][i])) + (4.0 * tend3dyh_21[swindex[j]][i])) * dt) / 6.0));
            state_newu_21[j][i] = (state_oldu_21[swindex[j]][i] + (((((((tend1dxu_21[swindex[j]][i] + tend1dyu_21[swindex[j]][i]) + tend2dxu_21[swindex[j]][i]) + tend2dyu_21[swindex[j]][i]) + (4.0 * tend3dxu_21[swindex[j]][i])) + (4.0 * tend3dyu_21[swindex[j]][i])) * dt) / 6.0));
            state_newv_21[j][i] = (state_oldv_21[swindex[j]][i] + (((((((tend1dxv_21[swindex[j]][i] + tend1dyv_21[swindex[j]][i]) + tend2dxv_21[swindex[j]][i]) + tend2dyv_21[swindex[j]][i]) + (4.0 * tend3dxv_21[swindex[j]][i])) + (4.0 * tend3dyv_21[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_21+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_21[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_4(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_4_info *data = (update_rk3_0_4_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_22 = data->state_oldh_22;
  __attribute__ ((aligned (32))) double state_oldh_22[4][130];
  
  double *_tend1dxh_22 = data->tend1dxh_22;
  __attribute__ ((aligned (32))) double tend1dxh_22[4][130];
  
  double *_tend1dyh_22 = data->tend1dyh_22;
  __attribute__ ((aligned (32))) double tend1dyh_22[4][130];
  
  double *_state_oldu_22 = data->state_oldu_22;
  __attribute__ ((aligned (32))) double state_oldu_22[4][130];
  
  double *_tend1dxu_22 = data->tend1dxu_22;
  __attribute__ ((aligned (32))) double tend1dxu_22[4][130];
  
  double *_tend1dyu_22 = data->tend1dyu_22;
  __attribute__ ((aligned (32))) double tend1dyu_22[4][130];
  
  double *_state_oldv_22 = data->state_oldv_22;
  __attribute__ ((aligned (32))) double state_oldv_22[4][130];
  
  double *_tend1dxv_22 = data->tend1dxv_22;
  __attribute__ ((aligned (32))) double tend1dxv_22[4][130];
  
  double *_tend1dyv_22 = data->tend1dyv_22;
  __attribute__ ((aligned (32))) double tend1dyv_22[4][130];
  
  double *_tend2dxh_22 = data->tend2dxh_22;
  __attribute__ ((aligned (32))) double tend2dxh_22[4][130];
  
  double *_tend2dyh_22 = data->tend2dyh_22;
  __attribute__ ((aligned (32))) double tend2dyh_22[4][130];
  
  double *_tend3dxh_22 = data->tend3dxh_22;
  __attribute__ ((aligned (32))) double tend3dxh_22[4][130];
  
  double *_tend3dyh_22 = data->tend3dyh_22;
  __attribute__ ((aligned (32))) double tend3dyh_22[4][130];
  
  double *_tend2dxu_22 = data->tend2dxu_22;
  __attribute__ ((aligned (32))) double tend2dxu_22[4][130];
  
  double *_tend2dyu_22 = data->tend2dyu_22;
  __attribute__ ((aligned (32))) double tend2dyu_22[4][130];
  
  double *_tend3dxu_22 = data->tend3dxu_22;
  __attribute__ ((aligned (32))) double tend3dxu_22[4][130];
  
  double *_tend3dyu_22 = data->tend3dyu_22;
  __attribute__ ((aligned (32))) double tend3dyu_22[4][130];
  
  double *_tend2dxv_22 = data->tend2dxv_22;
  __attribute__ ((aligned (32))) double tend2dxv_22[4][130];
  
  double *_tend2dyv_22 = data->tend2dyv_22;
  __attribute__ ((aligned (32))) double tend2dyv_22[4][130];
  
  double *_tend3dxv_22 = data->tend3dxv_22;
  __attribute__ ((aligned (32))) double tend3dxv_22[4][130];
  
  double *_tend3dyv_22 = data->tend3dyv_22;
  __attribute__ ((aligned (32))) double tend3dyv_22[4][130];
  
  double *_state_newh_22 = data->state_newh_22;
  __attribute__ ((aligned (32))) double state_newh_22[4][130];
  
  double *_state_newu_22 = data->state_newu_22;
  __attribute__ ((aligned (32))) double state_newu_22[4][130];
  
  double *_state_newv_22 = data->state_newv_22;
  __attribute__ ((aligned (32))) double state_newv_22[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_22+(j*fnumx+ib), &tend2dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_22+(j*fnumx+ib), &tend2dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_22+(j*fnumx+ib), &tend3dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_22+(j*fnumx+ib), &tend3dxh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_22+(j*fnumx+ib), &tend3dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_22+(j*fnumx+ib), &tend3dyh_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_22+(j*fnumx+ib), &tend2dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_22+(j*fnumx+ib), &tend2dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_22+(j*fnumx+ib), &tend3dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_22+(j*fnumx+ib), &tend3dxu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_22+(j*fnumx+ib), &tend3dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_22+(j*fnumx+ib), &tend3dyu_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_22+(j*fnumx+ib), &tend2dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_22+(j*fnumx+ib), &tend2dyv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_22+(j*fnumx+ib), &tend3dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_22+(j*fnumx+ib), &tend3dxv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_22+(j*fnumx+ib), &tend3dyv_22[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_22+(j*fnumx+ib), &tend3dyv_22[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_22+(j*fnumx+ib), &state_oldh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_22+(j*fnumx+ib), &tend1dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_22+(j*fnumx+ib), &tend1dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_22+(j*fnumx+ib), &state_oldu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_22+(j*fnumx+ib), &tend1dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_22+(j*fnumx+ib), &tend1dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_22+(j*fnumx+ib), &state_oldv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_22+(j*fnumx+ib), &tend1dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_22+(j*fnumx+ib), &tend1dyv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_22+(j*fnumx+ib), &tend2dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_22+(j*fnumx+ib), &tend2dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_22+(j*fnumx+ib), &tend3dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_22+(j*fnumx+ib), &tend3dxh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_22+(j*fnumx+ib), &tend3dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_22+(j*fnumx+ib), &tend3dyh_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_22+(j*fnumx+ib), &tend2dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_22+(j*fnumx+ib), &tend2dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_22+(j*fnumx+ib), &tend3dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_22+(j*fnumx+ib), &tend3dxu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_22+(j*fnumx+ib), &tend3dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_22+(j*fnumx+ib), &tend3dyu_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_22+(j*fnumx+ib), &tend2dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_22+(j*fnumx+ib), &tend2dyv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_22+(j*fnumx+ib), &tend3dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_22+(j*fnumx+ib), &tend3dxv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_22+(j*fnumx+ib), &tend3dyv_22[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_22+(j*fnumx+ib), &tend3dyv_22[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_22[j][i] = (state_oldh_22[swindex[j]][i] + (((((((tend1dxh_22[swindex[j]][i] + tend1dyh_22[swindex[j]][i]) + tend2dxh_22[swindex[j]][i]) + tend2dyh_22[swindex[j]][i]) + (4.0 * tend3dxh_22[swindex[j]][i])) + (4.0 * tend3dyh_22[swindex[j]][i])) * dt) / 6.0));
            state_newu_22[j][i] = (state_oldu_22[swindex[j]][i] + (((((((tend1dxu_22[swindex[j]][i] + tend1dyu_22[swindex[j]][i]) + tend2dxu_22[swindex[j]][i]) + tend2dyu_22[swindex[j]][i]) + (4.0 * tend3dxu_22[swindex[j]][i])) + (4.0 * tend3dyu_22[swindex[j]][i])) * dt) / 6.0));
            state_newv_22[j][i] = (state_oldv_22[swindex[j]][i] + (((((((tend1dxv_22[swindex[j]][i] + tend1dyv_22[swindex[j]][i]) + tend2dxv_22[swindex[j]][i]) + tend2dyv_22[swindex[j]][i]) + (4.0 * tend3dxv_22[swindex[j]][i])) + (4.0 * tend3dyv_22[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_22+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_22[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_5(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_5_info *data = (update_rk3_0_5_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_23 = data->state_oldh_23;
  __attribute__ ((aligned (32))) double state_oldh_23[4][130];
  
  double *_tend1dxh_23 = data->tend1dxh_23;
  __attribute__ ((aligned (32))) double tend1dxh_23[4][130];
  
  double *_tend1dyh_23 = data->tend1dyh_23;
  __attribute__ ((aligned (32))) double tend1dyh_23[4][130];
  
  double *_state_oldu_23 = data->state_oldu_23;
  __attribute__ ((aligned (32))) double state_oldu_23[4][130];
  
  double *_tend1dxu_23 = data->tend1dxu_23;
  __attribute__ ((aligned (32))) double tend1dxu_23[4][130];
  
  double *_tend1dyu_23 = data->tend1dyu_23;
  __attribute__ ((aligned (32))) double tend1dyu_23[4][130];
  
  double *_state_oldv_23 = data->state_oldv_23;
  __attribute__ ((aligned (32))) double state_oldv_23[4][130];
  
  double *_tend1dxv_23 = data->tend1dxv_23;
  __attribute__ ((aligned (32))) double tend1dxv_23[4][130];
  
  double *_tend1dyv_23 = data->tend1dyv_23;
  __attribute__ ((aligned (32))) double tend1dyv_23[4][130];
  
  double *_tend2dxh_23 = data->tend2dxh_23;
  __attribute__ ((aligned (32))) double tend2dxh_23[4][130];
  
  double *_tend2dyh_23 = data->tend2dyh_23;
  __attribute__ ((aligned (32))) double tend2dyh_23[4][130];
  
  double *_tend3dxh_23 = data->tend3dxh_23;
  __attribute__ ((aligned (32))) double tend3dxh_23[4][130];
  
  double *_tend3dyh_23 = data->tend3dyh_23;
  __attribute__ ((aligned (32))) double tend3dyh_23[4][130];
  
  double *_tend2dxu_23 = data->tend2dxu_23;
  __attribute__ ((aligned (32))) double tend2dxu_23[4][130];
  
  double *_tend2dyu_23 = data->tend2dyu_23;
  __attribute__ ((aligned (32))) double tend2dyu_23[4][130];
  
  double *_tend3dxu_23 = data->tend3dxu_23;
  __attribute__ ((aligned (32))) double tend3dxu_23[4][130];
  
  double *_tend3dyu_23 = data->tend3dyu_23;
  __attribute__ ((aligned (32))) double tend3dyu_23[4][130];
  
  double *_tend2dxv_23 = data->tend2dxv_23;
  __attribute__ ((aligned (32))) double tend2dxv_23[4][130];
  
  double *_tend2dyv_23 = data->tend2dyv_23;
  __attribute__ ((aligned (32))) double tend2dyv_23[4][130];
  
  double *_tend3dxv_23 = data->tend3dxv_23;
  __attribute__ ((aligned (32))) double tend3dxv_23[4][130];
  
  double *_tend3dyv_23 = data->tend3dyv_23;
  __attribute__ ((aligned (32))) double tend3dyv_23[4][130];
  
  double *_state_newh_23 = data->state_newh_23;
  __attribute__ ((aligned (32))) double state_newh_23[4][130];
  
  double *_state_newu_23 = data->state_newu_23;
  __attribute__ ((aligned (32))) double state_newu_23[4][130];
  
  double *_state_newv_23 = data->state_newv_23;
  __attribute__ ((aligned (32))) double state_newv_23[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_23+(j*fnumx+ib), &tend2dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_23+(j*fnumx+ib), &tend2dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_23+(j*fnumx+ib), &tend3dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_23+(j*fnumx+ib), &tend3dxh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_23+(j*fnumx+ib), &tend3dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_23+(j*fnumx+ib), &tend3dyh_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_23+(j*fnumx+ib), &tend2dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_23+(j*fnumx+ib), &tend2dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_23+(j*fnumx+ib), &tend3dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_23+(j*fnumx+ib), &tend3dxu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_23+(j*fnumx+ib), &tend3dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_23+(j*fnumx+ib), &tend3dyu_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_23+(j*fnumx+ib), &tend2dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_23+(j*fnumx+ib), &tend2dyv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_23+(j*fnumx+ib), &tend3dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_23+(j*fnumx+ib), &tend3dxv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_23+(j*fnumx+ib), &tend3dyv_23[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_23+(j*fnumx+ib), &tend3dyv_23[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_23+(j*fnumx+ib), &state_oldh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_23+(j*fnumx+ib), &tend1dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_23+(j*fnumx+ib), &tend1dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_23+(j*fnumx+ib), &state_oldu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_23+(j*fnumx+ib), &tend1dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_23+(j*fnumx+ib), &tend1dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_23+(j*fnumx+ib), &state_oldv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_23+(j*fnumx+ib), &tend1dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_23+(j*fnumx+ib), &tend1dyv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_23+(j*fnumx+ib), &tend2dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_23+(j*fnumx+ib), &tend2dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_23+(j*fnumx+ib), &tend3dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_23+(j*fnumx+ib), &tend3dxh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_23+(j*fnumx+ib), &tend3dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_23+(j*fnumx+ib), &tend3dyh_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_23+(j*fnumx+ib), &tend2dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_23+(j*fnumx+ib), &tend2dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_23+(j*fnumx+ib), &tend3dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_23+(j*fnumx+ib), &tend3dxu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_23+(j*fnumx+ib), &tend3dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_23+(j*fnumx+ib), &tend3dyu_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_23+(j*fnumx+ib), &tend2dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_23+(j*fnumx+ib), &tend2dyv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_23+(j*fnumx+ib), &tend3dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_23+(j*fnumx+ib), &tend3dxv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_23+(j*fnumx+ib), &tend3dyv_23[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_23+(j*fnumx+ib), &tend3dyv_23[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_23[j][i] = (state_oldh_23[swindex[j]][i] + (((((((tend1dxh_23[swindex[j]][i] + tend1dyh_23[swindex[j]][i]) + tend2dxh_23[swindex[j]][i]) + tend2dyh_23[swindex[j]][i]) + (4.0 * tend3dxh_23[swindex[j]][i])) + (4.0 * tend3dyh_23[swindex[j]][i])) * dt) / 6.0));
            state_newu_23[j][i] = (state_oldu_23[swindex[j]][i] + (((((((tend1dxu_23[swindex[j]][i] + tend1dyu_23[swindex[j]][i]) + tend2dxu_23[swindex[j]][i]) + tend2dyu_23[swindex[j]][i]) + (4.0 * tend3dxu_23[swindex[j]][i])) + (4.0 * tend3dyu_23[swindex[j]][i])) * dt) / 6.0));
            state_newv_23[j][i] = (state_oldv_23[swindex[j]][i] + (((((((tend1dxv_23[swindex[j]][i] + tend1dyv_23[swindex[j]][i]) + tend2dxv_23[swindex[j]][i]) + tend2dyv_23[swindex[j]][i]) + (4.0 * tend3dxv_23[swindex[j]][i])) + (4.0 * tend3dyv_23[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_23+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_23[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_6(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_6_info *data = (update_rk3_0_6_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_31 = data->state_oldh_31;
  __attribute__ ((aligned (32))) double state_oldh_31[4][130];
  
  double *_tend1dxh_31 = data->tend1dxh_31;
  __attribute__ ((aligned (32))) double tend1dxh_31[4][130];
  
  double *_tend1dyh_31 = data->tend1dyh_31;
  __attribute__ ((aligned (32))) double tend1dyh_31[4][130];
  
  double *_state_oldu_31 = data->state_oldu_31;
  __attribute__ ((aligned (32))) double state_oldu_31[4][130];
  
  double *_tend1dxu_31 = data->tend1dxu_31;
  __attribute__ ((aligned (32))) double tend1dxu_31[4][130];
  
  double *_tend1dyu_31 = data->tend1dyu_31;
  __attribute__ ((aligned (32))) double tend1dyu_31[4][130];
  
  double *_state_oldv_31 = data->state_oldv_31;
  __attribute__ ((aligned (32))) double state_oldv_31[4][130];
  
  double *_tend1dxv_31 = data->tend1dxv_31;
  __attribute__ ((aligned (32))) double tend1dxv_31[4][130];
  
  double *_tend1dyv_31 = data->tend1dyv_31;
  __attribute__ ((aligned (32))) double tend1dyv_31[4][130];
  
  double *_tend2dxh_31 = data->tend2dxh_31;
  __attribute__ ((aligned (32))) double tend2dxh_31[4][130];
  
  double *_tend2dyh_31 = data->tend2dyh_31;
  __attribute__ ((aligned (32))) double tend2dyh_31[4][130];
  
  double *_tend3dxh_31 = data->tend3dxh_31;
  __attribute__ ((aligned (32))) double tend3dxh_31[4][130];
  
  double *_tend3dyh_31 = data->tend3dyh_31;
  __attribute__ ((aligned (32))) double tend3dyh_31[4][130];
  
  double *_tend2dxu_31 = data->tend2dxu_31;
  __attribute__ ((aligned (32))) double tend2dxu_31[4][130];
  
  double *_tend2dyu_31 = data->tend2dyu_31;
  __attribute__ ((aligned (32))) double tend2dyu_31[4][130];
  
  double *_tend3dxu_31 = data->tend3dxu_31;
  __attribute__ ((aligned (32))) double tend3dxu_31[4][130];
  
  double *_tend3dyu_31 = data->tend3dyu_31;
  __attribute__ ((aligned (32))) double tend3dyu_31[4][130];
  
  double *_tend2dxv_31 = data->tend2dxv_31;
  __attribute__ ((aligned (32))) double tend2dxv_31[4][130];
  
  double *_tend2dyv_31 = data->tend2dyv_31;
  __attribute__ ((aligned (32))) double tend2dyv_31[4][130];
  
  double *_tend3dxv_31 = data->tend3dxv_31;
  __attribute__ ((aligned (32))) double tend3dxv_31[4][130];
  
  double *_tend3dyv_31 = data->tend3dyv_31;
  __attribute__ ((aligned (32))) double tend3dyv_31[4][130];
  
  double *_state_newh_31 = data->state_newh_31;
  __attribute__ ((aligned (32))) double state_newh_31[4][130];
  
  double *_state_newu_31 = data->state_newu_31;
  __attribute__ ((aligned (32))) double state_newu_31[4][130];
  
  double *_state_newv_31 = data->state_newv_31;
  __attribute__ ((aligned (32))) double state_newv_31[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_31+(j*fnumx+ib), &tend2dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_31+(j*fnumx+ib), &tend2dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_31+(j*fnumx+ib), &tend3dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_31+(j*fnumx+ib), &tend3dxh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_31+(j*fnumx+ib), &tend3dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_31+(j*fnumx+ib), &tend3dyh_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_31+(j*fnumx+ib), &tend2dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_31+(j*fnumx+ib), &tend2dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_31+(j*fnumx+ib), &tend3dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_31+(j*fnumx+ib), &tend3dxu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_31+(j*fnumx+ib), &tend3dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_31+(j*fnumx+ib), &tend3dyu_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_31+(j*fnumx+ib), &tend2dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_31+(j*fnumx+ib), &tend2dyv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_31+(j*fnumx+ib), &tend3dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_31+(j*fnumx+ib), &tend3dxv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_31+(j*fnumx+ib), &tend3dyv_31[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_31+(j*fnumx+ib), &tend3dyv_31[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_31+(j*fnumx+ib), &state_oldh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_31+(j*fnumx+ib), &tend1dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_31+(j*fnumx+ib), &tend1dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_31+(j*fnumx+ib), &state_oldu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_31+(j*fnumx+ib), &tend1dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_31+(j*fnumx+ib), &tend1dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_31+(j*fnumx+ib), &state_oldv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_31+(j*fnumx+ib), &tend1dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_31+(j*fnumx+ib), &tend1dyv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_31+(j*fnumx+ib), &tend2dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_31+(j*fnumx+ib), &tend2dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_31+(j*fnumx+ib), &tend3dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_31+(j*fnumx+ib), &tend3dxh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_31+(j*fnumx+ib), &tend3dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_31+(j*fnumx+ib), &tend3dyh_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_31+(j*fnumx+ib), &tend2dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_31+(j*fnumx+ib), &tend2dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_31+(j*fnumx+ib), &tend3dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_31+(j*fnumx+ib), &tend3dxu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_31+(j*fnumx+ib), &tend3dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_31+(j*fnumx+ib), &tend3dyu_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_31+(j*fnumx+ib), &tend2dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_31+(j*fnumx+ib), &tend2dyv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_31+(j*fnumx+ib), &tend3dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_31+(j*fnumx+ib), &tend3dxv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_31+(j*fnumx+ib), &tend3dyv_31[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_31+(j*fnumx+ib), &tend3dyv_31[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_31[j][i] = (state_oldh_31[swindex[j]][i] + (((((((tend1dxh_31[swindex[j]][i] + tend1dyh_31[swindex[j]][i]) + tend2dxh_31[swindex[j]][i]) + tend2dyh_31[swindex[j]][i]) + (4.0 * tend3dxh_31[swindex[j]][i])) + (4.0 * tend3dyh_31[swindex[j]][i])) * dt) / 6.0));
            state_newu_31[j][i] = (state_oldu_31[swindex[j]][i] + (((((((tend1dxu_31[swindex[j]][i] + tend1dyu_31[swindex[j]][i]) + tend2dxu_31[swindex[j]][i]) + tend2dyu_31[swindex[j]][i]) + (4.0 * tend3dxu_31[swindex[j]][i])) + (4.0 * tend3dyu_31[swindex[j]][i])) * dt) / 6.0));
            state_newv_31[j][i] = (state_oldv_31[swindex[j]][i] + (((((((tend1dxv_31[swindex[j]][i] + tend1dyv_31[swindex[j]][i]) + tend2dxv_31[swindex[j]][i]) + tend2dyv_31[swindex[j]][i]) + (4.0 * tend3dxv_31[swindex[j]][i])) + (4.0 * tend3dyv_31[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_31+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_31[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_7(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_7_info *data = (update_rk3_0_7_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_32 = data->state_oldh_32;
  __attribute__ ((aligned (32))) double state_oldh_32[4][130];
  
  double *_tend1dxh_32 = data->tend1dxh_32;
  __attribute__ ((aligned (32))) double tend1dxh_32[4][130];
  
  double *_tend1dyh_32 = data->tend1dyh_32;
  __attribute__ ((aligned (32))) double tend1dyh_32[4][130];
  
  double *_state_oldu_32 = data->state_oldu_32;
  __attribute__ ((aligned (32))) double state_oldu_32[4][130];
  
  double *_tend1dxu_32 = data->tend1dxu_32;
  __attribute__ ((aligned (32))) double tend1dxu_32[4][130];
  
  double *_tend1dyu_32 = data->tend1dyu_32;
  __attribute__ ((aligned (32))) double tend1dyu_32[4][130];
  
  double *_state_oldv_32 = data->state_oldv_32;
  __attribute__ ((aligned (32))) double state_oldv_32[4][130];
  
  double *_tend1dxv_32 = data->tend1dxv_32;
  __attribute__ ((aligned (32))) double tend1dxv_32[4][130];
  
  double *_tend1dyv_32 = data->tend1dyv_32;
  __attribute__ ((aligned (32))) double tend1dyv_32[4][130];
  
  double *_tend2dxh_32 = data->tend2dxh_32;
  __attribute__ ((aligned (32))) double tend2dxh_32[4][130];
  
  double *_tend2dyh_32 = data->tend2dyh_32;
  __attribute__ ((aligned (32))) double tend2dyh_32[4][130];
  
  double *_tend3dxh_32 = data->tend3dxh_32;
  __attribute__ ((aligned (32))) double tend3dxh_32[4][130];
  
  double *_tend3dyh_32 = data->tend3dyh_32;
  __attribute__ ((aligned (32))) double tend3dyh_32[4][130];
  
  double *_tend2dxu_32 = data->tend2dxu_32;
  __attribute__ ((aligned (32))) double tend2dxu_32[4][130];
  
  double *_tend2dyu_32 = data->tend2dyu_32;
  __attribute__ ((aligned (32))) double tend2dyu_32[4][130];
  
  double *_tend3dxu_32 = data->tend3dxu_32;
  __attribute__ ((aligned (32))) double tend3dxu_32[4][130];
  
  double *_tend3dyu_32 = data->tend3dyu_32;
  __attribute__ ((aligned (32))) double tend3dyu_32[4][130];
  
  double *_tend2dxv_32 = data->tend2dxv_32;
  __attribute__ ((aligned (32))) double tend2dxv_32[4][130];
  
  double *_tend2dyv_32 = data->tend2dyv_32;
  __attribute__ ((aligned (32))) double tend2dyv_32[4][130];
  
  double *_tend3dxv_32 = data->tend3dxv_32;
  __attribute__ ((aligned (32))) double tend3dxv_32[4][130];
  
  double *_tend3dyv_32 = data->tend3dyv_32;
  __attribute__ ((aligned (32))) double tend3dyv_32[4][130];
  
  double *_state_newh_32 = data->state_newh_32;
  __attribute__ ((aligned (32))) double state_newh_32[4][130];
  
  double *_state_newu_32 = data->state_newu_32;
  __attribute__ ((aligned (32))) double state_newu_32[4][130];
  
  double *_state_newv_32 = data->state_newv_32;
  __attribute__ ((aligned (32))) double state_newv_32[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_32+(j*fnumx+ib), &tend2dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_32+(j*fnumx+ib), &tend2dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_32+(j*fnumx+ib), &tend3dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_32+(j*fnumx+ib), &tend3dxh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_32+(j*fnumx+ib), &tend3dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_32+(j*fnumx+ib), &tend3dyh_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_32+(j*fnumx+ib), &tend2dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_32+(j*fnumx+ib), &tend2dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_32+(j*fnumx+ib), &tend3dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_32+(j*fnumx+ib), &tend3dxu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_32+(j*fnumx+ib), &tend3dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_32+(j*fnumx+ib), &tend3dyu_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_32+(j*fnumx+ib), &tend2dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_32+(j*fnumx+ib), &tend2dyv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_32+(j*fnumx+ib), &tend3dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_32+(j*fnumx+ib), &tend3dxv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_32+(j*fnumx+ib), &tend3dyv_32[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_32+(j*fnumx+ib), &tend3dyv_32[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_32+(j*fnumx+ib), &state_oldh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_32+(j*fnumx+ib), &tend1dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_32+(j*fnumx+ib), &tend1dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_32+(j*fnumx+ib), &state_oldu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_32+(j*fnumx+ib), &tend1dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_32+(j*fnumx+ib), &tend1dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_32+(j*fnumx+ib), &state_oldv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_32+(j*fnumx+ib), &tend1dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_32+(j*fnumx+ib), &tend1dyv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_32+(j*fnumx+ib), &tend2dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_32+(j*fnumx+ib), &tend2dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_32+(j*fnumx+ib), &tend3dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_32+(j*fnumx+ib), &tend3dxh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_32+(j*fnumx+ib), &tend3dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_32+(j*fnumx+ib), &tend3dyh_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_32+(j*fnumx+ib), &tend2dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_32+(j*fnumx+ib), &tend2dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_32+(j*fnumx+ib), &tend3dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_32+(j*fnumx+ib), &tend3dxu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_32+(j*fnumx+ib), &tend3dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_32+(j*fnumx+ib), &tend3dyu_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_32+(j*fnumx+ib), &tend2dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_32+(j*fnumx+ib), &tend2dyv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_32+(j*fnumx+ib), &tend3dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_32+(j*fnumx+ib), &tend3dxv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_32+(j*fnumx+ib), &tend3dyv_32[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_32+(j*fnumx+ib), &tend3dyv_32[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_32[j][i] = (state_oldh_32[swindex[j]][i] + (((((((tend1dxh_32[swindex[j]][i] + tend1dyh_32[swindex[j]][i]) + tend2dxh_32[swindex[j]][i]) + tend2dyh_32[swindex[j]][i]) + (4.0 * tend3dxh_32[swindex[j]][i])) + (4.0 * tend3dyh_32[swindex[j]][i])) * dt) / 6.0));
            state_newu_32[j][i] = (state_oldu_32[swindex[j]][i] + (((((((tend1dxu_32[swindex[j]][i] + tend1dyu_32[swindex[j]][i]) + tend2dxu_32[swindex[j]][i]) + tend2dyu_32[swindex[j]][i]) + (4.0 * tend3dxu_32[swindex[j]][i])) + (4.0 * tend3dyu_32[swindex[j]][i])) * dt) / 6.0));
            state_newv_32[j][i] = (state_oldv_32[swindex[j]][i] + (((((((tend1dxv_32[swindex[j]][i] + tend1dyv_32[swindex[j]][i]) + tend2dxv_32[swindex[j]][i]) + tend2dyv_32[swindex[j]][i]) + (4.0 * tend3dxv_32[swindex[j]][i])) + (4.0 * tend3dyv_32[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_32+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_32[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void update_rk3_0_8(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  update_rk3_0_8_info *data = (update_rk3_0_8_info*)(_ptr);
  
  int lx = data->lx;
  int ly = data->ly;
  int lz = data->lz;
  int hx = data->hx;
  int hy = data->hy;
  int hz = data->hz;
  int ox = data->ox;
  int oy = data->oy;
  int oz = data->oz;
  int bx = data->bx;
  int by = data->by;
  int bz = data->bz;
  int mx = data->mx;
  int my = data->my;
  int mz = data->mz;
  
  
  double *_state_oldh_33 = data->state_oldh_33;
  __attribute__ ((aligned (32))) double state_oldh_33[4][130];
  
  double *_tend1dxh_33 = data->tend1dxh_33;
  __attribute__ ((aligned (32))) double tend1dxh_33[4][130];
  
  double *_tend1dyh_33 = data->tend1dyh_33;
  __attribute__ ((aligned (32))) double tend1dyh_33[4][130];
  
  double *_state_oldu_33 = data->state_oldu_33;
  __attribute__ ((aligned (32))) double state_oldu_33[4][130];
  
  double *_tend1dxu_33 = data->tend1dxu_33;
  __attribute__ ((aligned (32))) double tend1dxu_33[4][130];
  
  double *_tend1dyu_33 = data->tend1dyu_33;
  __attribute__ ((aligned (32))) double tend1dyu_33[4][130];
  
  double *_state_oldv_33 = data->state_oldv_33;
  __attribute__ ((aligned (32))) double state_oldv_33[4][130];
  
  double *_tend1dxv_33 = data->tend1dxv_33;
  __attribute__ ((aligned (32))) double tend1dxv_33[4][130];
  
  double *_tend1dyv_33 = data->tend1dyv_33;
  __attribute__ ((aligned (32))) double tend1dyv_33[4][130];
  
  double *_tend2dxh_33 = data->tend2dxh_33;
  __attribute__ ((aligned (32))) double tend2dxh_33[4][130];
  
  double *_tend2dyh_33 = data->tend2dyh_33;
  __attribute__ ((aligned (32))) double tend2dyh_33[4][130];
  
  double *_tend3dxh_33 = data->tend3dxh_33;
  __attribute__ ((aligned (32))) double tend3dxh_33[4][130];
  
  double *_tend3dyh_33 = data->tend3dyh_33;
  __attribute__ ((aligned (32))) double tend3dyh_33[4][130];
  
  double *_tend2dxu_33 = data->tend2dxu_33;
  __attribute__ ((aligned (32))) double tend2dxu_33[4][130];
  
  double *_tend2dyu_33 = data->tend2dyu_33;
  __attribute__ ((aligned (32))) double tend2dyu_33[4][130];
  
  double *_tend3dxu_33 = data->tend3dxu_33;
  __attribute__ ((aligned (32))) double tend3dxu_33[4][130];
  
  double *_tend3dyu_33 = data->tend3dyu_33;
  __attribute__ ((aligned (32))) double tend3dyu_33[4][130];
  
  double *_tend2dxv_33 = data->tend2dxv_33;
  __attribute__ ((aligned (32))) double tend2dxv_33[4][130];
  
  double *_tend2dyv_33 = data->tend2dyv_33;
  __attribute__ ((aligned (32))) double tend2dyv_33[4][130];
  
  double *_tend3dxv_33 = data->tend3dxv_33;
  __attribute__ ((aligned (32))) double tend3dxv_33[4][130];
  
  double *_tend3dyv_33 = data->tend3dyv_33;
  __attribute__ ((aligned (32))) double tend3dyv_33[4][130];
  
  double *_state_newh_33 = data->state_newh_33;
  __attribute__ ((aligned (32))) double state_newh_33[4][130];
  
  double *_state_newu_33 = data->state_newu_33;
  __attribute__ ((aligned (32))) double state_newu_33[4][130];
  
  double *_state_newv_33 = data->state_newv_33;
  __attribute__ ((aligned (32))) double state_newv_33[4][130];
  
  
  double dt = data->dt;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[by+2*hy];
  const int ws  = by+2*hy;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_33+(j*fnumx+ib), &tend2dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_33+(j*fnumx+ib), &tend2dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_33+(j*fnumx+ib), &tend3dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_33+(j*fnumx+ib), &tend3dxh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_33+(j*fnumx+ib), &tend3dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_33+(j*fnumx+ib), &tend3dyh_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_33+(j*fnumx+ib), &tend2dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_33+(j*fnumx+ib), &tend2dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_33+(j*fnumx+ib), &tend3dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_33+(j*fnumx+ib), &tend3dxu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_33+(j*fnumx+ib), &tend3dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_33+(j*fnumx+ib), &tend3dyu_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_33+(j*fnumx+ib), &tend2dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_33+(j*fnumx+ib), &tend2dyv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_33+(j*fnumx+ib), &tend3dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_33+(j*fnumx+ib), &tend3dxv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_33+(j*fnumx+ib), &tend3dyv_33[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_33+(j*fnumx+ib), &tend3dyv_33[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_oldh_33+(j*fnumx+ib), &state_oldh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxh_33+(j*fnumx+ib), &tend1dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyh_33+(j*fnumx+ib), &tend1dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldu_33+(j*fnumx+ib), &state_oldu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxu_33+(j*fnumx+ib), &tend1dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyu_33+(j*fnumx+ib), &tend1dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_oldv_33+(j*fnumx+ib), &state_oldv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dxv_33+(j*fnumx+ib), &tend1dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend1dyv_33+(j*fnumx+ib), &tend1dyv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxh_33+(j*fnumx+ib), &tend2dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyh_33+(j*fnumx+ib), &tend2dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_33+(j*fnumx+ib), &tend3dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxh_33+(j*fnumx+ib), &tend3dxh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_33+(j*fnumx+ib), &tend3dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyh_33+(j*fnumx+ib), &tend3dyh_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxu_33+(j*fnumx+ib), &tend2dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyu_33+(j*fnumx+ib), &tend2dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_33+(j*fnumx+ib), &tend3dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxu_33+(j*fnumx+ib), &tend3dxu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_33+(j*fnumx+ib), &tend3dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyu_33+(j*fnumx+ib), &tend3dyu_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dxv_33+(j*fnumx+ib), &tend2dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend2dyv_33+(j*fnumx+ib), &tend2dyv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_33+(j*fnumx+ib), &tend3dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dxv_33+(j*fnumx+ib), &tend3dxv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_33+(j*fnumx+ib), &tend3dyv_33[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tend3dyv_33+(j*fnumx+ib), &tend3dyv_33[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_newh_33[j][i] = (state_oldh_33[swindex[j]][i] + (((((((tend1dxh_33[swindex[j]][i] + tend1dyh_33[swindex[j]][i]) + tend2dxh_33[swindex[j]][i]) + tend2dyh_33[swindex[j]][i]) + (4.0 * tend3dxh_33[swindex[j]][i])) + (4.0 * tend3dyh_33[swindex[j]][i])) * dt) / 6.0));
            state_newu_33[j][i] = (state_oldu_33[swindex[j]][i] + (((((((tend1dxu_33[swindex[j]][i] + tend1dyu_33[swindex[j]][i]) + tend2dxu_33[swindex[j]][i]) + tend2dyu_33[swindex[j]][i]) + (4.0 * tend3dxu_33[swindex[j]][i])) + (4.0 * tend3dyu_33[swindex[j]][i])) * dt) / 6.0));
            state_newv_33[j][i] = (state_oldv_33[swindex[j]][i] + (((((((tend1dxv_33[swindex[j]][i] + tend1dyv_33[swindex[j]][i]) + tend2dxv_33[swindex[j]][i]) + tend2dyv_33[swindex[j]][i]) + (4.0 * tend3dxv_33[swindex[j]][i])) + (4.0 * tend3dyv_33[swindex[j]][i])) * dt) / 6.0));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_newh_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newh_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newu_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newu_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_newv_33+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_newv_33[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

