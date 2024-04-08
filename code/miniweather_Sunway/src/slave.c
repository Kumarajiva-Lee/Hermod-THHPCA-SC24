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

struct InitVector slave_thermal(double x, double z)
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









void StateInitial_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  StateInitial_0_0_info *data = (StateInitial_0_0_info*)(_ptr);
  
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
  
  
  double *_staterho = data->staterho;
  __attribute__ ((aligned (32))) double staterho[by+2*hy][bx+2*hx];
  
  double *_stateu = data->stateu;
  __attribute__ ((aligned (32))) double stateu[by+2*hy][bx+2*hx];
  
  double *_statew = data->statew;
  __attribute__ ((aligned (32))) double statew[by+2*hy][bx+2*hx];
  
  double *_statept = data->statept;
  __attribute__ ((aligned (32))) double statept[by+2*hy][bx+2*hx];
  
  double *_xp = data->xp;
  __attribute__ ((aligned (32))) double xp[fnx+2*hx];
  
  double *_zpf = data->zpf;
  __attribute__ ((aligned (32))) double zpf[fnz+2*hz];
  
  int kk;
  int ii;
  double qpi;
  double qwi;
  double qpk;
  double qwk;
  double x;
  double z;
  struct InitVector InitV;
  
  
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
  
  dmasize = fnx + 2*hx;
  DMA_READ(_xp+(ghx-hx), &xp[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_zpf+(ghz-hz), &zpf[0], dmasize* sizeof(double));
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            for (kk=0;kk<3;kk+=1){
              for (ii=0;ii<3;ii+=1){
                if ((ii == 0))
                {
                  qpi = 0.11270166537925831;
                  qwi = 0.2777777777777778;
                }
                else
                {
                  if ((ii == 1))
                  {
                    qpi = 0.5;
                    qwi = 0.4444444444444444;
                  }
                  else
                  {
                    if ((ii == 2))
                    {
                      qpi = 0.8872983346207417;
                      qwi = 0.2777777777777778;
                    }
                  }
                }
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
                x = (xp[i+ib+ox] + ((qpi - 0.5) * 50.0));
                z = (zpf[kb+oz] + ((qpk - 0.5) * 50.0));
                InitV = (struct InitVector){0.0,0.0,0.0,0.0,0.0,0.0};
                InitV = slave_thermal(x,z);
                staterho[j][i] = (staterho[j][i] + ((InitV.r * qwi) * qwk));
                stateu[j][i] = (stateu[j][i] + ((((InitV.r + InitV.hr) * InitV.u) * qwi) * qwk));
                statew[j][i] = (statew[j][i] + ((((InitV.r + InitV.hr) * InitV.w) * qwi) * qwk));
                statept[j][i] = (statept[j][i] + (((((InitV.r + InitV.hr) * (InitV.t + InitV.ht)) - (InitV.hr * InitV.ht)) * qwi) * qwk));
              }
            }
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_staterho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &staterho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_stateu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &stateu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_statew+(kb*fnumy*fnumx+j*fnumx+ib+hx), &statew[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_statept+(kb*fnumy*fnumx+j*fnumx+ib+hx), &statept[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void StateInitial_0_1(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  StateInitial_0_1_info *data = (StateInitial_0_1_info*)(_ptr);
  
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
  
  
  double *_staterho = data->staterho;
  __attribute__ ((aligned (32))) double staterho[8][260];
  
  double *_stateu = data->stateu;
  __attribute__ ((aligned (32))) double stateu[8][260];
  
  double *_statew = data->statew;
  __attribute__ ((aligned (32))) double statew[8][260];
  
  double *_statept = data->statept;
  __attribute__ ((aligned (32))) double statept[8][260];
  
  double *_state_tmprho = data->state_tmprho;
  __attribute__ ((aligned (32))) double state_tmprho[8][260];
  
  double *_state_tmpu = data->state_tmpu;
  __attribute__ ((aligned (32))) double state_tmpu[8][260];
  
  double *_state_tmpw = data->state_tmpw;
  __attribute__ ((aligned (32))) double state_tmpw[8][260];
  
  double *_state_tmppt = data->state_tmppt;
  __attribute__ ((aligned (32))) double state_tmppt[8][260];
  
  
  
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
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_tmprho[j][i] = staterho[swindex[j]][i];
            state_tmpu[j][i] = stateu[swindex[j]][i];
            state_tmpw[j][i] = statew[swindex[j]][i];
            state_tmppt[j][i] = statept[swindex[j]][i];
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_tmprho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_tmprho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_tmpu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_tmpu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_tmpw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_tmpw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_tmppt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_tmppt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

void ComputeTendenciesX_0_20(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesX_0_20_info *data = (ComputeTendenciesX_0_20_info*)(_ptr);
  
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
  
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_state_forcingu = data->state_forcingu;
  __attribute__ ((aligned (32))) double state_forcingu[8][260];
  
  double *_state_forcingw = data->state_forcingw;
  __attribute__ ((aligned (32))) double state_forcingw[8][260];
  
  double *_state_forcingpt = data->state_forcingpt;
  __attribute__ ((aligned (32))) double state_forcingpt[8][260];
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_hy_dens_cell = data->hy_dens_cell;
  __attribute__ ((aligned (32))) double hy_dens_cell[fnz+2*hz];
  
  double *_hy_dens_theta_cell = data->hy_dens_theta_cell;
  __attribute__ ((aligned (32))) double hy_dens_theta_cell[fnz+2*hz];
  
  struct Vector4 vals;
  struct Vector4 d3_vals;
  double r;
  double u;
  double w;
  double t;
  double p;
  
  double hv_coef = data->hv_coef;
  
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
  
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_cell+(ghz-hz), &hy_dens_cell[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_theta_cell+(ghz-hz), &hy_dens_theta_cell[0], dmasize* sizeof(double));
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_forcingrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingpt[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_forcingrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vals.r = ((((-state_forcingrho[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingrho[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingrho[swindex[j]][i]) / 12)) - (state_forcingrho[swindex[j]][(i + 1)] / 12));
            d3_vals.r = (((-state_forcingrho[swindex[j]][(i - 2)] + (3 * state_forcingrho[swindex[j]][(i - 1)])) - (3 * state_forcingrho[swindex[j]][i])) + state_forcingrho[swindex[j]][(i + 1)]);
            vals.u = ((((-state_forcingu[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingu[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingu[swindex[j]][i]) / 12)) - (state_forcingu[swindex[j]][(i + 1)] / 12));
            d3_vals.u = (((-state_forcingu[swindex[j]][(i - 2)] + (3 * state_forcingu[swindex[j]][(i - 1)])) - (3 * state_forcingu[swindex[j]][i])) + state_forcingu[swindex[j]][(i + 1)]);
            vals.w = ((((-state_forcingw[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingw[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingw[swindex[j]][i]) / 12)) - (state_forcingw[swindex[j]][(i + 1)] / 12));
            d3_vals.w = (((-state_forcingw[swindex[j]][(i - 2)] + (3 * state_forcingw[swindex[j]][(i - 1)])) - (3 * state_forcingw[swindex[j]][i])) + state_forcingw[swindex[j]][(i + 1)]);
            vals.t = ((((-state_forcingpt[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingpt[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingpt[swindex[j]][i]) / 12)) - (state_forcingpt[swindex[j]][(i + 1)] / 12));
            d3_vals.t = (((-state_forcingpt[swindex[j]][(i - 2)] + (3 * state_forcingpt[swindex[j]][(i - 1)])) - (3 * state_forcingpt[swindex[j]][i])) + state_forcingpt[swindex[j]][(i + 1)]);
            r = (vals.r + hy_dens_cell[kb+oz]);
            u = (vals.u / r);
            w = (vals.w / r);
            t = ((vals.t + hy_dens_theta_cell[kb+oz]) / r);
            p = (27.562941092972594 * pow((r * t),1.400278940027894));
            fluxfrho[j][i] = ((r * u) - (hv_coef * d3_vals.r));
            fluxfu[j][i] = ((((r * u) * u) + p) - (hv_coef * d3_vals.u));
            fluxfw[j][i] = (((r * u) * w) - (hv_coef * d3_vals.w));
            fluxfpt[j][i] = (((r * u) * t) - (hv_coef * d3_vals.t));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfrho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfpt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesX_0_21(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesX_0_21_info *data = (ComputeTendenciesX_0_21_info*)(_ptr);
  
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
  
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][258];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][258];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][258];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][258];
  
  double *_tenddrho = data->tenddrho;
  __attribute__ ((aligned (32))) double tenddrho[8][258];
  
  double *_tenddu = data->tenddu;
  __attribute__ ((aligned (32))) double tenddu[8][258];
  
  double *_tenddw = data->tenddw;
  __attribute__ ((aligned (32))) double tenddw[8][258];
  
  double *_tenddpt = data->tenddpt;
  __attribute__ ((aligned (32))) double tenddpt[8][258];
  
  
  
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
            DMA_IREAD(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfpt[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            tenddrho[j][i] = (-(fluxfrho[swindex[j]][(i + 1)] - fluxfrho[swindex[j]][i]) / 50.0);
            tenddu[j][i] = (-(fluxfu[swindex[j]][(i + 1)] - fluxfu[swindex[j]][i]) / 50.0);
            tenddw[j][i] = (-(fluxfw[swindex[j]][(i + 1)] - fluxfw[swindex[j]][i]) / 50.0);
            tenddpt[j][i] = (-(fluxfpt[swindex[j]][(i + 1)] - fluxfpt[swindex[j]][i]) / 50.0);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddrho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddrho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddpt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddpt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesZ_0_20(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesZ_0_20_info *data = (ComputeTendenciesZ_0_20_info*)(_ptr);
  
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
  
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_state_forcingu = data->state_forcingu;
  __attribute__ ((aligned (32))) double state_forcingu[8][260];
  
  double *_state_forcingw = data->state_forcingw;
  __attribute__ ((aligned (32))) double state_forcingw[8][260];
  
  double *_state_forcingpt = data->state_forcingpt;
  __attribute__ ((aligned (32))) double state_forcingpt[8][260];
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_hy_dens_int = data->hy_dens_int;
  __attribute__ ((aligned (32))) double hy_dens_int[hnz+2*hz];
  
  double *_hy_dens_theta_int = data->hy_dens_theta_int;
  __attribute__ ((aligned (32))) double hy_dens_theta_int[hnz+2*hz];
  
  double *_hy_pressure_int = data->hy_pressure_int;
  __attribute__ ((aligned (32))) double hy_pressure_int[hnz+2*hz];
  
  struct Vector4 vals;
  struct Vector4 d3_vals;
  double r;
  double u;
  double w;
  double t;
  double p;
  
  double hv_coef = data->hv_coef;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange,krange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[bz+2*hz];
  const int ws  = bz+2*hz;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j,k;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_int+(ghz-hz), &hy_dens_int[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_theta_int+(ghz-hz), &hy_dens_theta_int[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_pressure_int+(ghz-hz), &hy_pressure_int[0], dmasize* sizeof(double));
  
  for (jb = ybeg; jb < yend; jb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (kb= zbeg; kb < zend ; kb+= bz)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        krange = MIN(kb+bz,zend) - kb;
        dmasize = (irange + 2*hx);
        
        if (kb == zbeg){
          for (k = kb ; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingu[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingpt[k - kb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + bz) % ws;
          for (k = kb + hz; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            int pos = swindex[k - kb];
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (k = hz ; k < hz + krange ; k++)
          for (i = hx ; i < hx + irange ; i++){
            vals.r = ((((-state_forcingrho[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingrho[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingrho[swindex[k]][i]) / 12)) - (state_forcingrho[swindex[(k + 1)]][i] / 12));
            d3_vals.r = (((-state_forcingrho[swindex[(k - 2)]][i] + (3 * state_forcingrho[swindex[(k - 1)]][i])) - (3 * state_forcingrho[swindex[k]][i])) + state_forcingrho[swindex[(k + 1)]][i]);
            vals.u = ((((-state_forcingu[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingu[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingu[swindex[k]][i]) / 12)) - (state_forcingu[swindex[(k + 1)]][i] / 12));
            d3_vals.u = (((-state_forcingu[swindex[(k - 2)]][i] + (3 * state_forcingu[swindex[(k - 1)]][i])) - (3 * state_forcingu[swindex[k]][i])) + state_forcingu[swindex[(k + 1)]][i]);
            vals.w = ((((-state_forcingw[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingw[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingw[swindex[k]][i]) / 12)) - (state_forcingw[swindex[(k + 1)]][i] / 12));
            d3_vals.w = (((-state_forcingw[swindex[(k - 2)]][i] + (3 * state_forcingw[swindex[(k - 1)]][i])) - (3 * state_forcingw[swindex[k]][i])) + state_forcingw[swindex[(k + 1)]][i]);
            vals.t = ((((-state_forcingpt[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingpt[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingpt[swindex[k]][i]) / 12)) - (state_forcingpt[swindex[(k + 1)]][i] / 12));
            d3_vals.t = (((-state_forcingpt[swindex[(k - 2)]][i] + (3 * state_forcingpt[swindex[(k - 1)]][i])) - (3 * state_forcingpt[swindex[k]][i])) + state_forcingpt[swindex[(k + 1)]][i]);
            r = (vals.r + hy_dens_int[k+kb+oz]);
            u = (vals.u / r);
            w = (vals.w / r);
            t = ((vals.t + hy_dens_theta_int[k+kb+oz]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - hy_pressure_int[k+kb+oz]);
            fluxfrho[k][i] = ((r * w) - (hv_coef * d3_vals.r));
            fluxfu[k][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            fluxfw[k][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            fluxfpt[k][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        
        for (k = kb + hz  ; k < MIN(kb+ bz, zend) + hz ; k++){
          DMA_IWRITE(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfrho[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfu[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfw[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfpt[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesZ_0_21(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesZ_0_21_info *data = (ComputeTendenciesZ_0_21_info*)(_ptr);
  
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
  
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_tenddw = data->tenddw;
  __attribute__ ((aligned (32))) double tenddw[8][260];
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_tenddrho = data->tenddrho;
  __attribute__ ((aligned (32))) double tenddrho[8][260];
  
  double *_tenddu = data->tenddu;
  __attribute__ ((aligned (32))) double tenddu[8][260];
  
  double *_tenddpt = data->tenddpt;
  __attribute__ ((aligned (32))) double tenddpt[8][260];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange,krange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[bz+2*hz];
  const int ws  = bz+2*hz;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j,k;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (jb = ybeg; jb < yend; jb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (kb= zbeg; kb < zend ; kb+= bz)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        krange = MIN(kb+bz,zend) - kb;
        dmasize = (irange + 2*hx);
        
        if (kb == zbeg){
          for (k = kb ; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            DMA_IREAD(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfu[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib), &tenddw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfpt[k - kb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + bz) % ws;
          for (k = kb + hz; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            int pos = swindex[k - kb];
            DMA_IREAD(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib), &tenddw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (k = hz ; k < hz + krange ; k++)
          for (i = hx ; i < hx + irange ; i++){
            tenddrho[k][i] = (-(fluxfrho[swindex[(k + 1)]][i] - fluxfrho[swindex[k]][i]) / 50.0);
            tenddu[k][i] = (-(fluxfu[swindex[(k + 1)]][i] - fluxfu[swindex[k]][i]) / 50.0);
            tenddw[k][i] = (-(fluxfw[swindex[(k + 1)]][i] - fluxfw[swindex[k]][i]) / 50.0);
            tenddw[k][i] = (tenddw[k][i] - (state_forcingrho[swindex[k]][i] * 9.8));
            tenddpt[k][i] = (-(fluxfpt[swindex[(k + 1)]][i] - fluxfpt[swindex[k]][i]) / 50.0);
          }
        
        for (k = kb + hz  ; k < MIN(kb+ bz, zend) + hz ; k++){
          DMA_IWRITE(_tenddrho+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddrho[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddu+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddu[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddw[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddpt+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddpt[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesZ_0_22(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesZ_0_22_info *data = (ComputeTendenciesZ_0_22_info*)(_ptr);
  
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
  
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_state_forcingu = data->state_forcingu;
  __attribute__ ((aligned (32))) double state_forcingu[8][260];
  
  double *_state_forcingw = data->state_forcingw;
  __attribute__ ((aligned (32))) double state_forcingw[8][260];
  
  double *_state_forcingpt = data->state_forcingpt;
  __attribute__ ((aligned (32))) double state_forcingpt[8][260];
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_hy_dens_int = data->hy_dens_int;
  __attribute__ ((aligned (32))) double hy_dens_int[hnz+2*hz];
  
  double *_hy_dens_theta_int = data->hy_dens_theta_int;
  __attribute__ ((aligned (32))) double hy_dens_theta_int[hnz+2*hz];
  
  double *_hy_pressure_int = data->hy_pressure_int;
  __attribute__ ((aligned (32))) double hy_pressure_int[hnz+2*hz];
  
  struct Vector4 vals;
  struct Vector4 d3_vals;
  double r;
  double u;
  double w;
  double t;
  double p;
  
  double hv_coef = data->hv_coef;
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange,krange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[bz+2*hz];
  const int ws  = bz+2*hz;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j,k;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_int+(ghz-hz), &hy_dens_int[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_theta_int+(ghz-hz), &hy_dens_theta_int[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_pressure_int+(ghz-hz), &hy_pressure_int[0], dmasize* sizeof(double));
  
  for (jb = ybeg; jb < yend; jb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (kb= zbeg; kb < zend ; kb+= bz)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        krange = MIN(kb+bz,zend) - kb;
        dmasize = (irange + 2*hx);
        
        if (kb == zbeg){
          for (k = kb ; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingu[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingpt[k - kb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + bz) % ws;
          for (k = kb + hz; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            int pos = swindex[k - kb];
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (k = hz ; k < hz + krange ; k++)
          for (i = hx ; i < hx + irange ; i++){
            vals.r = ((((-state_forcingrho[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingrho[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingrho[swindex[k]][i]) / 12)) - (state_forcingrho[swindex[(k + 1)]][i] / 12));
            d3_vals.r = (((-state_forcingrho[swindex[(k - 2)]][i] + (3 * state_forcingrho[swindex[(k - 1)]][i])) - (3 * state_forcingrho[swindex[k]][i])) + state_forcingrho[swindex[(k + 1)]][i]);
            vals.u = ((((-state_forcingu[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingu[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingu[swindex[k]][i]) / 12)) - (state_forcingu[swindex[(k + 1)]][i] / 12));
            d3_vals.u = (((-state_forcingu[swindex[(k - 2)]][i] + (3 * state_forcingu[swindex[(k - 1)]][i])) - (3 * state_forcingu[swindex[k]][i])) + state_forcingu[swindex[(k + 1)]][i]);
            vals.w = ((((-state_forcingw[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingw[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingw[swindex[k]][i]) / 12)) - (state_forcingw[swindex[(k + 1)]][i] / 12));
            d3_vals.w = (((-state_forcingw[swindex[(k - 2)]][i] + (3 * state_forcingw[swindex[(k - 1)]][i])) - (3 * state_forcingw[swindex[k]][i])) + state_forcingw[swindex[(k + 1)]][i]);
            vals.t = ((((-state_forcingpt[swindex[(k - 2)]][i] / 12) + ((7 * state_forcingpt[swindex[(k - 1)]][i]) / 12)) + ((7 * state_forcingpt[swindex[k]][i]) / 12)) - (state_forcingpt[swindex[(k + 1)]][i] / 12));
            d3_vals.t = (((-state_forcingpt[swindex[(k - 2)]][i] + (3 * state_forcingpt[swindex[(k - 1)]][i])) - (3 * state_forcingpt[swindex[k]][i])) + state_forcingpt[swindex[(k + 1)]][i]);
            r = (vals.r + hy_dens_int[k+kb+oz]);
            u = (vals.u / r);
            w = (vals.w / r);
            t = ((vals.t + hy_dens_theta_int[k+kb+oz]) / r);
            p = ((27.562941092972594 * pow((r * t),1.400278940027894)) - hy_pressure_int[k+kb+oz]);
            fluxfrho[k][i] = ((r * w) - (hv_coef * d3_vals.r));
            fluxfu[k][i] = (((r * w) * u) - (hv_coef * d3_vals.u));
            fluxfw[k][i] = ((((r * w) * w) + p) - (hv_coef * d3_vals.w));
            fluxfpt[k][i] = (((r * w) * t) - (hv_coef * d3_vals.t));
          }
        
        for (k = kb + hz  ; k < MIN(kb+ bz, zend) + hz ; k++){
          DMA_IWRITE(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfrho[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfu[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfw[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib+hx), &fluxfpt[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesZ_0_23(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesZ_0_23_info *data = (ComputeTendenciesZ_0_23_info*)(_ptr);
  
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
  
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_tenddw = data->tenddw;
  __attribute__ ((aligned (32))) double tenddw[8][260];
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_tenddrho = data->tenddrho;
  __attribute__ ((aligned (32))) double tenddrho[8][260];
  
  double *_tenddu = data->tenddu;
  __attribute__ ((aligned (32))) double tenddu[8][260];
  
  double *_tenddpt = data->tenddpt;
  __attribute__ ((aligned (32))) double tenddpt[8][260];
  
  
  
  int zbeg,zend,zlen,ybeg,yend,ylen,xbeg,xend,xlen;
  int irange,jrange,krange;
  int fnumx = fnx + 2*ghx;
  int fnumy = fny + 2*ghy;
  int hnumx = hnx + 2*ghx;
  int hnumy = hny + 2*ghy;
  
  int swindex[bz+2*hz];
  const int ws  = bz+2*hz;
  
  int id,rid,cid,wid,eid,sid,nid;
  int ib,jb,kb,i,j,k;
  int dmasize;
  
  CalcID(&id, &rid, &cid, &wid, &eid, &sid, &nid, mx, my);
  CalcRange(id,rid,cid,lx,ly,lz,mx,my,mz,&xbeg,&xend,&xlen,&ybeg,&yend,&ylen,&zbeg,&zend,&zlen);
  
  DMA_DECL_DESC(READ_WRITE);
  CRTS_ssync_array();
  
  
  for (jb = ybeg; jb < yend; jb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (kb= zbeg; kb < zend ; kb+= bz)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        krange = MIN(kb+bz,zend) - kb;
        dmasize = (irange + 2*hx);
        
        if (kb == zbeg){
          for (k = kb ; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            DMA_IREAD(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfu[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib), &tenddw[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[k - kb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfpt[k - kb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + bz) % ws;
          for (k = kb + hz; k < MIN(kb + bz, zend) + 2 * hz ; k++){
            int pos = swindex[k - kb];
            DMA_IREAD(_fluxfrho+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib), &tenddw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingrho+(k*fnumy*fnumx+jb*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(k*fnumy*fnumx+jb*fnumx+ib), &fluxfpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (k = hz ; k < hz + krange ; k++)
          for (i = hx ; i < hx + irange ; i++){
            tenddrho[k][i] = (-(fluxfrho[swindex[(k + 1)]][i] - fluxfrho[swindex[k]][i]) / 50.0);
            tenddu[k][i] = (-(fluxfu[swindex[(k + 1)]][i] - fluxfu[swindex[k]][i]) / 50.0);
            tenddw[k][i] = (-(fluxfw[swindex[(k + 1)]][i] - fluxfw[swindex[k]][i]) / 50.0);
            tenddw[k][i] = (tenddw[k][i] - (state_forcingrho[swindex[k]][i] * 9.8));
            tenddpt[k][i] = (-(fluxfpt[swindex[(k + 1)]][i] - fluxfpt[swindex[k]][i]) / 50.0);
          }
        
        for (k = kb + hz  ; k < MIN(kb+ bz, zend) + hz ; k++){
          DMA_IWRITE(_tenddrho+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddrho[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddu+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddu[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddw+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddw[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddpt+(k*fnumy*fnumx+jb*fnumx+ib+hx), &tenddpt[k - kb][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesX_0_22(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesX_0_22_info *data = (ComputeTendenciesX_0_22_info*)(_ptr);
  
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
  
  
  double *_state_forcingrho = data->state_forcingrho;
  __attribute__ ((aligned (32))) double state_forcingrho[8][260];
  
  double *_state_forcingu = data->state_forcingu;
  __attribute__ ((aligned (32))) double state_forcingu[8][260];
  
  double *_state_forcingw = data->state_forcingw;
  __attribute__ ((aligned (32))) double state_forcingw[8][260];
  
  double *_state_forcingpt = data->state_forcingpt;
  __attribute__ ((aligned (32))) double state_forcingpt[8][260];
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][260];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][260];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][260];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][260];
  
  double *_hy_dens_cell = data->hy_dens_cell;
  __attribute__ ((aligned (32))) double hy_dens_cell[fnz+2*hz];
  
  double *_hy_dens_theta_cell = data->hy_dens_theta_cell;
  __attribute__ ((aligned (32))) double hy_dens_theta_cell[fnz+2*hz];
  
  struct Vector4 vals;
  struct Vector4 d3_vals;
  double r;
  double u;
  double w;
  double t;
  double p;
  
  double hv_coef = data->hv_coef;
  
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
  
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_cell+(ghz-hz), &hy_dens_cell[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_theta_cell+(ghz-hz), &hy_dens_theta_cell[0], dmasize* sizeof(double));
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_state_forcingrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingpt[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_forcingrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingu+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingw+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_forcingpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_forcingpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            vals.r = ((((-state_forcingrho[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingrho[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingrho[swindex[j]][i]) / 12)) - (state_forcingrho[swindex[j]][(i + 1)] / 12));
            d3_vals.r = (((-state_forcingrho[swindex[j]][(i - 2)] + (3 * state_forcingrho[swindex[j]][(i - 1)])) - (3 * state_forcingrho[swindex[j]][i])) + state_forcingrho[swindex[j]][(i + 1)]);
            vals.u = ((((-state_forcingu[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingu[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingu[swindex[j]][i]) / 12)) - (state_forcingu[swindex[j]][(i + 1)] / 12));
            d3_vals.u = (((-state_forcingu[swindex[j]][(i - 2)] + (3 * state_forcingu[swindex[j]][(i - 1)])) - (3 * state_forcingu[swindex[j]][i])) + state_forcingu[swindex[j]][(i + 1)]);
            vals.w = ((((-state_forcingw[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingw[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingw[swindex[j]][i]) / 12)) - (state_forcingw[swindex[j]][(i + 1)] / 12));
            d3_vals.w = (((-state_forcingw[swindex[j]][(i - 2)] + (3 * state_forcingw[swindex[j]][(i - 1)])) - (3 * state_forcingw[swindex[j]][i])) + state_forcingw[swindex[j]][(i + 1)]);
            vals.t = ((((-state_forcingpt[swindex[j]][(i - 2)] / 12) + ((7 * state_forcingpt[swindex[j]][(i - 1)]) / 12)) + ((7 * state_forcingpt[swindex[j]][i]) / 12)) - (state_forcingpt[swindex[j]][(i + 1)] / 12));
            d3_vals.t = (((-state_forcingpt[swindex[j]][(i - 2)] + (3 * state_forcingpt[swindex[j]][(i - 1)])) - (3 * state_forcingpt[swindex[j]][i])) + state_forcingpt[swindex[j]][(i + 1)]);
            r = (vals.r + hy_dens_cell[kb+oz]);
            u = (vals.u / r);
            w = (vals.w / r);
            t = ((vals.t + hy_dens_theta_cell[kb+oz]) / r);
            p = (27.562941092972594 * pow((r * t),1.400278940027894));
            fluxfrho[j][i] = ((r * u) - (hv_coef * d3_vals.r));
            fluxfu[j][i] = ((((r * u) * u) + p) - (hv_coef * d3_vals.u));
            fluxfw[j][i] = (((r * u) * w) - (hv_coef * d3_vals.w));
            fluxfpt[j][i] = (((r * u) * t) - (hv_coef * d3_vals.t));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfrho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &fluxfpt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void ComputeTendenciesX_0_23(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  ComputeTendenciesX_0_23_info *data = (ComputeTendenciesX_0_23_info*)(_ptr);
  
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
  
  
  double *_fluxfrho = data->fluxfrho;
  __attribute__ ((aligned (32))) double fluxfrho[8][258];
  
  double *_fluxfu = data->fluxfu;
  __attribute__ ((aligned (32))) double fluxfu[8][258];
  
  double *_fluxfw = data->fluxfw;
  __attribute__ ((aligned (32))) double fluxfw[8][258];
  
  double *_fluxfpt = data->fluxfpt;
  __attribute__ ((aligned (32))) double fluxfpt[8][258];
  
  double *_tenddrho = data->tenddrho;
  __attribute__ ((aligned (32))) double tenddrho[8][258];
  
  double *_tenddu = data->tenddu;
  __attribute__ ((aligned (32))) double tenddu[8][258];
  
  double *_tenddw = data->tenddw;
  __attribute__ ((aligned (32))) double tenddw[8][258];
  
  double *_tenddpt = data->tenddpt;
  __attribute__ ((aligned (32))) double tenddpt[8][258];
  
  
  
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
            DMA_IREAD(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfpt[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_fluxfrho+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfu+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfw+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_fluxfpt+(kb*fnumy*fnumx+j*fnumx+ib), &fluxfpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            tenddrho[j][i] = (-(fluxfrho[swindex[j]][(i + 1)] - fluxfrho[swindex[j]][i]) / 50.0);
            tenddu[j][i] = (-(fluxfu[swindex[j]][(i + 1)] - fluxfu[swindex[j]][i]) / 50.0);
            tenddw[j][i] = (-(fluxfw[swindex[j]][(i + 1)] - fluxfw[swindex[j]][i]) / 50.0);
            tenddpt[j][i] = (-(fluxfpt[swindex[j]][(i + 1)] - fluxfpt[swindex[j]][i]) / 50.0);
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_tenddrho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddrho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_tenddpt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &tenddpt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}
void UpdateState_0_5(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  UpdateState_0_5_info *data = (UpdateState_0_5_info*)(_ptr);
  
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
  
  
  double *_state_initrho = data->state_initrho;
  __attribute__ ((aligned (32))) double state_initrho[8][260];
  
  double *_tenddrho = data->tenddrho;
  __attribute__ ((aligned (32))) double tenddrho[8][260];
  
  double *_state_initu = data->state_initu;
  __attribute__ ((aligned (32))) double state_initu[8][260];
  
  double *_tenddu = data->tenddu;
  __attribute__ ((aligned (32))) double tenddu[8][260];
  
  double *_state_initw = data->state_initw;
  __attribute__ ((aligned (32))) double state_initw[8][260];
  
  double *_tenddw = data->tenddw;
  __attribute__ ((aligned (32))) double tenddw[8][260];
  
  double *_state_initpt = data->state_initpt;
  __attribute__ ((aligned (32))) double state_initpt[8][260];
  
  double *_tenddpt = data->tenddpt;
  __attribute__ ((aligned (32))) double tenddpt[8][260];
  
  double *_state_outrho = data->state_outrho;
  __attribute__ ((aligned (32))) double state_outrho[8][260];
  
  double *_state_outu = data->state_outu;
  __attribute__ ((aligned (32))) double state_outu[8][260];
  
  double *_state_outw = data->state_outw;
  __attribute__ ((aligned (32))) double state_outw[8][260];
  
  double *_state_outpt = data->state_outpt;
  __attribute__ ((aligned (32))) double state_outpt[8][260];
  
  
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
            DMA_IREAD(_state_initrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_initrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddrho+(kb*fnumy*fnumx+j*fnumx+ib), &tenddrho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initu+(kb*fnumy*fnumx+j*fnumx+ib), &state_initu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddu+(kb*fnumy*fnumx+j*fnumx+ib), &tenddu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initw+(kb*fnumy*fnumx+j*fnumx+ib), &state_initw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(kb*fnumy*fnumx+j*fnumx+ib), &tenddw[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_initpt[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddpt+(kb*fnumy*fnumx+j*fnumx+ib), &tenddpt[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_state_initrho+(kb*fnumy*fnumx+j*fnumx+ib), &state_initrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddrho+(kb*fnumy*fnumx+j*fnumx+ib), &tenddrho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initu+(kb*fnumy*fnumx+j*fnumx+ib), &state_initu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddu+(kb*fnumy*fnumx+j*fnumx+ib), &tenddu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initw+(kb*fnumy*fnumx+j*fnumx+ib), &state_initw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddw+(kb*fnumy*fnumx+j*fnumx+ib), &tenddw[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_state_initpt+(kb*fnumy*fnumx+j*fnumx+ib), &state_initpt[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_tenddpt+(kb*fnumy*fnumx+j*fnumx+ib), &tenddpt[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            state_outrho[j][i] = (state_initrho[swindex[j]][i] + (dt * tenddrho[swindex[j]][i]));
            state_outu[j][i] = (state_initu[swindex[j]][i] + (dt * tenddu[swindex[j]][i]));
            state_outw[j][i] = (state_initw[swindex[j]][i] + (dt * tenddw[swindex[j]][i]));
            state_outpt[j][i] = (state_initpt[swindex[j]][i] + (dt * tenddpt[swindex[j]][i]));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_state_outrho+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_outrho[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_outu+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_outu[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_outw+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_outw[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_state_outpt+(kb*fnumy*fnumx+j*fnumx+ib+hx), &state_outpt[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}


void OutputPrepare_0_0(void *_ptr){
  
  volatile int reply = 0;
  volatile int COUNT = 0;
  
  OutputPrepare_0_0_info *data = (OutputPrepare_0_0_info*)(_ptr);
  
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
  
  
  double *_staterho = data->staterho;
  __attribute__ ((aligned (32))) double staterho[8][260];
  
  double *_stateu = data->stateu;
  __attribute__ ((aligned (32))) double stateu[8][260];
  
  double *_statew = data->statew;
  __attribute__ ((aligned (32))) double statew[8][260];
  
  double *_statept = data->statept;
  __attribute__ ((aligned (32))) double statept[8][260];
  
  double *_stateoutdens = data->stateoutdens;
  __attribute__ ((aligned (32))) double stateoutdens[8][260];
  
  double *_stateoutuwnd = data->stateoutuwnd;
  __attribute__ ((aligned (32))) double stateoutuwnd[8][260];
  
  double *_stateoutwwnd = data->stateoutwwnd;
  __attribute__ ((aligned (32))) double stateoutwwnd[8][260];
  
  double *_stateouttheta = data->stateouttheta;
  __attribute__ ((aligned (32))) double stateouttheta[8][260];
  
  double *_hy_dens_cell = data->hy_dens_cell;
  __attribute__ ((aligned (32))) double hy_dens_cell[fnz+2*hz];
  
  double *_hy_dens_theta_cell = data->hy_dens_theta_cell;
  __attribute__ ((aligned (32))) double hy_dens_theta_cell[fnz+2*hz];
  
  
  
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
  
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_cell+(ghz-hz), &hy_dens_cell[0], dmasize* sizeof(double));
  dmasize = fnz + 2*hz;
  DMA_READ(_hy_dens_theta_cell+(ghz-hz), &hy_dens_theta_cell[0], dmasize* sizeof(double));
  
  for (kb = zbeg; kb < zend; kb++){
    for (i = 0; i < ws; i++) swindex[i] = i;
    
    for (jb= ybeg; jb < yend ; jb+= by)
      for (ib = xbeg; ib < xend; ib+= bx){
        irange = MIN(ib+bx,xend) - ib;
        jrange = MIN(jb+by,yend) - jb;
        dmasize = (irange + 2*hx);
        
        if (jb == ybeg){
          for (j = jb ; j < MIN(jb + by, yend) + 2 * hy ; j++){
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[j - jb][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[j - jb][0], dmasize* sizeof(double));
          }
        }
        else{
          for (i = 0; i < ws; i++)
            swindex[i] = (swindex[i] + by) % ws;
          for (j = jb + hy; j < MIN(jb + by, yend) + 2 * hy ; j++){
            int pos = swindex[j - jb];
            DMA_IREAD(_staterho+(kb*fnumy*fnumx+j*fnumx+ib), &staterho[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_stateu+(kb*fnumy*fnumx+j*fnumx+ib), &stateu[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statew+(kb*fnumy*fnumx+j*fnumx+ib), &statew[pos][0], dmasize* sizeof(double));
            DMA_IREAD(_statept+(kb*fnumy*fnumx+j*fnumx+ib), &statept[pos][0], dmasize* sizeof(double));
          }
        }
        
        DMA_WAIT_ALL;
        
        for (j = hy ; j < hy + jrange ; j++)
          for (i = hx ; i < hx + irange ; i++){
            stateoutdens[j][i] = staterho[swindex[j]][i];
            stateoutuwnd[j][i] = (stateu[swindex[j]][i] / (hy_dens_cell[kb+oz] + staterho[swindex[j]][i]));
            stateoutwwnd[j][i] = (statew[swindex[j]][i] / (hy_dens_cell[kb+oz] + staterho[swindex[j]][i]));
            stateouttheta[j][i] = (((statept[swindex[j]][i] + hy_dens_theta_cell[kb+oz]) / (hy_dens_cell[kb+oz] + staterho[swindex[j]][i])) - (hy_dens_theta_cell[kb+oz] / hy_dens_cell[kb+oz]));
          }
        
        for (j = jb + hy  ; j < MIN(jb+ by, yend) + hy ; j++){
          DMA_IWRITE(_stateoutdens+(kb*fnumy*fnumx+j*fnumx+ib+hx), &stateoutdens[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_stateoutuwnd+(kb*fnumy*fnumx+j*fnumx+ib+hx), &stateoutuwnd[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_stateoutwwnd+(kb*fnumy*fnumx+j*fnumx+ib+hx), &stateoutwwnd[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
          DMA_IWRITE(_stateouttheta+(kb*fnumy*fnumx+j*fnumx+ib+hx), &stateouttheta[j - jb ][hx], (dmasize - 2*hx) * sizeof(double));
        }
      }
  }
}

