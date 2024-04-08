#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

#include"../ExternLIB/filter.h"
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

struct Vector3 cartesian_transform(double lon, double lat)
{
  double cos_lat, x, y, z;
  cos_lat = cos(lat);
  x = ((6371220.0 * cos_lat) * cos(lon));
  y = ((6371220.0 * cos_lat) * sin(lon));
  z = (6371220.0 * sin(lat));
  return (struct Vector3){x,y,z};
}

struct Vector2 inverse_cartesian_transform(double x, double y, double z)
{
  double lon, lat;
  lon = atan2(y,x);
  lat = asin((z / 6371220.0));
  if ((lon < 0.0))
  {
    lon = (lon + 6.283185307179586);
  }
  return (struct Vector2){lon,lat};
}

struct Vector3 cross_product(struct Vector3* a, struct Vector3* b)
{
  struct Vector3 res;
  res = (struct Vector3){0,0,0};
  res.x = (((a->y) * (b->z)) - ((a->z) * (b->y)));
  res.y = (((a->z) * (b->x)) - ((a->x) * (b->z)));
  res.z = (((a->x) * (b->y)) - ((a->y) * (b->x)));
  return res;
}

struct Vector3 norm_vector(struct Vector3* v)
{
  struct Vector3 res;
  double sum;
  res = (struct Vector3){0,0,0};
  sum = ((((v->x) * (v->x)) + ((v->y) * (v->y))) + ((v->z) * (v->z)));
  sum = sqrt(sum);
  if ((sum != 0))
  {
    res.x = ((v->x) / sum);
    res.y = ((v->y) / sum);
    res.z = ((v->z) / sum);
  }
  return res;
}

double calc_sphere_angle(struct Vector3* a, struct Vector3* b, struct Vector3* c)
{
  struct Vector3 nab, nbc, ta, tb;
  double dot, res;
  nab = (struct Vector3){0,0,0};
  nbc = (struct Vector3){0,0,0};
  ta = (struct Vector3){0,0,0};
  tb = (struct Vector3){0,0,0};
  ta = cross_product(a,b);
  tb = cross_product(b,c);
  nab = norm_vector(&ta);
  nbc = norm_vector(&tb);
  dot = (((nab.x * nbc.x) + (nab.y * nbc.y)) + (nab.z * nbc.z));
  res = acos(-MAX(MIN(dot,1.0),-1.0));
  return res;
}

double calc_area(struct Vector3* x, struct Vector3* y, struct Vector3* z)
{
  int n;
  double res, angle;
  struct Vector3 vi, vip, vim;
  n = 3;
  res = 0.0;
  angle = 0.0;
  vi = (struct Vector3){(x->x),(y->x),(z->x)};
  vip = (struct Vector3){(x->y),(y->y),(z->y)};
  vim = (struct Vector3){(x->z),(y->z),(z->z)};
  angle = calc_sphere_angle(&vim,&vi,&vip);
  res = (res + angle);
  angle = calc_sphere_angle(&vi,&vip,&vim);
  res = (res + angle);
  angle = calc_sphere_angle(&vip,&vim,&vi);
  res = (res + angle);
  res = (pow(6371220.0,2) * (res - ((n - 2) * 3.141592653589793)));
  return res;
}

double calc_area_with_last_small_arc(struct Vector3* x, struct Vector3* y, struct Vector3* z)
{
  struct Vector3 xv, yv, zv;
  double res, lat0, lon1, lat1, lon2, dlon, area1, area2, area3;
  struct Vector2 lonlat0, lonlat1, lonlat2;
  xv = (struct Vector3){0,0,0};
  yv = (struct Vector3){0,0,0};
  zv = (struct Vector3){0,0,0};
  res = calc_area(x,y,z);
  lonlat0 = inverse_cartesian_transform((x->x),(y->x),(z->x));
  lat0 = lonlat0.y;
  lonlat1 = inverse_cartesian_transform((x->y),(y->y),(z->y));
  lon1 = lonlat1.x;
  lat1 = lonlat1.y;
  lonlat2 = inverse_cartesian_transform((x->z),(y->z),(z->z));
  lon2 = lonlat2.x;
  if ((lat1 == 0))
  {
    return res;
  }
  else
  {
    if ((lat0 > lat1))
    {
      dlon = (lon2 - lon1);
    }
    else
    {
      dlon = (lon1 - lon2);
    }
    if ((dlon < 0.0))
    {
      dlon = (dlon + 6.283185307179586);
    }
    xv.x = 0.0;
    yv.x = 0.0;
    if ((((lat0 * lat1) >= 0) && (fabs(lat0) > fabs(lat1))))
    {
      xv.y = (x->y);
      yv.y = (y->y);
      zv.y = (z->y);
      xv.z = (x->z);
      yv.z = (y->z);
      zv.z = (z->z);
    }
    else
    {
      xv.y = (x->z);
      yv.y = (y->z);
      zv.y = (z->z);
      xv.z = (x->y);
      yv.z = (y->y);
      zv.z = (z->y);
    }
    if ((lat1 > 0.0))
    {
      zv.x = 6371220.0;
      area1 = ((pow(6371220.0,2) * dlon) * (1.0 - sin(lat1)));
    }
    else
    {
      zv.x = -6371220.0;
      area1 = ((pow(6371220.0,2) * dlon) * (sin(lat1) + 1.0));
    }
    area2 = calc_area(&xv,&yv,&zv);
    area3 = (area1 - area2);
    if (((area3 < 0.0) && (fabs(area3) > 1e-10)))
    {
      area3 = 0;
    }
    if ((((lat0 * lat1) >= 0) && (fabs(lat0) > fabs(lat1))))
    {
      res = (res + area3);
    }
    else
    {
      res = (res - area3);
    }
    return res;
  }
}

double hybrid_coord_calc_ph(double hyam, double hybm, double phs)
{
  return (((hyam * (100000.0 - 219.4)) + (hybm * (phs - 219.4))) + 219.4);
}

double hybrid_coord_calc_ph_lev(double hyai, double hybi, double phs)
{
  return (((hyai * (100000.0 - 219.4)) + (hybi * (phs - 219.4))) + 219.4);
}

double hybrid_coord_calc_dphdt_lev(double hybi, double dphsdt)
{
  return (hybi * dphsdt);
}

double specific_humidity(double qv)
{
  return (qv / (1 + qv));
}

double mixing_ratio(double sh)
{
  return (sh / (1 - sh));
}

double sign(double x, double y)
{
  if ((y > 0.0))
  {
    return x;
  }
  else
  {
    return -x;
  }
}

double potential_temperature(double t, double p, double qv)
{
  return ((t * pow((100000.0 / p),0.2858964143426295)) * (1 + (1.607779403567447 * qv)));
}

double temperature(double pt, double p, double qv)
{
  return ((pt * pow((p / 100000.0),0.2858964143426295)) / (1 + (1.607779403567447 * qv)));
}

double virtual_temperature(double t, double sh)
{
  return (t * (1 + ((1.607779403567447 - 1) * sh)));
}

double dry_air_density(double pt, double p)
{
  return (((100000.0 / 287.04) / pt) * pow((p / 100000.0),0.7141434262948207));
}

double moist_air_density(double t, double p, double qv)
{
  return ((p / 287.04) / virtual_temperature(t,specific_humidity(qv)));
}

double upwind1(double dir, double wgt, double f1, double f2)
{
  return ((0.5 * (f2 + f1)) - (((0.5 * (f2 - f1)) * wgt) * dir));
}

double upwind3(double dir, double wgt, double f1, double f2, double f3, double f4)
{
  double c31, c32, c33;
  c31 = (7.0 / 12.0);
  c32 = (-1.0 / 12.0);
  c33 = (1.0 / 12.0);
  return (((c31 * (f3 + f2)) + (c32 * (f4 + f1))) + (((c33 * ((f4 - f1) - (3 * (f3 - f2)))) * wgt) * dir));
}

double slope(double fm1, double f, double fp1)
{
  double df, df_min, df_max;
  df = ((fp1 - fm1) * 0.5);
  df_min = (2 * (f - MIN(MIN(fm1,f),fp1)));
  df_max = (2 * (MAX(MAX(fm1,f),fp1) - f));
  return sign(MIN(fabs(df),MIN(df_min,df_max)),df);
}

struct Vector3 ppm(double fm2, double fm1, double f, double fp1, double fp2)
{
  double dfl, df, dfr, fl, fr, f6;
  dfl = slope(fm2,fm1,f);
  df = slope(fm1,f,fp1);
  dfr = slope(f,fp1,fp2);
  fl = ((0.5 * (fm1 + f)) + ((dfl - df) / 6.0));
  fr = ((0.5 * (fp1 + f)) + ((df - dfr) / 6.0));
  fl = (f - sign(MIN(fabs(df),fabs((fl - f))),df));
  fr = (f + sign(MIN(fabs(df),fabs((fr - f))),df));
  f6 = ((6 * f) - (3 * (fl + fr)));
  df = (fr - fl);
  return (struct Vector3){fl,df,f6};
}

void LatLonMeshInit(struct HybridMeshField* mesh)
{
  int t, i, j, k;
  struct MeshVector6 loopinfo;
  double dlon, dlat0;
  struct MeshVector3 v_inx;
  struct Vector3 v1, v2, v3, buf;
  t = 0;
  loopinfo = (struct MeshVector6){1,1,2400+2*Proc.lon_hw,0,0,-Proc.lon_hw + 1};
  dlon = (double)(((6.283185307179586 - 0) / 2400));
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=0; i<2400+2*Proc.lon_hw; i+=1){
    v_inx = disect(t,&loopinfo);
    t = (t + 1);
    mesh->full_lon[0][0][i] = (0 + ((v_inx.z - 1) * dlon));
    mesh->half_lon[0][0][i] = (mesh->full_lon[0][0][i] + (0.5 * dlon));
    mesh->full_lon_deg[0][0][i] = (mesh->full_lon[0][0][i] * 57.29577951308232);
    mesh->half_lon_deg[0][0][i] = (mesh->half_lon[0][0][i] * 57.29577951308232);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dlat0 = (double)(((1.5707963267948966 - -1.5707963267948966) / 1199));
  t = 0;
  loopinfo = (struct MeshVector6){1,1199,1,0,1,0};
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    v_inx = disect(t,&loopinfo);
    t = (t + 1);
    mesh->half_lat[0][j][0] = (-1.5707963267948966 + ((v_inx.y - 0.5) * dlat0));
    if ((fabs(mesh->half_lat[0][j][0]) < 1e-12))
    {
      mesh->half_lat[0][j][0] = 0.0;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->dlat[0][j][0] = dlat0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1+Proc.lat_hw; j+=1){
    mesh->full_lat[0][j][0] = -1.5707963267948966;
    mesh->full_lat_deg[0][j][0] = (-1.5707963267948966 * 57.29577951308232);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->full_lat[0][j][0] = (mesh->full_lat[0][(j - 1)][0] + mesh->dlat[0][(j - 1)][0]);
    if ((fabs(mesh->full_lat[0][j][0]) < 1e-12))
    {
      mesh->full_lat[0][j][0] = 0.0;
    }
    mesh->full_lat_deg[0][j][0] = (mesh->full_lat[0][j][0] * 57.29577951308232);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1199+Proc.lat_hw; j<1200+Proc.lat_hw; j+=1){
    mesh->full_lat[0][j][0] = 1.5707963267948966;
    mesh->full_lat_deg[0][j][0] = (1.5707963267948966 * 57.29577951308232);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    if (!(mesh->full_lat[0][j][0] == 1.5707963267948966))
    {
      mesh->half_lat[0][j][0] = (mesh->full_lat[0][j][0] + (0.5 * mesh->dlat[0][j][0]));
      if ((fabs(mesh->half_lat[0][j][0]) < 1e-12))
      {
        mesh->half_lat[0][j][0] = 0.0;
      }
      mesh->half_lat_deg[0][j][0] = (mesh->half_lat[0][j][0] * 57.29577951308232);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=0; i<2400+2*Proc.lon_hw; i+=1){
    mesh->full_cos_lon[0][0][i] = cos(mesh->full_lon[0][0][i]);
    mesh->full_sin_lon[0][0][i] = sin(mesh->full_lon[0][0][i]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(i=0; i<2401+2*Proc.lon_hw; i+=1){
    mesh->half_cos_lon[0][0][i] = cos(mesh->half_lon[0][0][i]);
    mesh->half_sin_lon[0][0][i] = sin(mesh->half_lon[0][0][i]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0; j<1199+2*Proc.lat_hw; j+=1){
    if (((mesh->half_lat[0][j][0] >= -1.5707963267948966) && (mesh->half_lat[0][j][0] <= 1.5707963267948966)))
    {
      mesh->half_cos_lat[0][j][0] = cos(mesh->half_lat[0][j][0]);
      mesh->half_sin_lat[0][j][0] = sin(mesh->half_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0; j<1200+2*Proc.lat_hw; j+=1){
    if (((mesh->full_lat[0][j][0] >= -1.5707963267948966) && (mesh->full_lat[0][j][0] <= 1.5707963267948966)))
    {
      mesh->full_cos_lat[0][j][0] = cos(mesh->full_lat[0][j][0]);
      mesh->full_sin_lat[0][j][0] = sin(mesh->full_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1+Proc.lat_hw; j+=1){
    mesh->full_cos_lat[0][j][0] = 0.0;
    mesh->full_sin_lat[0][j][0] = -1.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1199+Proc.lat_hw; j<1200+Proc.lat_hw; j+=1){
    mesh->full_cos_lat[0][j][0] = 0.0;
    mesh->full_sin_lat[0][j][0] = 1.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1+Proc.lat_hw; j+=1){
    mesh->area_cell[0][j][0] = ((pow(6371220.0,2) * dlon) * (mesh->half_sin_lat[0][j][0] + 1.0));
    mesh->area_subcell_1[0][j][0] = (((pow(6371220.0,2) * 0.5) * dlon) * (mesh->half_sin_lat[0][j][0] + 1.0));
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1199+Proc.lat_hw; j<1200+Proc.lat_hw; j+=1){
    mesh->area_cell[0][j][0] = ((pow(6371220.0,2) * dlon) * (1.0 - mesh->half_sin_lat[0][(j - 1)][0]));
    mesh->area_subcell_0[0][j][0] = (((pow(6371220.0,2) * 0.5) * dlon) * (1.0 - mesh->half_sin_lat[0][(j - 1)][0]));
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  v1 = (struct Vector3){0.0,0.0,0.0};
  v2 = (struct Vector3){0.0,0.0,0.0};
  v3 = (struct Vector3){0.0,0.0,0.0};
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    for(i=0+Proc.lon_hw; i<1+Proc.lon_hw; i+=1){
      mesh->area_cell[0][j][0] = ((pow(6371220.0,2) * dlon) * (mesh->half_sin_lat[0][j][0] - mesh->half_sin_lat[0][(j - 1)][0]));
      mesh->area_subcell_0[0][j][0] = (((pow(6371220.0,2) * 0.5) * dlon) * (mesh->full_sin_lat[0][j][0] - mesh->half_sin_lat[0][(j - 1)][0]));
      mesh->area_subcell_1[0][j][0] = (((pow(6371220.0,2) * 0.5) * dlon) * (mesh->half_sin_lat[0][j][0] - mesh->full_sin_lat[0][j][0]));
      buf = cartesian_transform(mesh->full_lon[0][0][i],mesh->full_lat[0][j][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][(j - 1)][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][j][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lon_west[0][j][0] = calc_area(&v1,&v2,&v3);
      mesh->area_lon_east[0][j][0] = mesh->area_lon_west[0][j][0];
      mesh->area_lon[0][j][0] = (mesh->area_lon_west[0][j][0] + mesh->area_lon_east[0][j][0]);
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][j][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][i],mesh->full_lat[0][j][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][j][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lon_north[0][j][0] = calc_area_with_last_small_arc(&v1,&v2,&v3);
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][(j - 1)][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][j][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][i],mesh->full_lat[0][j][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lon_south[0][j][0] = calc_area_with_last_small_arc(&v1,&v2,&v3);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    for(i=0+Proc.lon_hw; i<1+Proc.lon_hw; i+=1){
      mesh->area_vtx[0][j][0] = ((pow(6371220.0,2) * dlon) * (mesh->full_sin_lat[0][(j + 1)][0] - mesh->full_sin_lat[0][j][0]));
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][j][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][j][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][(j + 1)][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lat_west[0][j][0] = calc_area(&v1,&v2,&v3);
      mesh->area_lat_east[0][j][0] = mesh->area_lat_west[0][j][0];
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1198+Proc.lat_hw; j+=1){
    for(i=0+Proc.lon_hw; i<1+Proc.lon_hw; i+=1){
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][(j + 1)][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][j][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][(i + 1)],mesh->half_lat[0][j][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lat_north[0][j][0] = calc_area_with_last_small_arc(&v1,&v2,&v3);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1198+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->area_lat_north[0][j][0] = mesh->area_cell[0][(j + 1)][0];
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1+Proc.lat_hw; j+=1){
    mesh->area_lat_south[0][j][0] = mesh->area_cell[0][j][0];
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    for(i=0+Proc.lon_hw; i<1+Proc.lon_hw; i+=1){
      buf = cartesian_transform(mesh->full_lon[0][0][(i + 1)],mesh->full_lat[0][j][0]);
      v1.x = buf.x;
      v2.x = buf.y;
      v3.x = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][(i + 1)],mesh->half_lat[0][j][0]);
      v1.y = buf.x;
      v2.y = buf.y;
      v3.y = buf.z;
      buf = cartesian_transform(mesh->half_lon[0][0][i],mesh->half_lat[0][j][0]);
      v1.z = buf.x;
      v2.z = buf.y;
      v3.z = buf.z;
      mesh->area_lat_south[0][j][0] = calc_area_with_last_small_arc(&v1,&v2,&v3);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->area_lat[0][j][0] = (mesh->area_lat_north[0][j][0] + mesh->area_lat_south[0][j][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->de_lon[0][j][0] = ((6371220.0 * mesh->full_cos_lat[0][j][0]) * dlon);
    mesh->le_lon[0][j][0] = ((2.0 * mesh->area_lon[0][j][0]) / mesh->de_lon[0][j][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1+Proc.lat_hw; j+=1){
    mesh->le_lon[0][j][0] = 0.0;
    mesh->de_lon[0][j][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1199+Proc.lat_hw; j<1200+Proc.lat_hw; j+=1){
    mesh->le_lon[0][j][0] = 0.0;
    mesh->de_lon[0][j][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->le_lat[0][j][0] = ((6371220.0 * mesh->half_cos_lat[0][j][0]) * dlon);
    mesh->de_lat[0][j][0] = ((2.0 * mesh->area_lat[0][j][0]) / mesh->le_lat[0][j][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->full_tangent_wgt_0[0][j][0] = ((mesh->le_lat[0][(j - 1)][0] / mesh->de_lon[0][j][0]) * 0.25);
    mesh->full_tangent_wgt_1[0][j][0] = ((mesh->le_lat[0][j][0] / mesh->de_lon[0][j][0]) * 0.25);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->half_tangent_wgt_0[0][j][0] = ((mesh->le_lon[0][j][0] / mesh->de_lat[0][j][0]) * 0.25);
    mesh->half_tangent_wgt_1[0][j][0] = ((mesh->le_lon[0][(j + 1)][0] / mesh->de_lat[0][j][0]) * 0.25);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1200+Proc.lat_hw; j+=1){
    mesh->full_f[0][j][0] = ((2.0 * 7.27220521664304e-05) * mesh->full_sin_lat[0][j][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
    mesh->half_f[0][j][0] = ((2.0 * 7.27220521664304e-05) * mesh->half_sin_lat[0][j][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<2+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0027008056640625;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=2+Proc.lev_hw; k<3+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.007686834782362;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=3+Proc.lev_hw; k<4+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0158526937798342;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=4+Proc.lev_hw; k<5+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0276367470939261;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=5+Proc.lev_hw; k<6+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0414515552099706;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=6+Proc.lev_hw; k<7+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0562277844155789;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=7+Proc.lev_hw; k<8+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0724737972365225;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=8+Proc.lev_hw; k<9+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0893795253309864;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=9+Proc.lev_hw; k<10+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1063548465048298;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=10+Proc.lev_hw; k<11+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.125472525784964;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=11+Proc.lev_hw; k<12+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1479642011314251;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=12+Proc.lev_hw; k<13+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1744233770857367;
    mesh->hybi[k][0][0] = 0.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=13+Proc.lev_hw; k<14+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2054476619883455;
    mesh->hybi[k][0][0] = 0.0001056069805219;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=14+Proc.lev_hw; k<15+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2362396334209435;
    mesh->hybi[k][0][0] = 0.0059397906834944;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=15+Proc.lev_hw; k<16+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2617261684535087;
    mesh->hybi[k][0][0] = 0.0235374270991938;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=16+Proc.lev_hw; k<17+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2783057827853432;
    mesh->hybi[k][0][0] = 0.0576439176162314;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=17+Proc.lev_hw; k<18+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2816243248017372;
    mesh->hybi[k][0][0] = 0.1139578493851603;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=18+Proc.lev_hw; k<19+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2683280168286641;
    mesh->hybi[k][0][0] = 0.1932313936022966;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=19+Proc.lev_hw; k<20+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2429473354807992;
    mesh->hybi[k][0][0] = 0.2808764820852269;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=20+Proc.lev_hw; k<21+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.2108099294331058;
    mesh->hybi[k][0][0] = 0.3715629284436812;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=21+Proc.lev_hw; k<22+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1759998110016148;
    mesh->hybi[k][0][0] = 0.4612091687348547;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=22+Proc.lev_hw; k<23+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1415334037742616;
    mesh->hybi[k][0][0] = 0.5467954414732804;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=23+Proc.lev_hw; k<24+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.1095123294032484;
    mesh->hybi[k][0][0] = 0.6262219272258057;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=24+Proc.lev_hw; k<25+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0812761592574764;
    mesh->hybi[k][0][0] = 0.6981497395977654;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=25+Proc.lev_hw; k<26+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0575423961814687;
    mesh->hybi[k][0][0] = 0.7618575996945027;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=26+Proc.lev_hw; k<27+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0385285881906738;
    mesh->hybi[k][0][0] = 0.8171268253424389;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=27+Proc.lev_hw; k<28+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0240694984015316;
    mesh->hybi[k][0][0] = 0.8641240158900223;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=28+Proc.lev_hw; k<29+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0137226170696702;
    mesh->hybi[k][0][0] = 0.9032994538662682;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=29+Proc.lev_hw; k<30+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0068720184783044;
    mesh->hybi[k][0][0] = 0.935256044816174;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=30+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0027973066609278;
    mesh->hybi[k][0][0] = 0.9607227915792738;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0007580797514921;
    mesh->hybi[k][0][0] = 0.9804356127137464;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    mesh->hyai[k][0][0] = 0.0;
    mesh->hybi[k][0][0] = 1.0;
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    mesh->hyam[k][0][0] = (0.5 * (mesh->hyai[k][0][0] + mesh->hyai[(k + 1)][0][0]));
    mesh->hybm[k][0][0] = (0.5 * (mesh->hybi[k][0][0] + mesh->hybi[(k + 1)][0][0]));
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    mesh->full_lev[k][0][0] = (mesh->hyam[k][0][0] + mesh->hybm[k][0][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    mesh->half_lev[k][0][0] = (mesh->hyai[k][0][0] + mesh->hybi[k][0][0]);
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=1+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
      mesh->c_lon[k][j][0] = (((0.0078125 * MAX(1.0,(8 * (1 + tanh(log((219.4 / hybrid_coord_calc_ph(mesh->hyam[k][0][0],mesh->hybm[k][0][0],100000.0)))))))) * mesh->le_lon[0][j][0]) * mesh->de_lon[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=0+Proc.lat_hw; j<1199+Proc.lat_hw; j+=1){
      mesh->c_lat[k][j][0] = (((0.0078125 * MAX(1.0,(8 * (1 + tanh(log((219.4 / hybrid_coord_calc_ph(mesh->hyam[k][0][0],mesh->hybm[k][0][0],100000.0)))))))) * mesh->le_lat[0][j][0]) * mesh->de_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  LatLonMeshInitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void ncIO_Read_0(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int nid;
  nid = ncOpenFile("gmcore.restart.DSL.15.nc",0);
  ncGetVar(nid, "u_lon", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.half_nlon, Proc.lon_hw, &(state->u_lon[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "v_lat", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.half_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->v_lat[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  
  ncGetVar(nid, "pt", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->pt[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "phs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->phs[0][0][0]));
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  
  ncGetVar(nid, "gzs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(staticv->gzs[0][0][0]));
  UpdateHalo_2d_D(Proc, &staticv->gzs[0][0][0], &Proc.FieldReq[async_gzs], true, true, true, true, true, true, false);
  
  ncCloseFile(nid);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ncCloseFileTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void spaceOperatorInit_0(struct HybridStaticField* staticv, struct HybridMeshField* mesh)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlon[0][j][i] = (((staticv->gzs[0][j][(i + 1)] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lon[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlat[0][j][i] = (((staticv->gzs[0][(j + 1)][i] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void gzlevInit_0(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->gz_lev[k][j][i] = staticv->gzs[0][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  preparegzlevgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->gz_lev[32 + Proc.lev_hw][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void uvInit_0(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void operatorPrepareNull_0(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  int k, j, i, l;
  double b, dgz;
  double *tmp_state_tmpsum;
  struct Vector4 ke_vtx;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->t[k][j][i] = temperature(state->pt[k][j][i],state->ph[k][j][i],state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state->tmpsum[k][j][i] = advpt->vv[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
    tmp_state_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_state_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divy[k][j][i] = (((state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lon[k][j][i] = (state->m_lon[k][j][i] * state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lat[k][j][i] = (state->m_lat[k][j][i] * state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (state->mfx_lon[k][j][(i - 1)] + state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (state->mfx_lon[k][(j + 1)][(i - 1)] + state->mfx_lon[k][(j + 1)][i])));
        state->u_lat[k][j][i] = (state->mfx_lat[k][j][i] / state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (state->mfy_lat[k][(j - 1)][i] + state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (state->mfy_lat[k][j][i] + state->mfy_lat[k][j][(i + 1)])));
        state->v_lon[k][j][i] = (state->mfy_lon[k][j][i] / state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = pow(state->v_lat[k][j][i],2);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ketmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = (state->tmpsum[k][j][i] / 2400);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->vor[k][j][i] = ((((state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = (-state->u_lat[k][j][i] * mesh->le_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vortmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->vor[k][j][i] = ((state->tmpsum[k][j][i] / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pv[k][j][i] = ((state->vor[k][j][i] + mesh->half_f[0][j][0]) / state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->v_lon[k][j][i]) / (sqrt((pow(state->u_lon[k][j][i],2) + pow(state->v_lon[k][j][i],2))) + 1e-24));
        state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,state->v_lon[k][j][i]),1,state->pv[k][(j - 2)][i],state->pv[k][(j - 1)][i],state->pv[k][j][i],state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (state->pv[k][(j - 1)][i] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->u_lat[k][j][i]) / (sqrt((pow(state->u_lat[k][j][i],2) + pow(state->v_lat[k][j][i],2))) + 1e-24));
        state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,state->u_lat[k][j][i]),1,state->pv[k][j][(i - 2)],state->pv[k][j][(i - 1)],state->pv[k][j][i],state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (state->pv[k][j][(i - 1)] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->div[k][j][i] = ((((state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = state->v_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->div[k][j][i] = (((state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * state->t[l][j][i]) * log((state->ph_lev[(l + 1)][j][i] / state->ph_lev[l][j][i]))));
        }
        state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void advPrepare_0(struct HybridStateField* state, struct HybridAdvField* advm, struct HybridAdvField* advpt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advm->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void filterinit_0(struct HybridMeshField* mesh, double dt)
{
  Filter_Init(mesh,dt);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_InitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void stepForwardBackward_0(struct HybridStateField* old_state, struct HybridStateField* star_state, struct HybridStateField* new_state, struct HybridStaticField* staticv, struct HybridTendField* tend1, struct HybridTendField* tend2, struct HybridTendPara* tendPara, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double b, mfs, cursum, cf, s1, s2, ds1, ds2, ds3, dgz, tl, dph1, dph2, dgz1, dgz2;
  int k, j, i, l, ci;
  double *tmp_star_state_tmpsum;
  struct Vector4 ke_vtx;
  struct Vector3 buf;
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = star_state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = star_state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + star_state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + star_state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + star_state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + star_state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          star_state->tmpsum[k][j][i] = advpt->vv[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
    tmp_star_state_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_star_state_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divy[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lon[k][j][i] = (star_state->m_lon[k][j][i] * star_state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lat[k][j][i] = (star_state->m_lat[k][j][i] * star_state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = star_state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = star_state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (star_state->mfx_lon[k][j][(i - 1)] + star_state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (star_state->mfx_lon[k][(j + 1)][(i - 1)] + star_state->mfx_lon[k][(j + 1)][i])));
        star_state->u_lat[k][j][i] = (star_state->mfx_lat[k][j][i] / star_state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (star_state->mfy_lat[k][(j - 1)][i] + star_state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (star_state->mfy_lat[k][j][i] + star_state->mfy_lat[k][j][(i + 1)])));
        star_state->v_lon[k][j][i] = (star_state->mfy_lon[k][j][i] / star_state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(star_state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(star_state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(star_state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * star_state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = pow(star_state->v_lat[k][j][i],2);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ketmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = (star_state->tmpsum[k][j][i] / 2400);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->div[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (star_state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((star_state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = star_state->v_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->div[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->vor[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (star_state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((star_state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (star_state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = (-star_state->u_lat[k][j][i] * mesh->le_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vortmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->vor[k][j][i] = ((star_state->tmpsum[k][j][i] / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pv[k][j][i] = ((star_state->vor[k][j][i] + mesh->half_f[0][j][0]) / star_state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->v_lon[k][j][i]) / (sqrt((pow(star_state->u_lon[k][j][i],2) + pow(star_state->v_lon[k][j][i],2))) + 1e-24));
        star_state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,star_state->v_lon[k][j][i]),1,star_state->pv[k][(j - 2)][i],star_state->pv[k][(j - 1)][i],star_state->pv[k][j][i],star_state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (star_state->pv[k][(j - 1)][i] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->u_lat[k][j][i]) / (sqrt((pow(star_state->u_lat[k][j][i],2) + pow(star_state->v_lat[k][j][i],2))) + 1e-24));
        star_state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,star_state->u_lat[k][j][i]),1,star_state->pv[k][j][(i - 2)],star_state->pv[k][j][(i - 1)],star_state->pv[k][j][i],star_state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (star_state->pv[k][j][(i - 1)] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlon[k][j][i] = (((star_state->mfx_lon[k][j][i] - star_state->mfx_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlat[k][j][i] = (((star_state->mfy_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->mfy_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = star_state->mfy_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mftmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlat[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      tend1->dphs[0][j][i] = 0.0;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dphs[0][j][i] = ((tend1->dphs[0][j][i] - tend1->dmfdlon[k][j][i]) - tend1->dmfdlat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        mfs = 0.0;
        for (l=0+Proc.lev_hw;l<((k - 1) + 1);l+=1){
          mfs = ((mfs + tend1->dmfdlon[l][j][i]) + tend1->dmfdlat[l][j][i]);
        }
        star_state->we_lev[k][j][i] = (-hybrid_coord_calc_dphdt_lev(mesh->hybi[k][0][0],tend1->dphs[0][j][i]) - mfs);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->we_lev[0][0][0], &Proc.FieldReq[async_we_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = advpt->we0[k][j][i];
          advpt->mm[k][j][i] = advpt->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = star_state->we_lev[k][j][i];
          advpt->mm[k][j][i] = star_state->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = ((advpt->we[k][j][i] + star_state->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->mm[k][j][i] = ((advpt->mm[k][j][i] + star_state->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->we0[k][j][i] = star_state->we_lev[k][j][i];
            advpt->m0[k][j][i] = star_state->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = (advpt->we[k][j][i] + star_state->we_lev[k][j][i]);
            advpt->mm[k][j][i] = (advpt->mm[k][j][i] + star_state->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflz[k][j][i] = ((advpt->we[k][j][i] / advpt->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lon[k][j][i] = (((mesh->area_lon_west[0][j][0] * star_state->we_lev[k][j][i]) + (mesh->area_lon_east[0][j][0] * star_state->we_lev[k][j][(i + 1)])) / mesh->area_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lon_edgewe_lev_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lat[k][j][i] = (((mesh->area_lat_north[0][j][0] * star_state->we_lev[k][(j + 1)][i]) + (mesh->area_lat_south[0][j][0] * star_state->we_lev[k][j][i])) / mesh->area_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lat_edgewe_lev_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = ((((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) + (star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i]))) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = ((((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) + (star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i]))) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[k][j][(i - 2)],star_state->pt[k][j][(i - 1)],star_state->pt[k][j][i],star_state->pt[k][j][(i + 1)],star_state->pt[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(star_state->pt[k][(j - 2)][i],star_state->pt[k][(j - 1)][i],star_state->pt[k][j][i],star_state->pt[k][(j + 1)][i],star_state->pt[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + star_state->pt[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + star_state->pt[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->vv[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->vv[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->qx[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (advpt->divx[k][j][i] * star_state->pt[k][j][i]))) * dt));
        advpt->qy[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (advpt->divy[k][j][i] * star_state->pt[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = star_state->ptf_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->qx[k][j][i] = star_state->pt[k][j][i];
        advpt->qy[k][j][i] = (star_state->pt[k][j][i] + ((0.5 * ((((star_state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]) - (advpt->divy[k][j][i] * star_state->pt[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &advpt->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(advpt->qy[k][j][(i - 2)],advpt->qy[k][j][(i - 1)],advpt->qy[k][j][i],advpt->qy[k][j][(i + 1)],advpt->qy[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(advpt->qx[k][(j - 2)][i],advpt->qx[k][(j - 1)][i],advpt->qx[k][j][i],advpt->qx[k][(j + 1)][i],advpt->qx[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + advpt->qy[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + advpt->qy[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->mfy[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->mfy[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlon[k][j][i] = (((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlat[k][j][i] = (((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = star_state->ptf_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptftmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlat[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<34+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k - 1)][j][i]) - star_state->pt[(k - 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[(k - 2)][j][i],star_state->pt[(k - 1)][j][i],star_state->pt[k][j][i],star_state->pt[(k + 1)][j][i],star_state->pt[(k + 2)][j][i]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflz[k][j][i]);
        cf = (advpt->cflz[k][j][i] - ci);
        if ((advpt->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + star_state->pt[l][j][i]);
          }
          star_state->ptf_lev[k][j][i] = ((advpt->we[k][j][i] * (((cursum + (advpt->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * advpt->dqx[((k - ci) - 1)][j][i]) * ds2)) + (advpt->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
        }
        else
        {
          if ((advpt->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + star_state->pt[l][j][i]);
            }
            star_state->ptf_lev[k][j][i] = ((-advpt->we[k][j][i] * (((cursum + (advpt->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * advpt->dqx[(k - ci)][j][i]) * ds2)) + (advpt->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
          }
          else
          {
            star_state->ptf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmptf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlev[k][j][i] = (star_state->ptf_lev[(k + 1)][j][i] - star_state->ptf_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhu[k][j][i] = (((mesh->half_tangent_wgt_0[0][j][0] * ((star_state->mfx_lon[k][j][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][(i - 1)])) + (star_state->mfx_lon[k][j][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][i])))) + (mesh->half_tangent_wgt_1[0][j][0] * ((star_state->mfx_lon[k][(j + 1)][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][(i - 1)])) + (star_state->mfx_lon[k][(j + 1)][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][i]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhuTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhv[k][j][i] = (((mesh->full_tangent_wgt_0[0][j][0] * ((star_state->mfy_lat[k][(j - 1)][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][i])) + (star_state->mfy_lat[k][(j - 1)][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][(i + 1)])))) + (mesh->full_tangent_wgt_1[0][j][0] * ((star_state->mfy_lat[k][j][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][i])) + (star_state->mfy_lat[k][j][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][(i + 1)]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlon[k][j][i] = ((star_state->ke[k][j][(i + 1)] - star_state->ke[k][j][i]) / mesh->de_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlat[k][j][i] = ((star_state->ke[k][(j + 1)][i] - star_state->ke[k][j][i]) / mesh->de_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->du[k][j][i] = ((tend1->qhv[k][j][i] - tend1->dkedlon[k][j][i]) - tend1->wedudlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dv[k][j][i] = ((-tend1->qhu[k][j][i] - tend1->dkedlat[k][j][i]) - tend1->wedvdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dpt[k][j][i] = ((-tend1->dptfdlon[k][j][i] - tend1->dptfdlat[k][j][i]) - tend1->dptfdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->phs) = 1;
  (tendPara->pt) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend1->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend1->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend1->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend1->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend1->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend1->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend1->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend1->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend1->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->t[k][j][i] = temperature(new_state->pt[k][j][i],new_state->ph[k][j][i],new_state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * new_state->t[l][j][i]) * log((new_state->ph_lev[(l + 1)][j][i] / new_state->ph_lev[l][j][i]))));
        }
        new_state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &new_state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][j][(i + 1)])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][j][(i + 1)] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][j][(i + 1)]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][j][(i + 1)]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][j][(i + 1)]);
        tend2->pgf_lon[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lon[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][(j + 1)][i])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][(j + 1)][i] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][(j + 1)][i]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][(j + 1)][i]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][(j + 1)][i]);
        tend2->pgf_lat[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lat[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_du]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->du[k][j][i] = (tend1->du[k][j][i] - tend2->pgf_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->dv[k][j][i] = (tend1->dv[k][j][i] - tend2->pgf_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->u) = 1;
  (tendPara->v) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend2->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend2->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend2->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend2->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend2->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend2->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend2->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend2->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend2->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend2->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
}

void damp_run_0(struct HybridStateField* state, struct HybridTendField* tend, struct HybridMeshField* mesh, double dt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + ((mesh->c_lon[k][j][0] * (state->div[k][j][(i + 1)] - state->div[k][j][i])) / mesh->de_lon[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + ((mesh->c_lat[k][j][0] * (state->div[k][(j + 1)][i] - state->div[k][j][i])) / mesh->de_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_t[k][j][i] = (((state->u_lon[k][j][i] - state->u_lon[k][j][(i - 1)]) / mesh->de_lon[0][j][0]) - ((((state->v_lat[k][j][i] * mesh->half_cos_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_tTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_t[0][0][0], &Proc.FieldReq[async_smag_t], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_s[k][j][i] = (((state->v_lat[k][j][(i + 1)] - state->v_lat[k][j][i]) / mesh->le_lat[0][j][0]) + ((((state->u_lon[k][(j + 1)][i] * mesh->full_cos_lat[0][(j + 1)][0]) - (state->u_lon[k][j][i] * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_sTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_s[0][0][0], &Proc.FieldReq[async_smag_s], false, false, true, true, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lon[k][j][i] = ((0.1 / ((1.0 / pow(mesh->de_lon[0][j][0],2)) + (1 / pow(mesh->le_lon[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][j][(i + 1)],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][(j - 1)][i],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lat[k][j][i] = ((0.1 / ((1.0 / pow(mesh->le_lat[0][j][0],2)) + (1 / pow(mesh->de_lat[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][(j + 1)][i],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][j][(i - 1)],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dudt[k][j][i] = (state->kmh_lon[k][j][i] * ((((state->u_lon[k][j][(i - 1)] - (2 * state->u_lon[k][j][i])) + state->u_lon[k][j][(i + 1)]) / pow(mesh->de_lon[0][j][0],2)) + ((((((state->u_lon[k][(j + 1)][i] - state->u_lon[k][j][i]) / mesh->de_lat[0][j][0]) * mesh->half_cos_lat[0][j][0]) - (((state->u_lon[k][j][i] - state->u_lon[k][(j - 1)][i]) / mesh->de_lat[0][(j - 1)][0]) * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dudtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + (dt * tend->smag_dudt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dvdt[k][j][i] = (state->kmh_lat[k][j][i] * (((state->v_lat[k][j][(i - 1)] - (2 * state->v_lat[k][j][i])) + state->v_lat[k][j][(i + 1)]) / pow(mesh->le_lat[0][j][0],2)));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dvdtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1198)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dvdt[k][j][i] = (state->kmh_lat[k][j][i] * ((((state->v_lat[k][j][(i - 1)] - (2 * state->v_lat[k][j][i])) + state->v_lat[k][j][(i + 1)]) / pow(mesh->le_lat[0][j][0],2)) + ((((((state->v_lat[k][(j + 1)][i] - state->v_lat[k][j][i]) / mesh->le_lon[0][(j + 1)][0]) * mesh->full_cos_lat[0][(j + 1)][0]) - (((state->v_lat[k][j][i] - state->v_lat[k][(j - 1)][i]) / mesh->le_lon[0][j][0]) * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dvdtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + (dt * tend->smag_dvdt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] * state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  trickyptptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,0,&state->phs[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(2,1,&state->pt[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] / state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pole_damp_runptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lon_edge(3,&state->u_lon[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lat_edge(3,&state->v_lat[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void c2a_0(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void timeadvance_0()
{
  Time_Advance();
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Time_AdvanceTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void adv_run_0(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double cursum, cf, s1, s2, ds1, ds2, ds3;
  int k, j, i, ci, l;
  double *tmp_state_old_tmpsum, *tmp_state_new_tmpsum;
  struct Vector3 buf;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = adv->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = adv->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = state_old->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = state_old->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = ((adv->uu[k][j][i] + state_old->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = ((adv->vv[k][j][i] + state_old->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->u0[k][j][i] = adv->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->v0[k][j][i] = adv->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = (adv->uu[k][j][i] + state_old->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = (adv->vv[k][j][i] + state_old->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflx[k][j][i] = ((adv->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cfly[k][j][i] = ((adv->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->divx[k][j][i] = (((adv->uu[k][j][i] - adv->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          adv->divy[k][j][i] = (((adv->vv[k][j][i] * mesh->le_lat[0][j][0]) - (adv->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state_old->tmpsum[k][j][i] = adv->vv[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:state_old->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
    tmp_state_old_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_state_old_tmpsum[k-Proc.lev_hw]=tmp_state_old_tmpsum[k-Proc.lev_hw]+state_old->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_state_old_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state_old->tmpsum[k][j][i]=tmp_state_old_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_state_old_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->divy[k][j][i] = (((state_old->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
  }
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = state_old->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = state_old->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = ((adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = ((adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = (adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = (adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = adv->we0[k][j][i];
          adv->mm[k][j][i] = adv->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = state_old->we_lev[k][j][i];
          adv->mm[k][j][i] = state_old->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = ((adv->we[k][j][i] + state_old->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->mm[k][j][i] = ((adv->mm[k][j][i] + state_old->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->we0[k][j][i] = state_old->we_lev[k][j][i];
            adv->m0[k][j][i] = state_old->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = (adv->we[k][j][i] + state_old->we_lev[k][j][i]);
            adv->mm[k][j][i] = (adv->mm[k][j][i] + state_old->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflz[k][j][i] = ((adv->we[k][j][i] / adv->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_old->qv[k][j][(i - 2)],state_old->qv[k][j][(i - 1)],state_old->qv[k][j][i],state_old->qv[k][j][(i + 1)],state_old->qv[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(state_old->qv[k][(j - 2)][i],state_old->qv[k][(j - 1)][i],state_old->qv[k][j][i],state_old->qv[k][(j + 1)][i],state_old->qv[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + state_old->qv[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + state_old->qv[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->vv[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->vv[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->qx[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (adv->divx[k][j][i] * state_old->qv[k][j][i]))) * dt));
        adv->qy[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (adv->divy[k][j][i] * state_old->qv[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_old->tmpsum[k][j][i] = adv->qmf_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state_old->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_state_old_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_old_tmpsum[k-Proc.lev_hw]=tmp_state_old_tmpsum[k-Proc.lev_hw]+state_old->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_old_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_old->tmpsum[k][j][i]=tmp_state_old_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_old_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->qx[k][j][i] = state_old->qv[k][j][i];
        adv->qy[k][j][i] = (state_old->qv[k][j][i] + ((0.5 * ((((state_old->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]) - (adv->divy[k][j][i] * state_old->qv[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &adv->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(adv->qy[k][j][(i - 2)],adv->qy[k][j][(i - 1)],adv->qy[k][j][i],adv->qy[k][j][(i + 1)],adv->qy[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(adv->qx[k][(j - 2)][i],adv->qx[k][(j - 1)][i],adv->qx[k][j][i],adv->qx[k][(j + 1)][i],adv->qx[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + adv->qy[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + adv->qy[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->mfy[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->mfy[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (((adv->old_m[k][j][i] * state_old->qv[k][j][i]) - ((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0])) + ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->tmpsum[k][j][i] = adv->qmf_lat[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state_new->tmpsum.[0, 0, 32, 0, 0, 1, 1, 0, 2400]
  tmp_state_new_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_new_tmpsum[k-Proc.lev_hw]=tmp_state_new_tmpsum[k-Proc.lev_hw]+state_new->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_new_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->tmpsum[k][j][i]=tmp_state_new_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_new_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = ((adv->old_m[k][j][i] * state_old->qv[k][j][i]) - ((((state_new->tmpsum[k][j][i] * mesh->le_lat[0][j][0]) / 2400) / mesh->area_cell[0][j][0]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_new->qv[(k - 2)][j][i],state_new->qv[(k - 1)][j][i],state_new->qv[k][j][i],state_new->qv[(k + 1)][j][i],state_new->qv[(k + 2)][j][i]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflz[k][j][i]);
        cf = (adv->cflz[k][j][i] - ci);
        if ((adv->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + state_new->qv[l][j][i]);
          }
          adv->qmf_lev[k][j][i] = ((adv->we[k][j][i] * (((cursum + (adv->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * adv->dqx[((k - ci) - 1)][j][i]) * ds2)) + (adv->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
        }
        else
        {
          if ((adv->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + state_new->qv[l][j][i]);
            }
            adv->qmf_lev[k][j][i] = ((-adv->we[k][j][i] * (((cursum + (adv->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * adv->dqx[(k - ci)][j][i]) * ds2)) + (adv->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
          }
          else
          {
            adv->qmf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqmf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = ((state_new->qv[k][j][i] * state_old->m[k][j][i]) - ((adv->qmf_lev[(k + 1)][j][i] - adv->qmf_lev[k][j][i]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,1,&state_new->qv[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->old_m[k][j][i] = state_old->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void Diagnose_State_0(struct HybridStateField* state)
{
  Diagnose(state,15);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  DiagnoseTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void ncIO_Read_1(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int nid;
  nid = ncOpenFile("gmcore.restart.DSL.15.nc",0);
  ncGetVar(nid, "u_lon", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.half_nlon, Proc.lon_hw, &(state->u_lon[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "v_lat", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.half_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->v_lat[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  
  ncGetVar(nid, "pt", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->pt[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "phs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->phs[0][0][0]));
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  
  ncGetVar(nid, "gzs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(staticv->gzs[0][0][0]));
  UpdateHalo_2d_D(Proc, &staticv->gzs[0][0][0], &Proc.FieldReq[async_gzs], true, true, true, true, true, true, false);
  
  ncCloseFile(nid);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ncCloseFileTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void spaceOperatorInit_1(struct HybridStaticField* staticv, struct HybridMeshField* mesh)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlon[0][j][i] = (((staticv->gzs[0][j][(i + 1)] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lon[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlat[0][j][i] = (((staticv->gzs[0][(j + 1)][i] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void gzlevInit_1(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->gz_lev[k][j][i] = staticv->gzs[0][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  preparegzlevgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->gz_lev[32 + Proc.lev_hw][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void uvInit_1(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void operatorPrepareNull_1(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  int k, j, i, l;
  double b, dgz;
  struct Vector4 ke_vtx;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->t[k][j][i] = temperature(state->pt[k][j][i],state->ph[k][j][i],state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lon[k][j][i] = (state->m_lon[k][j][i] * state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lat[k][j][i] = (state->m_lat[k][j][i] * state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (state->mfx_lon[k][j][(i - 1)] + state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (state->mfx_lon[k][(j + 1)][(i - 1)] + state->mfx_lon[k][(j + 1)][i])));
        state->u_lat[k][j][i] = (state->mfx_lat[k][j][i] / state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (state->mfy_lat[k][(j - 1)][i] + state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (state->mfy_lat[k][j][i] + state->mfy_lat[k][j][(i + 1)])));
        state->v_lon[k][j][i] = (state->mfy_lon[k][j][i] / state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->vor[k][j][i] = ((((state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pv[k][j][i] = ((state->vor[k][j][i] + mesh->half_f[0][j][0]) / state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->v_lon[k][j][i]) / (sqrt((pow(state->u_lon[k][j][i],2) + pow(state->v_lon[k][j][i],2))) + 1e-24));
        state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,state->v_lon[k][j][i]),1,state->pv[k][(j - 2)][i],state->pv[k][(j - 1)][i],state->pv[k][j][i],state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (state->pv[k][(j - 1)][i] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->u_lat[k][j][i]) / (sqrt((pow(state->u_lat[k][j][i],2) + pow(state->v_lat[k][j][i],2))) + 1e-24));
        state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,state->u_lat[k][j][i]),1,state->pv[k][j][(i - 2)],state->pv[k][j][(i - 1)],state->pv[k][j][i],state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (state->pv[k][j][(i - 1)] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->div[k][j][i] = ((((state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * state->t[l][j][i]) * log((state->ph_lev[(l + 1)][j][i] / state->ph_lev[l][j][i]))));
        }
        state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void advPrepare_1(struct HybridStateField* state, struct HybridAdvField* advm, struct HybridAdvField* advpt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advm->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void filterinit_1(struct HybridMeshField* mesh, double dt)
{
  Filter_Init(mesh,dt);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_InitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void stepForwardBackward_1(struct HybridStateField* old_state, struct HybridStateField* star_state, struct HybridStateField* new_state, struct HybridStaticField* staticv, struct HybridTendField* tend1, struct HybridTendField* tend2, struct HybridTendPara* tendPara, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double b, mfs, cursum, cf, s1, s2, ds1, ds2, ds3, dgz, tl, dph1, dph2, dgz1, dgz2;
  int k, j, i, l, ci;
  struct Vector4 ke_vtx;
  struct Vector3 buf;
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = star_state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = star_state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + star_state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + star_state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + star_state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + star_state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lon[k][j][i] = (star_state->m_lon[k][j][i] * star_state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lat[k][j][i] = (star_state->m_lat[k][j][i] * star_state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = star_state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = star_state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (star_state->mfx_lon[k][j][(i - 1)] + star_state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (star_state->mfx_lon[k][(j + 1)][(i - 1)] + star_state->mfx_lon[k][(j + 1)][i])));
        star_state->u_lat[k][j][i] = (star_state->mfx_lat[k][j][i] / star_state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (star_state->mfy_lat[k][(j - 1)][i] + star_state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (star_state->mfy_lat[k][j][i] + star_state->mfy_lat[k][j][(i + 1)])));
        star_state->v_lon[k][j][i] = (star_state->mfy_lon[k][j][i] / star_state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(star_state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(star_state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(star_state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * star_state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->div[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (star_state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((star_state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->vor[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (star_state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((star_state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (star_state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pv[k][j][i] = ((star_state->vor[k][j][i] + mesh->half_f[0][j][0]) / star_state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->v_lon[k][j][i]) / (sqrt((pow(star_state->u_lon[k][j][i],2) + pow(star_state->v_lon[k][j][i],2))) + 1e-24));
        star_state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,star_state->v_lon[k][j][i]),1,star_state->pv[k][(j - 2)][i],star_state->pv[k][(j - 1)][i],star_state->pv[k][j][i],star_state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (star_state->pv[k][(j - 1)][i] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->u_lat[k][j][i]) / (sqrt((pow(star_state->u_lat[k][j][i],2) + pow(star_state->v_lat[k][j][i],2))) + 1e-24));
        star_state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,star_state->u_lat[k][j][i]),1,star_state->pv[k][j][(i - 2)],star_state->pv[k][j][(i - 1)],star_state->pv[k][j][i],star_state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (star_state->pv[k][j][(i - 1)] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlon[k][j][i] = (((star_state->mfx_lon[k][j][i] - star_state->mfx_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlat[k][j][i] = (((star_state->mfy_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->mfy_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      tend1->dphs[0][j][i] = 0.0;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dphs[0][j][i] = ((tend1->dphs[0][j][i] - tend1->dmfdlon[k][j][i]) - tend1->dmfdlat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        mfs = 0.0;
        for (l=0+Proc.lev_hw;l<((k - 1) + 1);l+=1){
          mfs = ((mfs + tend1->dmfdlon[l][j][i]) + tend1->dmfdlat[l][j][i]);
        }
        star_state->we_lev[k][j][i] = (-hybrid_coord_calc_dphdt_lev(mesh->hybi[k][0][0],tend1->dphs[0][j][i]) - mfs);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->we_lev[0][0][0], &Proc.FieldReq[async_we_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = advpt->we0[k][j][i];
          advpt->mm[k][j][i] = advpt->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = star_state->we_lev[k][j][i];
          advpt->mm[k][j][i] = star_state->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = ((advpt->we[k][j][i] + star_state->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->mm[k][j][i] = ((advpt->mm[k][j][i] + star_state->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->we0[k][j][i] = star_state->we_lev[k][j][i];
            advpt->m0[k][j][i] = star_state->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = (advpt->we[k][j][i] + star_state->we_lev[k][j][i]);
            advpt->mm[k][j][i] = (advpt->mm[k][j][i] + star_state->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflz[k][j][i] = ((advpt->we[k][j][i] / advpt->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lon[k][j][i] = (((mesh->area_lon_west[0][j][0] * star_state->we_lev[k][j][i]) + (mesh->area_lon_east[0][j][0] * star_state->we_lev[k][j][(i + 1)])) / mesh->area_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lon_edgewe_lev_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lat[k][j][i] = (((mesh->area_lat_north[0][j][0] * star_state->we_lev[k][(j + 1)][i]) + (mesh->area_lat_south[0][j][0] * star_state->we_lev[k][j][i])) / mesh->area_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lat_edgewe_lev_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = ((((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) + (star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i]))) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = ((((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) + (star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i]))) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[k][j][(i - 2)],star_state->pt[k][j][(i - 1)],star_state->pt[k][j][i],star_state->pt[k][j][(i + 1)],star_state->pt[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(star_state->pt[k][(j - 2)][i],star_state->pt[k][(j - 1)][i],star_state->pt[k][j][i],star_state->pt[k][(j + 1)][i],star_state->pt[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + star_state->pt[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + star_state->pt[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->vv[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->vv[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->qx[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (advpt->divx[k][j][i] * star_state->pt[k][j][i]))) * dt));
        advpt->qy[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (advpt->divy[k][j][i] * star_state->pt[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &advpt->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(advpt->qy[k][j][(i - 2)],advpt->qy[k][j][(i - 1)],advpt->qy[k][j][i],advpt->qy[k][j][(i + 1)],advpt->qy[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(advpt->qx[k][(j - 2)][i],advpt->qx[k][(j - 1)][i],advpt->qx[k][j][i],advpt->qx[k][(j + 1)][i],advpt->qx[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + advpt->qy[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + advpt->qy[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->mfy[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->mfy[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlon[k][j][i] = (((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlat[k][j][i] = (((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<34+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k - 1)][j][i]) - star_state->pt[(k - 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[(k - 2)][j][i],star_state->pt[(k - 1)][j][i],star_state->pt[k][j][i],star_state->pt[(k + 1)][j][i],star_state->pt[(k + 2)][j][i]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflz[k][j][i]);
        cf = (advpt->cflz[k][j][i] - ci);
        if ((advpt->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + star_state->pt[l][j][i]);
          }
          star_state->ptf_lev[k][j][i] = ((advpt->we[k][j][i] * (((cursum + (advpt->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * advpt->dqx[((k - ci) - 1)][j][i]) * ds2)) + (advpt->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
        }
        else
        {
          if ((advpt->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + star_state->pt[l][j][i]);
            }
            star_state->ptf_lev[k][j][i] = ((-advpt->we[k][j][i] * (((cursum + (advpt->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * advpt->dqx[(k - ci)][j][i]) * ds2)) + (advpt->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
          }
          else
          {
            star_state->ptf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmptf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlev[k][j][i] = (star_state->ptf_lev[(k + 1)][j][i] - star_state->ptf_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhu[k][j][i] = (((mesh->half_tangent_wgt_0[0][j][0] * ((star_state->mfx_lon[k][j][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][(i - 1)])) + (star_state->mfx_lon[k][j][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][i])))) + (mesh->half_tangent_wgt_1[0][j][0] * ((star_state->mfx_lon[k][(j + 1)][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][(i - 1)])) + (star_state->mfx_lon[k][(j + 1)][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][i]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhuTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhv[k][j][i] = (((mesh->full_tangent_wgt_0[0][j][0] * ((star_state->mfy_lat[k][(j - 1)][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][i])) + (star_state->mfy_lat[k][(j - 1)][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][(i + 1)])))) + (mesh->full_tangent_wgt_1[0][j][0] * ((star_state->mfy_lat[k][j][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][i])) + (star_state->mfy_lat[k][j][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][(i + 1)]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlon[k][j][i] = ((star_state->ke[k][j][(i + 1)] - star_state->ke[k][j][i]) / mesh->de_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlat[k][j][i] = ((star_state->ke[k][(j + 1)][i] - star_state->ke[k][j][i]) / mesh->de_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->du[k][j][i] = ((tend1->qhv[k][j][i] - tend1->dkedlon[k][j][i]) - tend1->wedudlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dv[k][j][i] = ((-tend1->qhu[k][j][i] - tend1->dkedlat[k][j][i]) - tend1->wedvdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dpt[k][j][i] = ((-tend1->dptfdlon[k][j][i] - tend1->dptfdlat[k][j][i]) - tend1->dptfdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->phs) = 1;
  (tendPara->pt) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend1->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend1->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend1->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend1->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend1->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend1->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend1->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend1->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend1->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->t[k][j][i] = temperature(new_state->pt[k][j][i],new_state->ph[k][j][i],new_state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * new_state->t[l][j][i]) * log((new_state->ph_lev[(l + 1)][j][i] / new_state->ph_lev[l][j][i]))));
        }
        new_state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &new_state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][j][(i + 1)])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][j][(i + 1)] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][j][(i + 1)]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][j][(i + 1)]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][j][(i + 1)]);
        tend2->pgf_lon[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lon[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][(j + 1)][i])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][(j + 1)][i] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][(j + 1)][i]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][(j + 1)][i]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][(j + 1)][i]);
        tend2->pgf_lat[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lat[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_du]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->du[k][j][i] = (tend1->du[k][j][i] - tend2->pgf_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->dv[k][j][i] = (tend1->dv[k][j][i] - tend2->pgf_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->u) = 1;
  (tendPara->v) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend2->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend2->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend2->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend2->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend2->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend2->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend2->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend2->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend2->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend2->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
}

void damp_run_1(struct HybridStateField* state, struct HybridTendField* tend, struct HybridMeshField* mesh, double dt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + ((mesh->c_lon[k][j][0] * (state->div[k][j][(i + 1)] - state->div[k][j][i])) / mesh->de_lon[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + ((mesh->c_lat[k][j][0] * (state->div[k][(j + 1)][i] - state->div[k][j][i])) / mesh->de_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_t[k][j][i] = (((state->u_lon[k][j][i] - state->u_lon[k][j][(i - 1)]) / mesh->de_lon[0][j][0]) - ((((state->v_lat[k][j][i] * mesh->half_cos_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_tTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_t[0][0][0], &Proc.FieldReq[async_smag_t], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_s[k][j][i] = (((state->v_lat[k][j][(i + 1)] - state->v_lat[k][j][i]) / mesh->le_lat[0][j][0]) + ((((state->u_lon[k][(j + 1)][i] * mesh->full_cos_lat[0][(j + 1)][0]) - (state->u_lon[k][j][i] * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_sTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_s[0][0][0], &Proc.FieldReq[async_smag_s], false, false, true, true, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lon[k][j][i] = ((0.1 / ((1.0 / pow(mesh->de_lon[0][j][0],2)) + (1 / pow(mesh->le_lon[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][j][(i + 1)],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][(j - 1)][i],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lat[k][j][i] = ((0.1 / ((1.0 / pow(mesh->le_lat[0][j][0],2)) + (1 / pow(mesh->de_lat[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][(j + 1)][i],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][j][(i - 1)],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dudt[k][j][i] = (state->kmh_lon[k][j][i] * ((((state->u_lon[k][j][(i - 1)] - (2 * state->u_lon[k][j][i])) + state->u_lon[k][j][(i + 1)]) / pow(mesh->de_lon[0][j][0],2)) + ((((((state->u_lon[k][(j + 1)][i] - state->u_lon[k][j][i]) / mesh->de_lat[0][j][0]) * mesh->half_cos_lat[0][j][0]) - (((state->u_lon[k][j][i] - state->u_lon[k][(j - 1)][i]) / mesh->de_lat[0][(j - 1)][0]) * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dudtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + (dt * tend->smag_dudt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1198)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dvdt[k][j][i] = (state->kmh_lat[k][j][i] * ((((state->v_lat[k][j][(i - 1)] - (2 * state->v_lat[k][j][i])) + state->v_lat[k][j][(i + 1)]) / pow(mesh->le_lat[0][j][0],2)) + ((((((state->v_lat[k][(j + 1)][i] - state->v_lat[k][j][i]) / mesh->le_lon[0][(j + 1)][0]) * mesh->full_cos_lat[0][(j + 1)][0]) - (((state->v_lat[k][j][i] - state->v_lat[k][(j - 1)][i]) / mesh->le_lon[0][j][0]) * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dvdtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + (dt * tend->smag_dvdt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] * state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  trickyptptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,0,&state->phs[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(2,1,&state->pt[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] / state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pole_damp_runptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lon_edge(3,&state->u_lon[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lat_edge(3,&state->v_lat[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void c2a_1(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void timeadvance_1()
{
  Time_Advance();
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Time_AdvanceTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void adv_run_1(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double cursum, cf, s1, s2, ds1, ds2, ds3;
  int k, j, i, ci, l;
  struct Vector3 buf;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = adv->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = adv->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = state_old->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = state_old->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = ((adv->uu[k][j][i] + state_old->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = ((adv->vv[k][j][i] + state_old->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->u0[k][j][i] = adv->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->v0[k][j][i] = adv->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = (adv->uu[k][j][i] + state_old->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = (adv->vv[k][j][i] + state_old->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflx[k][j][i] = ((adv->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cfly[k][j][i] = ((adv->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->divx[k][j][i] = (((adv->uu[k][j][i] - adv->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          adv->divy[k][j][i] = (((adv->vv[k][j][i] * mesh->le_lat[0][j][0]) - (adv->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    
  }
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = state_old->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = state_old->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = ((adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = ((adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = (adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = (adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = adv->we0[k][j][i];
          adv->mm[k][j][i] = adv->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = state_old->we_lev[k][j][i];
          adv->mm[k][j][i] = state_old->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = ((adv->we[k][j][i] + state_old->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->mm[k][j][i] = ((adv->mm[k][j][i] + state_old->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->we0[k][j][i] = state_old->we_lev[k][j][i];
            adv->m0[k][j][i] = state_old->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = (adv->we[k][j][i] + state_old->we_lev[k][j][i]);
            adv->mm[k][j][i] = (adv->mm[k][j][i] + state_old->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflz[k][j][i] = ((adv->we[k][j][i] / adv->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_old->qv[k][j][(i - 2)],state_old->qv[k][j][(i - 1)],state_old->qv[k][j][i],state_old->qv[k][j][(i + 1)],state_old->qv[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(state_old->qv[k][(j - 2)][i],state_old->qv[k][(j - 1)][i],state_old->qv[k][j][i],state_old->qv[k][(j + 1)][i],state_old->qv[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + state_old->qv[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + state_old->qv[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->vv[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->vv[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->qx[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (adv->divx[k][j][i] * state_old->qv[k][j][i]))) * dt));
        adv->qy[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (adv->divy[k][j][i] * state_old->qv[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &adv->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(adv->qy[k][j][(i - 2)],adv->qy[k][j][(i - 1)],adv->qy[k][j][i],adv->qy[k][j][(i + 1)],adv->qy[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(adv->qx[k][(j - 2)][i],adv->qx[k][(j - 1)][i],adv->qx[k][j][i],adv->qx[k][(j + 1)][i],adv->qx[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + adv->qy[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + adv->qy[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->mfy[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->mfy[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (((adv->old_m[k][j][i] * state_old->qv[k][j][i]) - ((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0])) + ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_new->qv[(k - 2)][j][i],state_new->qv[(k - 1)][j][i],state_new->qv[k][j][i],state_new->qv[(k + 1)][j][i],state_new->qv[(k + 2)][j][i]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflz[k][j][i]);
        cf = (adv->cflz[k][j][i] - ci);
        if ((adv->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + state_new->qv[l][j][i]);
          }
          adv->qmf_lev[k][j][i] = ((adv->we[k][j][i] * (((cursum + (adv->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * adv->dqx[((k - ci) - 1)][j][i]) * ds2)) + (adv->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
        }
        else
        {
          if ((adv->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + state_new->qv[l][j][i]);
            }
            adv->qmf_lev[k][j][i] = ((-adv->we[k][j][i] * (((cursum + (adv->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * adv->dqx[(k - ci)][j][i]) * ds2)) + (adv->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
          }
          else
          {
            adv->qmf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqmf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = ((state_new->qv[k][j][i] * state_old->m[k][j][i]) - ((adv->qmf_lev[(k + 1)][j][i] - adv->qmf_lev[k][j][i]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,1,&state_new->qv[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->old_m[k][j][i] = state_old->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void Diagnose_State_1(struct HybridStateField* state)
{
  Diagnose(state,15);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  DiagnoseTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void ncIO_Read_2(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int nid;
  nid = ncOpenFile("gmcore.restart.DSL.15.nc",0);
  ncGetVar(nid, "u_lon", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.half_nlon, Proc.lon_hw, &(state->u_lon[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "v_lat", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.half_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->v_lat[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  
  ncGetVar(nid, "pt", 0, 32, Proc.lev_hw, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->pt[0][0][0]));
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  
  ncGetVar(nid, "phs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(state->phs[0][0][0]));
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  
  ncGetVar(nid, "gzs", 0, 1, 0, Proc.lat_beg, Proc.full_nlat, Proc.lat_hw, Proc.lon_beg, Proc.full_nlon, Proc.lon_hw, &(staticv->gzs[0][0][0]));
  UpdateHalo_2d_D(Proc, &staticv->gzs[0][0][0], &Proc.FieldReq[async_gzs], true, true, true, true, true, true, false);
  
  ncCloseFile(nid);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ncCloseFileTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void spaceOperatorInit_2(struct HybridStaticField* staticv, struct HybridMeshField* mesh)
{
  int j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlon[0][j][i] = (((staticv->gzs[0][j][(i + 1)] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lon[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      staticv->dzsdlat[0][j][i] = (((staticv->gzs[0][(j + 1)][i] - staticv->gzs[0][j][i]) / 9.80616) / mesh->de_lat[0][j][0]);
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  prepareStaticdzsdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void gzlevInit_2(struct HybridStateField* state, struct HybridStaticField* staticv)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->gz_lev[k][j][i] = staticv->gzs[0][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  preparegzlevgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->gz_lev[32 + Proc.lev_hw][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void uvInit_2(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void operatorPrepareNull_2(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  int k, j, i, l;
  double b, dgz;
  double *tmp_state_tmpsum;
  struct Vector4 ke_vtx;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->t[k][j][i] = temperature(state->pt[k][j][i],state->ph[k][j][i],state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state->tmpsum[k][j][i] = -advpt->vv[k][(j - 1)][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
    tmp_state_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_state_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divy[k][j][i] = (((state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lon[k][j][i] = (state->m_lon[k][j][i] * state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lat[k][j][i] = (state->m_lat[k][j][i] * state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (state->mfx_lon[k][j][(i - 1)] + state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (state->mfx_lon[k][(j + 1)][(i - 1)] + state->mfx_lon[k][(j + 1)][i])));
        state->u_lat[k][j][i] = (state->mfx_lat[k][j][i] / state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (state->mfy_lat[k][(j - 1)][i] + state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (state->mfy_lat[k][j][i] + state->mfy_lat[k][j][(i + 1)])));
        state->v_lon[k][j][i] = (state->mfy_lon[k][j][i] / state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = pow(state->v_lat[k][(j - 1)][i],2);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ketmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ke[k][j][i] = (state->tmpsum[k][j][i] / 2400);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->vor[k][j][i] = ((((state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = (state->u_lat[k][j][i] * mesh->le_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vortmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 1198, 1199, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->vor[k][j][i] = ((state->tmpsum[k][j][i] / 2400) / mesh->area_cell[0][(j + 1)][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pv[k][j][i] = ((state->vor[k][j][i] + mesh->half_f[0][j][0]) / state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->v_lon[k][j][i]) / (sqrt((pow(state->u_lon[k][j][i],2) + pow(state->v_lon[k][j][i],2))) + 1e-24));
        state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,state->v_lon[k][j][i]),1,state->pv[k][(j - 2)][i],state->pv[k][(j - 1)][i],state->pv[k][j][i],state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (state->pv[k][(j - 1)][i] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(state->u_lat[k][j][i]) / (sqrt((pow(state->u_lat[k][j][i],2) + pow(state->v_lat[k][j][i],2))) + 1e-24));
        state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,state->u_lat[k][j][i]),1,state->pv[k][j][(i - 2)],state->pv[k][j][(i - 1)],state->pv[k][j][i],state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (state->pv[k][j][(i - 1)] + state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->div[k][j][i] = ((((state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i] = -state->v_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_tmpsum[k-Proc.lev_hw]=tmp_state_tmpsum[k-Proc.lev_hw]+state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->tmpsum[k][j][i]=tmp_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->div[k][j][i] = (((state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * state->t[l][j][i]) * log((state->ph_lev[(l + 1)][j][i] / state->ph_lev[l][j][i]))));
        }
        state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void advPrepare_2(struct HybridStateField* state, struct HybridAdvField* advm, struct HybridAdvField* advpt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advm->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->old_m[k][j][i] = state->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void filterinit_2(struct HybridMeshField* mesh, double dt)
{
  Filter_Init(mesh,dt);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_InitTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void stepForwardBackward_2(struct HybridStateField* old_state, struct HybridStateField* star_state, struct HybridStateField* new_state, struct HybridStaticField* staticv, struct HybridTendField* tend1, struct HybridTendField* tend2, struct HybridTendPara* tendPara, struct HybridAdvField* advpt, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double b, mfs, cursum, cf, s1, s2, ds1, ds2, ds3, dgz, tl, dph1, dph2, dgz1, dgz2;
  int k, j, i, l, ci;
  double *tmp_star_state_tmpsum;
  struct Vector4 ke_vtx;
  struct Vector3 buf;
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = advpt->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = advpt->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->uu[k][j][i] = star_state->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->vv[k][j][i] = star_state->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = ((advpt->uu[k][j][i] + star_state->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = ((advpt->vv[k][j][i] + star_state->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->u0[k][j][i] = advpt->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->v0[k][j][i] = advpt->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->uu[k][j][i] = (advpt->uu[k][j][i] + star_state->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->vv[k][j][i] = (advpt->vv[k][j][i] + star_state->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &advpt->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflx[k][j][i] = ((advpt->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cfly[k][j][i] = ((advpt->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divx[k][j][i] = (((advpt->uu[k][j][i] - advpt->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          advpt->divy[k][j][i] = (((advpt->vv[k][j][i] * mesh->le_lat[0][j][0]) - (advpt->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          star_state->tmpsum[k][j][i] = -advpt->vv[k][(j - 1)][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
    tmp_star_state_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_star_state_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->divy[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lon[k][j][i] = (star_state->m_lon[k][j][i] * star_state->u_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfx_lon[0][0][0], &Proc.FieldReq[async_mfx_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lat[k][j][i] = (star_state->m_lat[k][j][i] * star_state->v_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->mfy_lat[0][0][0], &Proc.FieldReq[async_mfy_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfx[k][j][i] = star_state->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->mfy[k][j][i] = star_state->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = ((advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = ((advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfx[k][j][i] = (advpt->mfx[k][j][i] + star_state->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->mfy[k][j][i] = (advpt->mfy[k][j][i] + star_state->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfx_lat[k][j][i] = ((mesh->half_tangent_wgt_0[0][j][0] * (star_state->mfx_lon[k][j][(i - 1)] + star_state->mfx_lon[k][j][i])) + (mesh->half_tangent_wgt_1[0][j][0] * (star_state->mfx_lon[k][(j + 1)][(i - 1)] + star_state->mfx_lon[k][(j + 1)][i])));
        star_state->u_lat[k][j][i] = (star_state->mfx_lat[k][j][i] / star_state->m_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfx_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->mfy_lon[k][j][i] = ((mesh->full_tangent_wgt_0[0][j][0] * (star_state->mfy_lat[k][(j - 1)][i] + star_state->mfy_lat[k][(j - 1)][(i + 1)])) + (mesh->full_tangent_wgt_1[0][j][0] * (star_state->mfy_lat[k][j][i] + star_state->mfy_lat[k][j][(i + 1)])));
        star_state->v_lon[k][j][i] = (star_state->mfy_lon[k][j][i] / star_state->m_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mfmfy_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  ke_vtx = (struct Vector4){0.0,0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = (((((mesh->area_lon_west[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2)) + (mesh->area_lon_east[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lat_north[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lat_south[0][j][0] * pow(star_state->v_lat[k][j][i],2))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = ((((1.0 - 0.5) * ((((((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][(i - 1)],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][i],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) + (mesh->area_lon_south[0][(j + 1)][0] * pow(star_state->u_lon[k][(j + 1)][(i - 1)],2))) / mesh->area_vtx[0][j][0]) + (((((mesh->area_lat_east[0][j][0] * pow(star_state->v_lat[k][j][i],2)) + (mesh->area_lat_west[0][j][0] * pow(star_state->v_lat[k][j][(i + 1)],2))) + (mesh->area_lon_north[0][j][0] * pow(star_state->u_lon[k][j][i],2))) + (mesh->area_lon_south[0][0][(j + 1)] * pow(star_state->u_lon[k][(j + 1)][i],2))) / mesh->area_vtx[0][j][0])) * mesh->area_subcell_1[0][j][0]) + (((((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i - 1)],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][(i - 1)],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][(i - 1)],2))) / mesh->area_vtx[0][(j - 1)][0]) + (((((mesh->area_lat_east[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][i],2)) + (mesh->area_lat_west[0][(j - 1)][0] * pow(star_state->v_lat[k][(j - 1)][(i + 1)],2))) + (mesh->area_lon_north[0][(j - 1)][0] * pow(star_state->u_lon[k][(j - 1)][i],2))) + (mesh->area_lon_south[0][j][0] * pow(star_state->u_lon[k][j][i],2))) / mesh->area_vtx[0][(j - 1)][0])) * mesh->area_subcell_0[0][j][0]))) / mesh->area_cell[0][j][0]) + (0.5 * star_state->ke[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = pow(star_state->v_lat[k][(j - 1)][i],2);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ketmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->ke[k][j][i] = (star_state->tmpsum[k][j][i] / 2400);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_kekeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ke[0][0][0], &Proc.FieldReq[async_ke], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->div[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->le_lon[0][j][0]) - (star_state->u_lon[k][j][(i - 1)] * mesh->le_lon[0][j][0])) + ((star_state->v_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->v_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0]))) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = -star_state->v_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->div[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_divdivTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->div[0][0][0], &Proc.FieldReq[async_div], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->vor[k][j][i] = ((((star_state->u_lon[k][j][i] * mesh->de_lon[0][j][0]) - (star_state->u_lon[k][(j + 1)][i] * mesh->de_lon[0][(j + 1)][0])) + ((star_state->v_lat[k][j][(i + 1)] * mesh->de_lat[0][j][0]) - (star_state->v_lat[k][j][i] * mesh->de_lat[0][j][0]))) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = (star_state->u_lat[k][j][i] * mesh->le_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vortmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1198, 1199, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->vor[k][j][i] = ((star_state->tmpsum[k][j][i] / 2400) / mesh->area_cell[0][(j + 1)][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_vorvorTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pv[k][j][i] = ((star_state->vor[k][j][i] + mesh->half_f[0][j][0]) / star_state->m_vtx[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_pvpvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv[0][0][0], &Proc.FieldReq[async_pv], false, false, true, true, true, true, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->v_lon[k][j][i]) / (sqrt((pow(star_state->u_lon[k][j][i],2) + pow(star_state->v_lon[k][j][i],2))) + 1e-24));
        star_state->pv_lon[k][j][i] = ((b * upwind3(sign(1.0,star_state->v_lon[k][j][i]),1,star_state->pv[k][(j - 2)][i],star_state->pv[k][(j - 1)][i],star_state->pv[k][j][i],star_state->pv[k][(j + 1)][i])) + (((1 - b) * 0.5) * (star_state->pv[k][(j - 1)][i] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lon[0][0][0], &Proc.FieldReq[async_pv_lon], false, true, true, true, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        b = (fabs(star_state->u_lat[k][j][i]) / (sqrt((pow(star_state->u_lat[k][j][i],2) + pow(star_state->v_lat[k][j][i],2))) + 1e-24));
        star_state->pv_lat[k][j][i] = ((b * upwind3(sign(1.0,star_state->u_lat[k][j][i]),1,star_state->pv[k][j][(i - 2)],star_state->pv[k][j][(i - 1)],star_state->pv[k][j][i],star_state->pv[k][j][(i + 1)])) + (((1 - b) * 0.5) * (star_state->pv[k][j][(i - 1)] + star_state->pv[k][j][i])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_pv_upwindpv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pv_lat[0][0][0], &Proc.FieldReq[async_pv_lat], true, false, true, false, true, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlon[k][j][i] = (((star_state->mfx_lon[k][j][i] - star_state->mfx_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlat[k][j][i] = (((star_state->mfy_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->mfy_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = -star_state->mfy_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mftmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dmfdlat[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_mfdmfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
    for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
      tend1->dphs[0][j][i] = 0.0;
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dphs[0][j][i] = ((tend1->dphs[0][j][i] - tend1->dmfdlon[k][j][i]) - tend1->dmfdlat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_dphsdtdphsTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev[k][j][i] = 0.0;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dphs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        mfs = 0.0;
        for (l=0+Proc.lev_hw;l<((k - 1) + 1);l+=1){
          mfs = ((mfs + tend1->dmfdlon[l][j][i]) + tend1->dmfdlat[l][j][i]);
        }
        star_state->we_lev[k][j][i] = (-hybrid_coord_calc_dphdt_lev(mesh->hybi[k][0][0],tend1->dphs[0][j][i]) - mfs);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_we_levwe_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->we_lev[0][0][0], &Proc.FieldReq[async_we_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = advpt->we0[k][j][i];
          advpt->mm[k][j][i] = advpt->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->we[k][j][i] = star_state->we_lev[k][j][i];
          advpt->mm[k][j][i] = star_state->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = ((advpt->we[k][j][i] + star_state->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->mm[k][j][i] = ((advpt->mm[k][j][i] + star_state->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            advpt->we0[k][j][i] = star_state->we_lev[k][j][i];
            advpt->m0[k][j][i] = star_state->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            advpt->we[k][j][i] = (advpt->we[k][j][i] + star_state->we_lev[k][j][i]);
            advpt->mm[k][j][i] = (advpt->mm[k][j][i] + star_state->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          advpt->cflz[k][j][i] = ((advpt->we[k][j][i] / advpt->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lon[k][j][i] = (((mesh->area_lon_west[0][j][0] * star_state->we_lev[k][j][i]) + (mesh->area_lon_east[0][j][0] * star_state->we_lev[k][j][(i + 1)])) / mesh->area_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lon_edgewe_lev_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->we_lev_lat[k][j][i] = (((mesh->area_lat_north[0][j][0] * star_state->we_lev[k][(j + 1)][i]) + (mesh->area_lat_south[0][j][0] * star_state->we_lev[k][j][i])) / mesh->area_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interp_lev_edge_to_lev_lat_edgewe_lev_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = ((((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) + (star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i]))) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[(k + 1)][j][i] * (star_state->u_lon[(k + 1)][j][i] - star_state->u_lon[k][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedudlev[k][j][i] = (((star_state->we_lev_lon[k][j][i] * (star_state->u_lon[k][j][i] - star_state->u_lon[(k - 1)][j][i])) / star_state->m_lon[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedudlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<31+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = ((((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) + (star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i]))) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[(k + 1)][j][i] * (star_state->v_lat[(k + 1)][j][i] - star_state->v_lat[k][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=31+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->wedvdlev[k][j][i] = (((star_state->we_lev_lat[k][j][i] * (star_state->v_lat[k][j][i] - star_state->v_lat[(k - 1)][j][i])) / star_state->m_lat[k][j][i]) / 2.0);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_wedudlev_wedvdlevwedvdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[k][j][(i - 2)],star_state->pt[k][j][(i - 1)],star_state->pt[k][j][i],star_state->pt[k][j][(i + 1)],star_state->pt[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(star_state->pt[k][(j - 2)][i],star_state->pt[k][(j - 1)][i],star_state->pt[k][j][i],star_state->pt[k][(j + 1)][i],star_state->pt[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + star_state->pt[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + star_state->pt[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->uu[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->vv[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->vv[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->qx[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (advpt->divx[k][j][i] * star_state->pt[k][j][i]))) * dt));
        advpt->qy[k][j][i] = (star_state->pt[k][j][i] - ((0.5 * ((((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (advpt->divy[k][j][i] * star_state->pt[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = star_state->ptf_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        advpt->qx[k][j][i] = star_state->pt[k][j][i];
        advpt->qy[k][j][i] = (star_state->pt[k][j][i] + ((0.5 * ((((star_state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]) - (advpt->divy[k][j][i] * star_state->pt[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &advpt->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(advpt->qy[k][j][(i - 2)],advpt->qy[k][j][(i - 1)],advpt->qy[k][j][i],advpt->qy[k][j][(i + 1)],advpt->qy[k][j][(i + 2)]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
        buf = ppm(advpt->qx[k][(j - 2)][i],advpt->qx[k][(j - 1)][i],advpt->qx[k][j][i],advpt->qx[k][(j + 1)][i],advpt->qx[k][(j + 2)][i]);
        advpt->qly[k][j][i] = buf.x;
        advpt->dqy[k][j][i] = buf.y;
        advpt->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &advpt->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflx[k][j][i]);
        cf = (advpt->cflx[k][j][i] - ci);
        if ((advpt->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + advpt->qy[k][j][l]);
          }
          star_state->ptf_lon[k][j][i] = ((advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * advpt->dqx[k][j][(i - ci)]) * ds2)) + (advpt->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
        }
        else
        {
          if ((advpt->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + advpt->qy[k][j][l]);
            }
            star_state->ptf_lon[k][j][i] = ((-advpt->mfx[k][j][i] * (((cursum + (advpt->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * advpt->dqx[k][j][((i - ci) + 1)]) * ds2)) + (advpt->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflx[k][j][i]);
          }
          else
          {
            star_state->ptf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lon[0][0][0], &Proc.FieldReq[async_ptf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((advpt->cfly[k][j][i] > 0))
        {
          s1 = (1 - advpt->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          star_state->ptf_lat[k][j][i] = ((advpt->mfy[k][j][i] * (((advpt->qly[k][j][i] * ds1) + ((0.5 * advpt->dqy[k][j][i]) * ds2)) + (advpt->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
        }
        else
        {
          if ((advpt->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -advpt->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            star_state->ptf_lat[k][j][i] = ((-advpt->mfy[k][j][i] * (((advpt->qly[k][(j + 1)][i] * ds1) + ((0.5 * advpt->dqy[k][(j + 1)][i]) * ds2)) + (advpt->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cfly[k][j][i]);
          }
          else
          {
            star_state->ptf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerptf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->ptf_lat[0][0][0], &Proc.FieldReq[async_ptf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlon[k][j][i] = (((star_state->ptf_lon[k][j][i] - star_state->ptf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlat[k][j][i] = (((star_state->ptf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (star_state->ptf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ptf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i] = -star_state->ptf_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptftmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:star_state->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_star_state_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_star_state_tmpsum[k-Proc.lev_hw]=tmp_star_state_tmpsum[k-Proc.lev_hw]+star_state->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_star_state_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->tmpsum[k][j][i]=tmp_star_state_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_star_state_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlat[k][j][i] = (((star_state->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-1+Proc.lev_hw; k<0+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=-2+Proc.lev_hw; k<-1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k + 1)][j][i]) - star_state->pt[(k + 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<34+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        star_state->pt[k][j][i] = ((2 * star_state->pt[(k - 1)][j][i]) - star_state->pt[(k - 2)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &star_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(star_state->pt[(k - 2)][j][i],star_state->pt[(k - 1)][j][i],star_state->pt[k][j][i],star_state->pt[(k + 1)][j][i],star_state->pt[(k + 2)][j][i]);
        advpt->qlx[k][j][i] = buf.x;
        advpt->dqx[k][j][i] = buf.y;
        advpt->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &advpt->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &advpt->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(advpt->cflz[k][j][i]);
        cf = (advpt->cflz[k][j][i] - ci);
        if ((advpt->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + star_state->pt[l][j][i]);
          }
          star_state->ptf_lev[k][j][i] = ((advpt->we[k][j][i] * (((cursum + (advpt->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * advpt->dqx[((k - ci) - 1)][j][i]) * ds2)) + (advpt->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
        }
        else
        {
          if ((advpt->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + star_state->pt[l][j][i]);
            }
            star_state->ptf_lev[k][j][i] = ((-advpt->we[k][j][i] * (((cursum + (advpt->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * advpt->dqx[(k - ci)][j][i]) * ds2)) + (advpt->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / advpt->cflz[k][j][i]);
          }
          else
          {
            star_state->ptf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmptf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dptfdlev[k][j][i] = (star_state->ptf_lev[(k + 1)][j][i] - star_state->ptf_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_ptfdptfdlevTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhu[k][j][i] = (((mesh->half_tangent_wgt_0[0][j][0] * ((star_state->mfx_lon[k][j][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][(i - 1)])) + (star_state->mfx_lon[k][j][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][j][i])))) + (mesh->half_tangent_wgt_1[0][j][0] * ((star_state->mfx_lon[k][(j + 1)][(i - 1)] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][(i - 1)])) + (star_state->mfx_lon[k][(j + 1)][i] * (star_state->pv_lat[k][j][i] + star_state->pv_lon[k][(j + 1)][i]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhuTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_pv_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->qhv[k][j][i] = (((mesh->full_tangent_wgt_0[0][j][0] * ((star_state->mfy_lat[k][(j - 1)][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][i])) + (star_state->mfy_lat[k][(j - 1)][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][(j - 1)][(i + 1)])))) + (mesh->full_tangent_wgt_1[0][j][0] * ((star_state->mfy_lat[k][j][i] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][i])) + (star_state->mfy_lat[k][j][(i + 1)] * (star_state->pv_lon[k][j][i] + star_state->pv_lat[k][j][(i + 1)]))))) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_coriolisqhvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlon[k][j][i] = ((star_state->ke[k][j][(i + 1)] - star_state->ke[k][j][i]) / mesh->de_lon[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_ke]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dkedlat[k][j][i] = ((star_state->ke[k][(j + 1)][i] - star_state->ke[k][j][i]) / mesh->de_lat[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_grad_kedkedlatTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->du[k][j][i] = ((tend1->qhv[k][j][i] - tend1->dkedlon[k][j][i]) - tend1->wedudlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dv[k][j][i] = ((-tend1->qhu[k][j][i] - tend1->dkedlat[k][j][i]) - tend1->wedvdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend1->dpt[k][j][i] = ((-tend1->dptfdlon[k][j][i] - tend1->dptfdlat[k][j][i]) - tend1->dptfdlev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_forwarddptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->phs) = 1;
  (tendPara->pt) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend1->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend1->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend1->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend1->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend1->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend1->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend1->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend1->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend1->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend1->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend1->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  tendPara = tendPara;
  advPara = advPara;
  dt = dt;
  tendPara = tendPara;
  (tendPara->u) = 0;
  (tendPara->v) = 0;
  (tendPara->pt) = 0;
  (tendPara->gz) = 0;
  (tendPara->phs) = 0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->t[k][j][i] = temperature(new_state->pt[k][j][i],new_state->ph[k][j][i],new_state->qv[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_ttTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_gzs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        dgz = 0.0;
        for (l=k;l<((32 - 1) + 1)+Proc.lev_hw;l+=1){
          dgz = (dgz + ((287.04 * new_state->t[l][j][i]) * log((new_state->ph_lev[(l + 1)][j][i] / new_state->ph_lev[l][j][i]))));
        }
        new_state->gz_lev[k][j][i] = (staticv->gzs[0][j][i] + dgz);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_gz_levgz_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &new_state->gz_lev[0][0][0], &Proc.FieldReq[async_gz_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][j][(i + 1)])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][j][(i + 1)] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][j][(i + 1)]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][j][(i + 1)]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][j][(i + 1)]);
        tend2->pgf_lon[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lon[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qm]);
  HaloWait(Proc,&Proc.FieldReq[async_ph_exn_lev]);
  HaloWait(Proc,&Proc.FieldReq[async_gz_lev]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tl = (1 + (0.5 * (new_state->qm[k][j][i] + new_state->qm[k][(j + 1)][i])));
        dph1 = (new_state->ph_exn_lev[(k + 1)][(j + 1)][i] - new_state->ph_exn_lev[k][j][i]);
        dph2 = (new_state->ph_exn_lev[(k + 1)][j][i] - new_state->ph_exn_lev[k][(j + 1)][i]);
        dgz1 = (new_state->gz_lev[(k + 1)][j][i] - new_state->gz_lev[k][(j + 1)][i]);
        dgz2 = (new_state->gz_lev[k][j][i] - new_state->gz_lev[(k + 1)][(j + 1)][i]);
        tend2->pgf_lat[k][j][i] = (((-((dph1 * dgz1) + (dph2 * dgz2)) / mesh->de_lat[0][j][0]) / (dph1 + dph2)) / tl);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pgf_lin97pgf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  tendPara = tendPara;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_du]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->du[k][j][i] = (tend1->du[k][j][i] - tend2->pgf_lon[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwardduTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_dv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend2->dv[k][j][i] = (tend1->dv[k][j][i] - tend2->pgf_lat[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_tend_backwarddvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  (tendPara->u) = 1;
  (tendPara->v) = 1;
  tendPara = tendPara;
  dt = dt;
  if ((tendPara->phs))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,0,&tend2->dphs[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &tend2->dphs[0][0][0], &Proc.FieldReq[async_dphs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    HaloWait(Proc,&Proc.FieldReq[async_dphs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        new_state->phs[0][j][i] = (old_state->phs[0][j][i] + (dt * tend2->dphs[0][j][i]));
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statephsTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_2d_D(Proc, &new_state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_phs]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],new_state->phs[0][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph_exn_lev[k][j][i] = pow(new_state->ph_lev[k][j][i],0.2858964143426295);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phph_exn_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->ph[k][j][i] = (0.5 * (new_state->ph_lev[k][j][i] + new_state->ph_lev[(k + 1)][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_phphTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->pt))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_On_Cell(0,1,&tend2->dpt[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_On_CellTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dpt[0][0][0], &Proc.FieldReq[async_dpt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_pt]);
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    HaloWait(Proc,&Proc.FieldReq[async_dpt]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->pt[k][j][i] = (((old_state->pt[k][j][i] * old_state->m[k][j][i]) + (dt * tend2->dpt[k][j][i])) / new_state->m[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateptTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if ((tendPara->gz))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->gz[k][j][i] = (old_state->gz[k][j][i] + (dt * tend2->dgz[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stategzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m[k][j][i] = (new_state->ph_lev[(k + 1)][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mmTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph[k][j][i] - new_state->ph_lev[k][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lev[k][j][i] = (new_state->ph_lev[k][j][i] - new_state->ph[(k - 1)][j][i]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    calc_mm_levTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lon[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_lat[k][j][i] = ((new_state->m[k][j][i] + new_state->m[k][(j + 1)][i]) * 0.5);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_m]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->m_vtx[k][j][i] = ((((new_state->m[k][j][i] + new_state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((new_state->m[k][(j + 1)][i] + new_state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  if (((tendPara->u) && (tendPara->v)))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lon_edge(0,&tend2->du[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->du[0][0][0], &Proc.FieldReq[async_du], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    Filter_on_lat_edge(0,&tend2->dv[0][0][0]);
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &tend2->dv[0][0][0], &Proc.FieldReq[async_dv], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    HaloWait(Proc,&Proc.FieldReq[async_du]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->u_lon[k][j][i] = (old_state->u_lon[k][j][i] + (dt * tend2->du[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_stateu_lonTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    HaloWait(Proc,&Proc.FieldReq[async_dv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          new_state->v_lat[k][j][i] = (old_state->v_lat[k][j][i] + (dt * tend2->dv[k][j][i]));
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    update_statev_latTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &new_state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
}

void damp_run_2(struct HybridStateField* state, struct HybridTendField* tend, struct HybridMeshField* mesh, double dt)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + ((mesh->c_lon[k][j][0] * (state->div[k][j][(i + 1)] - state->div[k][j][i])) / mesh->de_lon[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_div]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + ((mesh->c_lat[k][j][0] * (state->div[k][(j + 1)][i] - state->div[k][j][i])) / mesh->de_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  div_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  dt = dt;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_t[k][j][i] = (((state->u_lon[k][j][i] - state->u_lon[k][j][(i - 1)]) / mesh->de_lon[0][j][0]) - ((((state->v_lat[k][j][i] * mesh->half_cos_lat[0][j][0]) - (state->v_lat[k][(j - 1)][i] * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_tTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_t[0][0][0], &Proc.FieldReq[async_smag_t], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->smag_s[k][j][i] = (((state->v_lat[k][j][(i + 1)] - state->v_lat[k][j][i]) / mesh->le_lat[0][j][0]) + ((((state->u_lon[k][(j + 1)][i] * mesh->full_cos_lat[0][(j + 1)][0]) - (state->u_lon[k][j][i] * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_sTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->smag_s[0][0][0], &Proc.FieldReq[async_smag_s], false, false, true, true, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lon[k][j][i] = ((0.1 / ((1.0 / pow(mesh->de_lon[0][j][0],2)) + (1 / pow(mesh->le_lon[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][j][(i + 1)],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][(j - 1)][i],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_smag_t]);
  HaloWait(Proc,&Proc.FieldReq[async_smag_s]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->kmh_lat[k][j][i] = ((0.1 / ((1.0 / pow(mesh->le_lat[0][j][0],2)) + (1 / pow(mesh->de_lat[0][j][0],2)))) * sqrt(((0.5 * (pow(state->smag_t[k][j][i],2) + pow(state->smag_t[k][(j + 1)][i],2))) + (0.5 * (pow(state->smag_s[k][j][i],2) + pow(state->smag_s[k][j][(i - 1)],2))))));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runkmh_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dudt[k][j][i] = (state->kmh_lon[k][j][i] * ((((state->u_lon[k][j][(i - 1)] - (2 * state->u_lon[k][j][i])) + state->u_lon[k][j][(i + 1)]) / pow(mesh->de_lon[0][j][0],2)) + ((((((state->u_lon[k][(j + 1)][i] - state->u_lon[k][j][i]) / mesh->de_lat[0][j][0]) * mesh->half_cos_lat[0][j][0]) - (((state->u_lon[k][j][i] - state->u_lon[k][(j - 1)][i]) / mesh->de_lat[0][(j - 1)][0]) * mesh->half_cos_lat[0][(j - 1)][0])) / mesh->le_lon[0][j][0]) / mesh->full_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dudtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u_lon[k][j][i] = (state->u_lon[k][j][i] + (dt * tend->smag_dudt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runu_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1198)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dvdt[k][j][i] = (state->kmh_lat[k][j][i] * ((((state->v_lat[k][j][(i - 1)] - (2 * state->v_lat[k][j][i])) + state->v_lat[k][j][(i + 1)]) / pow(mesh->le_lat[0][j][0],2)) + ((((((state->v_lat[k][(j + 1)][i] - state->v_lat[k][j][i]) / mesh->le_lon[0][(j + 1)][0]) * mesh->full_cos_lat[0][(j + 1)][0]) - (((state->v_lat[k][j][i] - state->v_lat[k][(j - 1)][i]) / mesh->le_lon[0][j][0]) * mesh->full_cos_lat[0][j][0])) / mesh->de_lat[0][j][0]) / mesh->half_cos_lat[0][j][0])));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dvdtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1198)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tend->smag_dvdt[k][j][i] = (state->kmh_lat[k][j][i] * (((state->v_lat[k][j][(i - 1)] - (2 * state->v_lat[k][j][i])) + state->v_lat[k][j][(i + 1)]) / pow(mesh->le_lat[0][j][0],2)));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runsmag_dvdtTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->v_lat[k][j][i] = (state->v_lat[k][j][i] + (dt * tend->smag_dvdt[k][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  smag_damp_runv_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] * state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  trickyptptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,0,&state->phs[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_2d_D(Proc, &state->phs[0][0][0], &Proc.FieldReq[async_phs], true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_phs]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_lev[k][j][i] = hybrid_coord_calc_ph_lev(mesh->hyai[k][0][0],mesh->hybi[k][0][0],state->phs[0][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph_exn_lev[k][j][i] = pow(state->ph_lev[k][j][i],0.2858964143426295);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phph_exn_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->ph_exn_lev[0][0][0], &Proc.FieldReq[async_ph_exn_lev], true, true, false, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->ph[k][j][i] = (0.5 * (state->ph_lev[k][j][i] + state->ph_lev[(k + 1)][j][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_phphTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m[k][j][i] = (state->ph_lev[(k + 1)][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mmTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->m[0][0][0], &Proc.FieldReq[async_m], true, true, true, false, true, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<1+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph[k][j][i] - state->ph_lev[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=32+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lev[k][j][i] = (state->ph_lev[k][j][i] - state->ph[(k - 1)][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  calc_mm_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lon[k][j][i] = ((state->m[k][j][i] + state->m[k][j][(i + 1)]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLonEdgem_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_lat[k][j][i] = ((state->m[k][j][i] + state->m[k][(j + 1)][i]) * 0.5);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  averageCellToLatEdgem_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->m_vtx[k][j][i] = ((((state->m[k][j][i] + state->m[k][j][(i + 1)]) * mesh->area_subcell_1[0][j][0]) + ((state->m[k][(j + 1)][i] + state->m[k][(j + 1)][(i + 1)]) * mesh->area_subcell_0[0][(j + 1)][0])) / mesh->area_vtx[0][j][0]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  interpCellToVtxm_vtxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(2,1,&state->pt[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_pt]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->pt[k][j][i] = (state->pt[k][j][i] / state->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  pole_damp_runptTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->pt[0][0][0], &Proc.FieldReq[async_pt], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lon_edge(3,&state->u_lon[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lon_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->u_lon[0][0][0], &Proc.FieldReq[async_u_lon], false, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_on_lat_edge(3,&state->v_lat[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_on_lat_edgeTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state->v_lat[0][0][0], &Proc.FieldReq[async_v_lat], true, false, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void c2a_2(struct HybridStateField* state)
{
  int k, j, i;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2399)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state->u[k][j][i] = (0.5 * (state->u_lon[k][j][i] + state->u_lon[k][j][(i - 1)]));
        state->v[k][j][i] = (0.5 * (state->v_lat[k][j][i] + state->v_lat[k][(j - 1)][i]));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  c2auTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void timeadvance_2()
{
  Time_Advance();
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Time_AdvanceTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void adv_run_2(struct HybridStateField* state_old, struct HybridStateField* state_new, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  double cursum, cf, s1, s2, ds1, ds2, ds3;
  int k, j, i, ci, l;
  double *tmp_state_old_tmpsum, *tmp_state_new_tmpsum;
  struct Vector3 buf;
  advPara = advPara;
  dt = dt;
  advPara = advPara;
  dt = dt;
  if (((advPara->uv_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = adv->u0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = adv->v0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->uv_step) = 1;
  }
  if (((advPara->uv_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->uu[k][j][i] = state_old->u_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celluuTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->vv[k][j][i] = state_old->v_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellvvTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->uv_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = ((adv->uu[k][j][i] + state_old->u_lon[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = ((adv->vv[k][j][i] + state_old->v_lat[k][j][i]) / ((advPara->nstep) + 1));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->u0[k][j][i] = adv->uu[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellu0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->v0[k][j][i] = adv->vv[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellv0Time += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_uu]);
      HaloWait(Proc,&Proc.FieldReq[async_u_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->uu[k][j][i] = (adv->uu[k][j][i] + state_old->u_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_celluuTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->uu[0][0][0], &Proc.FieldReq[async_uu], false, true, true, true, false, false, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_vv]);
      HaloWait(Proc,&Proc.FieldReq[async_v_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->vv[k][j][i] = (adv->vv[k][j][i] + state_old->v_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_uv_cellvvTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      UpdateHalo_3d_D(Proc, &adv->vv[0][0][0], &Proc.FieldReq[async_vv], true, false, true, false, false, true, false, true);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->uv_step) = 0;
  }
  else
  {
    (advPara->uv_step) = ((advPara->uv_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->uv_step) > (advPara->nstep))))
  {
    if (!(advPara->dynamic))
    {
      (advPara->uv_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflx[k][j][i] = ((adv->uu[k][j][i] * dt) / mesh->de_lon[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cfly[k][j][i] = ((adv->vv[k][j][i] * dt) / mesh->de_lat[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_cellcflyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_uu]);
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->divx[k][j][i] = (((adv->uu[k][j][i] - adv->uu[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]);
          adv->divy[k][j][i] = (((adv->vv[k][j][i] * mesh->le_lat[0][j][0]) - (adv->vv[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_vv]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state_old->tmpsum[k][j][i] = -adv->vv[k][(j - 1)][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celltmpsumTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    // SumCall:state_old->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
    tmp_state_old_tmpsum = allocate_1d_array_D(32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          tmp_state_old_tmpsum[k-Proc.lev_hw]=tmp_state_old_tmpsum[k-Proc.lev_hw]+state_old->tmpsum[k][j][i];
        }
      }
    }
    Zonal_Sum_1d_D(Proc, &tmp_state_old_tmpsum,32);
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          state_old->tmpsum[k][j][i]=tmp_state_old_tmpsum[k-Proc.lev_hw];
        }
      }
    }
    free_1d_array_D(tmp_state_old_tmpsum,32);
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->divy[k][j][i] = (((state_old->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_uv_celldivyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  advPara = advPara;
  if (((advPara->mf_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = 0.0;
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->mf_step) = 1;
  }
  if (((advPara->mf_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfx[k][j][i] = state_old->mfx_lon[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->mfy[k][j][i] = state_old->mfy_lat[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->mf_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = ((adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = ((adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]) / (advPara->nstep));
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfx_lon]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfx[k][j][i] = (adv->mfx[k][j][i] + state_old->mfx_lon[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfxTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_mfy_lat]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->mfy[k][j][i] = (adv->mfy[k][j][i] + state_old->mfy_lat[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_mf_cellmfyTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->mf_step) = 0;
  }
  else
  {
    (advPara->mf_step) = ((advPara->mf_step) + 1);
  }
  if ((!(advPara->dynamic) && ((advPara->mf_step) > (advPara->nstep))))
  {
    (advPara->mf_step) = -1;
  }
  advPara = advPara;
  dt = dt;
  if (((advPara->we_step) == -1))
  {
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = adv->we0[k][j][i];
          adv->mm[k][j][i] = adv->m0[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
    (advPara->we_step) = 1;
  }
  if (((advPara->we_step) == 0))
  {
    TimeBeg = MPI_Wtime();
    
    HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->we[k][j][i] = state_old->we_lev[k][j][i];
          adv->mm[k][j][i] = state_old->m_lev[k][j][i];
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levweTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  else
  {
    if (((advPara->we_step) == (advPara->nstep)))
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = ((adv->we[k][j][i] + state_old->we_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->mm[k][j][i] = ((adv->mm[k][j][i] + state_old->m_lev[k][j][i]) / ((advPara->nstep) + 1));
            adv->we0[k][j][i] = state_old->we_lev[k][j][i];
            adv->m0[k][j][i] = state_old->m_lev[k][j][i];
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
    else
    {
      TimeBeg = MPI_Wtime();
      
      HaloWait(Proc,&Proc.FieldReq[async_we_lev]);
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      for(k=0+Proc.lev_hw; k<33+Proc.lev_hw; k+=1){
        for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
          for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
            adv->we[k][j][i] = (adv->we[k][j][i] + state_old->we_lev[k][j][i]);
            adv->mm[k][j][i] = (adv->mm[k][j][i] + state_old->m_lev[k][j][i]);
          }
        }
      }
      TimeEnd = MPI_Wtime();
      CompTime += (TimeEnd - TimeBeg);
      accum_we_levweTime += (TimeEnd - TimeBeg);
      TimeBeg = MPI_Wtime();
      
      TimeEnd = MPI_Wtime();
      CommTime += (TimeEnd - TimeBeg);
      
      
    }
  }
  if ((advPara->dynamic))
  {
    (advPara->we_step) = 0;
  }
  else
  {
    (advPara->we_step) = ((advPara->we_step) + 1);
  }
  if (((advPara->dynamic) || ((advPara->we_step) > (advPara->nstep))))
  {
    if (((advPara->dynamic) == 0))
    {
      (advPara->we_step) = -1;
    }
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
      for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
        for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
          adv->cflz[k][j][i] = ((adv->we[k][j][i] / adv->mm[k][j][i]) * dt);
        }
      }
    }
    TimeEnd = MPI_Wtime();
    CompTime += (TimeEnd - TimeBeg);
    accum_we_levcflzTime += (TimeEnd - TimeBeg);
    TimeBeg = MPI_Wtime();
    
    TimeEnd = MPI_Wtime();
    CommTime += (TimeEnd - TimeBeg);
    
    
  }
  dt = dt;
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_old->qv[k][j][(i - 2)],state_old->qv[k][j][(i - 1)],state_old->qv[k][j][i],state_old->qv[k][j][(i + 1)],state_old->qv[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(state_old->qv[k][(j - 2)][i],state_old->qv[k][(j - 1)][i],state_old->qv[k][j][i],state_old->qv[k][(j + 1)][i],state_old->qv[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_uu]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + state_old->qv[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + state_old->qv[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->uu[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_vv]);
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->vv[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->vv[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_innerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->qx[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0]) / mesh->area_cell[0][j][0]) - (adv->divx[k][j][i] * state_old->qv[k][j][i]))) * dt));
        adv->qy[k][j][i] = (state_old->qv[k][j][i] - ((0.5 * ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) - (adv->divy[k][j][i] * state_old->qv[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_old->tmpsum[k][j][i] = adv->qmf_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state_old->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_state_old_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_old_tmpsum[k-Proc.lev_hw]=tmp_state_old_tmpsum[k-Proc.lev_hw]+state_old->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_old_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_old->tmpsum[k][j][i]=tmp_state_old_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_old_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->qx[k][j][i] = state_old->qv[k][j][i];
        adv->qy[k][j][i] = (state_old->qv[k][j][i] + ((0.5 * ((((state_old->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]) - (adv->divy[k][j][i] * state_old->qv[k][j][i]))) * dt));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  ffsl_calc_tracer_hflxqxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qx[0][0][0], &Proc.FieldReq[async_qx], true, true, true, false, false, true, true, true);
  UpdateHalo_3d_D(Proc, &adv->qy[0][0][0], &Proc.FieldReq[async_qy], true, true, true, true, true, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qx]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(adv->qy[k][j][(i - 2)],adv->qy[k][j][(i - 1)],adv->qy[k][j][i],adv->qy[k][j][(i + 1)],adv->qy[k][j][(i + 2)]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
        buf = ppm(adv->qx[k][(j - 2)][i],adv->qx[k][(j - 1)][i],adv->qx[k][j][i],adv->qx[k][(j + 1)][i],adv->qx[k][(j + 2)][i]);
        adv->qly[k][j][i] = buf.x;
        adv->dqy[k][j][i] = buf.y;
        adv->q6y[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->qly[0][0][0], &Proc.FieldReq[async_qly], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->dqy[0][0][0], &Proc.FieldReq[async_dqy], true, true, true, false, false, false, true, true);
  UpdateHalo_3d_D(Proc, &adv->q6y[0][0][0], &Proc.FieldReq[async_q6y], true, true, true, false, false, false, true, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qy]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflx[k][j][i]);
        cf = (adv->cflx[k][j][i] - ci);
        if ((adv->cflx[k][j][i] > 0.0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=((i + 1) - ci);l<(i + 1);l+=1){
            cursum = (cursum + adv->qy[k][j][l]);
          }
          adv->qmf_lon[k][j][i] = ((adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][(i - ci)] * ds1)) + ((0.5 * adv->dqx[k][j][(i - ci)]) * ds2)) + (adv->q6x[k][j][(i - ci)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
        }
        else
        {
          if ((adv->cflx[k][j][i] < 0.0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            cursum = 0.0;
            for (l=(i + 1);l<((i - ci) + 1);l+=1){
              cursum = (cursum + adv->qy[k][j][l]);
            }
            adv->qmf_lon[k][j][i] = ((-adv->mfx[k][j][i] * (((cursum + (adv->qlx[k][j][((i - ci) + 1)] * ds1)) + ((0.5 * adv->dqx[k][j][((i - ci) + 1)]) * ds2)) + (adv->q6x[k][j][((i - ci) + 1)] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflx[k][j][i]);
          }
          else
          {
            adv->qmf_lon[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_lonTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lon[0][0][0], &Proc.FieldReq[async_qmf_lon], false, true, true, true, false, false, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qly]);
  HaloWait(Proc,&Proc.FieldReq[async_dqy]);
  HaloWait(Proc,&Proc.FieldReq[async_q6y]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        if ((adv->cfly[k][j][i] > 0))
        {
          s1 = (1 - adv->cfly[k][j][i]);
          s2 = 1;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          adv->qmf_lat[k][j][i] = ((adv->mfy[k][j][i] * (((adv->qly[k][j][i] * ds1) + ((0.5 * adv->dqy[k][j][i]) * ds2)) + (adv->q6y[k][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
        }
        else
        {
          if ((adv->cfly[k][j][i] < 0))
          {
            s1 = 0;
            s2 = -adv->cfly[k][j][i];
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            adv->qmf_lat[k][j][i] = ((-adv->mfy[k][j][i] * (((adv->qly[k][(j + 1)][i] * ds1) + ((0.5 * adv->dqy[k][(j + 1)][i]) * ds2)) + (adv->q6y[k][(j + 1)][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cfly[k][j][i]);
          }
          else
          {
            adv->qmf_lat[k][j][i] = 0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  hflx_ppm_outerqmf_latTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qmf_lat[0][0][0], &Proc.FieldReq[async_qmf_lat], true, false, true, false, false, true, false, true);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lon]);
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1199)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (((adv->old_m[k][j][i] * state_old->qv[k][j][i]) - ((adv->qmf_lon[k][j][i] - adv->qmf_lon[k][j][(i - 1)]) * mesh->le_lon[0][j][0])) + ((((adv->qmf_lat[k][j][i] * mesh->le_lat[0][j][0]) - (adv->qmf_lat[k][(j - 1)][i] * mesh->le_lat[0][(j - 1)][0])) / mesh->area_cell[0][j][0]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qmf_lat]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->tmpsum[k][j][i] = adv->qmf_lat[k][(j - 1)][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runtmpsumTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  // SumCall:state_new->tmpsum.[0, 0, 32, 0, 1199, 1200, 1, 0, 2400]
  tmp_state_new_tmpsum = allocate_1d_array_D(32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        tmp_state_new_tmpsum[k-Proc.lev_hw]=tmp_state_new_tmpsum[k-Proc.lev_hw]+state_new->tmpsum[k][j][i];
      }
    }
  }
  Zonal_Sum_1d_D(Proc, &tmp_state_new_tmpsum,32);
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->tmpsum[k][j][i]=tmp_state_new_tmpsum[k-Proc.lev_hw];
      }
    }
  }
  free_1d_array_D(tmp_state_new_tmpsum,32);
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 1199)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = ((adv->old_m[k][j][i] * state_old->qv[k][j][i]) + ((((state_new->tmpsum[k][j][i] * mesh->le_lat[0][(j - 1)][0]) / 2400) / mesh->area_cell[0][j][0]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  buf = (struct Vector3){0.0,0.0,0.0};
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        buf = ppm(state_new->qv[(k - 2)][j][i],state_new->qv[(k - 1)][j][i],state_new->qv[k][j][i],state_new->qv[(k + 1)][j][i],state_new->qv[(k + 2)][j][i]);
        adv->qlx[k][j][i] = buf.x;
        adv->dqx[k][j][i] = buf.y;
        adv->q6x[k][j][i] = buf.z;
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqlxTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &adv->qlx[0][0][0], &Proc.FieldReq[async_qlx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->dqx[0][0][0], &Proc.FieldReq[async_dqx], true, true, true, true, true, false, false, false);
  UpdateHalo_3d_D(Proc, &adv->q6x[0][0][0], &Proc.FieldReq[async_q6x], true, true, true, true, true, false, false, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  cursum = 0.0;
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_qlx]);
  HaloWait(Proc,&Proc.FieldReq[async_dqx]);
  HaloWait(Proc,&Proc.FieldReq[async_q6x]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=1+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        ci = (int)(adv->cflz[k][j][i]);
        cf = (adv->cflz[k][j][i] - ci);
        if ((adv->cflz[k][j][i] > 0))
        {
          s1 = (1.0 - cf);
          s2 = 1.0;
          ds1 = (s2 - s1);
          ds2 = (pow(s2,2) - pow(s1,2));
          ds3 = (pow(s2,3) - pow(s1,3));
          cursum = 0.0;
          for (l=(k - ci);l<((k - 1) + 1);l+=1){
            cursum = (cursum + state_new->qv[l][j][i]);
          }
          adv->qmf_lev[k][j][i] = ((adv->we[k][j][i] * (((cursum + (adv->qlx[((k - ci) - 1)][j][i] * ds1)) + ((0.5 * adv->dqx[((k - ci) - 1)][j][i]) * ds2)) + (adv->q6x[((k - ci) - 1)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
        }
        else
        {
          if ((adv->cflz[k][j][i] < 0))
          {
            s1 = 0.0;
            s2 = -cf;
            ds1 = (s2 - s1);
            ds2 = (pow(s2,2) - pow(s1,2));
            ds3 = (pow(s2,3) - pow(s1,3));
            for (l=k;l<(((k - ci) - 1) + 1);l+=1){
              cursum = (cursum + state_new->qv[l][j][i]);
            }
            adv->qmf_lev[k][j][i] = ((-adv->we[k][j][i] * (((cursum + (adv->qlx[(k - ci)][j][i] * ds1)) + ((0.5 * adv->dqx[(k - ci)][j][i]) * ds2)) + (adv->q6x[(k - ci)][j][i] * ((ds2 / 2.0) - (ds3 / 3.0))))) / adv->cflz[k][j][i]);
          }
          else
          {
            adv->qmf_lev[k][j][i] = 0.0;
          }
        }
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  vflx_ppmqmf_levTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = ((state_new->qv[k][j][i] * state_old->m[k][j][i]) - ((adv->qmf_lev[(k + 1)][j][i] - adv->qmf_lev[k][j][i]) * 30));
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  Filter_On_Cell(1,1,&state_new->qv[0][0][0]);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  Filter_On_CellTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_qv]);
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        state_new->qv[k][j][i] = (state_new->qv[k][j][i] / state_old->m[k][j][i]);
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  adv_runqvTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  UpdateHalo_3d_D(Proc, &state_new->qv[0][0][0], &Proc.FieldReq[async_qv], true, true, true, true, true, true, true, false);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
  TimeBeg = MPI_Wtime();
  
  HaloWait(Proc,&Proc.FieldReq[async_m]);
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  for(k=0+Proc.lev_hw; k<32+Proc.lev_hw; k+=1){
    for(j=MAX(Proc.lat_beg, 0)-Proc.lat_beg+Proc.lat_hw; j<MIN(Proc.lat_end+1, 1200)-Proc.lat_beg+Proc.lat_hw; j+=1){
      for(i=MAX(Proc.lon_beg, 0)-Proc.lon_beg+Proc.lon_hw; i<MIN(Proc.lon_end+1, 2400)-Proc.lon_beg+Proc.lon_hw; i+=1){
        adv->old_m[k][j][i] = state_old->m[k][j][i];
      }
    }
  }
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  copy_old_mold_mTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void Diagnose_State_2(struct HybridStateField* state)
{
  Diagnose(state,15);
  TimeEnd = MPI_Wtime();
  CompTime += (TimeEnd - TimeBeg);
  DiagnoseTime += (TimeEnd - TimeBeg);
  TimeBeg = MPI_Wtime();
  
  TimeEnd = MPI_Wtime();
  CommTime += (TimeEnd - TimeBeg);
  
  
}

void MeshInit_0(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh)
{
  LatLonMeshInit(&global_mesh[0]);
  LatLonMeshInit_cp(&global_mesh[0], &mesh[0]);
}

void MeshInit_1(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh)
{
  LatLonMeshInit(&global_mesh[0]);
  LatLonMeshInit_cp(&global_mesh[0], &mesh[0]);
}

void MeshInit_2(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh)
{
  LatLonMeshInit(&global_mesh[0]);
  LatLonMeshInit_cp(&global_mesh[0], &mesh[0]);
}

void timeOperatorInit_0(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  ncIO_Read_0(&state[0], &staticv[0]);
  spaceOperatorInit_0(&staticv[0], &mesh[0]);
  gzlevInit_0(&state[Timeinfo.oldt], &staticv[0]);
  gzlevInit_0(&state[Timeinfo.newt], &staticv[0]);
  gzlevInit_0(&state[2], &staticv[0]);
  uvInit_0(&state[Timeinfo.oldt]);
  operatorPrepareNull_0(&state[Timeinfo.oldt], &staticv[0], &adv[0], advPara, &mesh[0], dt);
  advPrepare_0(&state[Timeinfo.oldt], &adv[1], &adv[0]);
  filterinit_0(&mesh[0], dt);
}

void timeOperatorInit_1(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  ncIO_Read_1(&state[0], &staticv[0]);
  spaceOperatorInit_1(&staticv[0], &mesh[0]);
  gzlevInit_1(&state[Timeinfo.oldt], &staticv[0]);
  gzlevInit_1(&state[Timeinfo.newt], &staticv[0]);
  gzlevInit_1(&state[2], &staticv[0]);
  uvInit_1(&state[Timeinfo.oldt]);
  operatorPrepareNull_1(&state[Timeinfo.oldt], &staticv[0], &adv[0], advPara, &mesh[0], dt);
  advPrepare_1(&state[Timeinfo.oldt], &adv[1], &adv[0]);
  filterinit_1(&mesh[0], dt);
}

void timeOperatorInit_2(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridAdvField* adv, struct HybridAdvPara* advPara, struct HybridMeshField* mesh, double dt)
{
  ncIO_Read_2(&state[0], &staticv[0]);
  spaceOperatorInit_2(&staticv[0], &mesh[0]);
  gzlevInit_2(&state[Timeinfo.oldt], &staticv[0]);
  gzlevInit_2(&state[Timeinfo.newt], &staticv[0]);
  gzlevInit_2(&state[2], &staticv[0]);
  uvInit_2(&state[Timeinfo.oldt]);
  operatorPrepareNull_2(&state[Timeinfo.oldt], &staticv[0], &adv[0], advPara, &mesh[0], dt);
  advPrepare_2(&state[Timeinfo.oldt], &adv[1], &adv[0]);
  filterinit_2(&mesh[0], dt);
}

void wrk3d_0(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridTendField* tend, struct HybridTendPara* tendPara, struct HybridAdvField* adv, struct HybridAdvPara* advptPara, struct HybridAdvPara* advmPara, struct HybridMeshField* mesh, double dt, double dtd2, double dtd3)
{
  stepForwardBackward_0(&state[Timeinfo.oldt], &state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd3);
  stepForwardBackward_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[2], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd2);
  stepForwardBackward_0(&state[Timeinfo.oldt], &state[2], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dt);
  damp_run_0(&state[Timeinfo.newt], &tend[Timeinfo.newt], &mesh[0], dt);
  c2a_0(&state[Timeinfo.newt]);
  timeadvance_0();
  adv_run_0(&state[Timeinfo.oldt], &state[Timeinfo.newt], &adv[1], advmPara, &mesh[0], dt);
  Diagnose_State_0(&state[Timeinfo.oldt]);
}

void wrk3d_1(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridTendField* tend, struct HybridTendPara* tendPara, struct HybridAdvField* adv, struct HybridAdvPara* advptPara, struct HybridAdvPara* advmPara, struct HybridMeshField* mesh, double dt, double dtd2, double dtd3)
{
  stepForwardBackward_1(&state[Timeinfo.oldt], &state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd3);
  stepForwardBackward_1(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[2], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd2);
  stepForwardBackward_1(&state[Timeinfo.oldt], &state[2], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dt);
  damp_run_1(&state[Timeinfo.newt], &tend[Timeinfo.newt], &mesh[0], dt);
  c2a_1(&state[Timeinfo.newt]);
  timeadvance_1();
  adv_run_1(&state[Timeinfo.oldt], &state[Timeinfo.newt], &adv[1], advmPara, &mesh[0], dt);
  Diagnose_State_1(&state[Timeinfo.oldt]);
}

void wrk3d_2(struct HybridStateField* state, struct HybridStaticField* staticv, struct HybridTendField* tend, struct HybridTendPara* tendPara, struct HybridAdvField* adv, struct HybridAdvPara* advptPara, struct HybridAdvPara* advmPara, struct HybridMeshField* mesh, double dt, double dtd2, double dtd3)
{
  stepForwardBackward_2(&state[Timeinfo.oldt], &state[Timeinfo.oldt], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd3);
  stepForwardBackward_2(&state[Timeinfo.oldt], &state[Timeinfo.newt], &state[2], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dtd2);
  stepForwardBackward_2(&state[Timeinfo.oldt], &state[2], &state[Timeinfo.newt], &staticv[0], &tend[Timeinfo.oldt], &tend[Timeinfo.newt], tendPara, &adv[0], advptPara, &mesh[0], dt);
  damp_run_2(&state[Timeinfo.newt], &tend[Timeinfo.newt], &mesh[0], dt);
  c2a_2(&state[Timeinfo.newt]);
  timeadvance_2();
  adv_run_2(&state[Timeinfo.oldt], &state[Timeinfo.newt], &adv[1], advmPara, &mesh[0], dt);
  Diagnose_State_2(&state[Timeinfo.oldt]);
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
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 5){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,337, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 6 && Proc.lat_end <= 11){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,33, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 12 && Proc.lat_end <= 17){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,18, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 18 && Proc.lat_end <= 23){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,12, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 24 && Proc.lat_end <= 29){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,10, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 30 && Proc.lat_end <= 35){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,8, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 36 && Proc.lat_end <= 41){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,7, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 42 && Proc.lat_end <= 47){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,6, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 48 && Proc.lat_end <= 59){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,5, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 60 && Proc.lat_end <= 83){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,4, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 84 && Proc.lat_end <= 143){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,3, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 144 && Proc.lat_end <= 1055){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,2, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1056 && Proc.lat_end <= 1115){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,3, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1116 && Proc.lat_end <= 1139){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,4, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1140 && Proc.lat_end <= 1151){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,5, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1152 && Proc.lat_end <= 1157){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,6, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1158 && Proc.lat_end <= 1163){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,7, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1164 && Proc.lat_end <= 1169){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,8, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1170 && Proc.lat_end <= 1175){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,10, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1176 && Proc.lat_end <= 1181){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,13, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1182 && Proc.lat_end <= 1187){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,19, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1188 && Proc.lat_end <= 1193){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,36, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  if (Proc.lat_beg >= 1194 && Proc.lat_end <= 1199){
    ProcInit_LonLat_Ngb(&Proc,MPI_COMM_WORLD,ProcLatNum,ProcLon,ProcLat,nlon,nlat,nlev,337, 2, 2, 2, 2, DTYPE_DOUBLE, 108);
  }
  
  PhysicalVariableInit();
  
  TimeInit();
  
  if (Proc.lat_beg >= 0 && Proc.lat_end <= 5)
  {
    MeshInit_0(global_mesh, mesh);
    timeOperatorInit_0(state, staticv, adv, &advptPara, mesh, 30);
    tottime_beg = MPI_Wtime();
    for (int t = 0; t < 30; t += 1){
      wrk3d_0(state, staticv, tend, &tendPara, adv, &advptPara, &advmPara, mesh, 30, 15.0, 10.0);
    }
    tottime_end = MPI_Wtime();
  }
  
  if (Proc.lat_beg >= 6 && Proc.lat_end <= 1193)
  {
    MeshInit_1(global_mesh, mesh);
    timeOperatorInit_1(state, staticv, adv, &advptPara, mesh, 30);
    tottime_beg = MPI_Wtime();
    for (int t = 0; t < 30; t += 1){
      wrk3d_1(state, staticv, tend, &tendPara, adv, &advptPara, &advmPara, mesh, 30, 15.0, 10.0);
    }
    tottime_end = MPI_Wtime();
  }
  
  if (Proc.lat_beg >= 1194 && Proc.lat_end <= 1199)
  {
    MeshInit_2(global_mesh, mesh);
    timeOperatorInit_2(state, staticv, adv, &advptPara, mesh, 30);
    tottime_beg = MPI_Wtime();
    for (int t = 0; t < 30; t += 1){
      wrk3d_2(state, staticv, tend, &tendPara, adv, &advptPara, &advmPara, mesh, 30, 15.0, 10.0);
    }
    tottime_end = MPI_Wtime();
  }
  
  ProfilingOutPut(Proc.id);
  if (Proc.id == 0) printf("Time is %.2f\n",tottime_end-tottime_beg);
}