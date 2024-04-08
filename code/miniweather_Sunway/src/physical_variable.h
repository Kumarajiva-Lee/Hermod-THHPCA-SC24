#ifndef PHYSICAL_VARIABLE_H_INCLUDED
#define PHYSICAL_VARIABLE_H_INCLUDED 1

#include<stdbool.h>
struct HybridMeshField{ 
  double*** xp;
  double*** zpf;
  double*** zph;
};
struct MeshVector6{
  int lx;
  int ly;
  int lz;
  int sx;
  int sy;
  int sz;
};
struct MeshVector3{
    int x;
    int y;
    int z;
};
struct HybridStateField{ 
  double*** u;
  double*** w;
  double*** rho;
  double*** pt;
};
struct HybridStaticField{ 
  double*** hy_dens_cell;
  double*** hy_dens_theta_cell;
  double*** hy_dens_int;
  double*** hy_dens_theta_int;
  double*** hy_pressure_int;
};
struct HybridDirPara{
  bool xd;
  bool zd;
  int stepnum;
};
struct InitVector{
  double r;
  double u;
  double w;
  double t;
  double hr;
  double ht;
};
struct HybridFluxField{ 
  double*** fu;
  double*** fw;
  double*** frho;
  double*** fpt;
};
struct HybridTendField{ 
  double*** du;
  double*** dw;
  double*** drho;
  double*** dpt;
};
struct Vector4{
      double r;
      double u;
      double w;
      double t;
};
struct HybridOutField{ 
  double*** dens;
  double*** uwnd;
  double*** wwnd;
  double*** theta;
};
struct HybridOutPara{
  int step;
  int interval;
};
struct Vector2{
  double x;
  double y;
};

struct HybridMeshField mesh[1];
struct HybridMeshField global_mesh[1];
struct HybridStateField state[2];
struct HybridStaticField staticv[1];
struct HybridStaticField global_staticv[1];
struct HybridDirPara DirPara;
struct HybridOutField stateout[1];
struct HybridTendField tend[1];
struct HybridFluxField flux[1];
struct HybridOutPara OutPara;

#define async_u 1
#define async_w 2
#define async_rho 3
#define async_pt 4
#define async_dens 5
#define async_uwnd 6
#define async_wwnd 7
#define async_theta 8
#define async_du 9
#define async_dw 10
#define async_drho 11
#define async_dpt 12
#define async_fu 13
#define async_fw 14
#define async_frho 15
#define async_fpt 16

extern double TimeBeg;
extern double TimeEnd;
extern double CommTime;
extern double CompTime;
extern double MeshInitTime;
extern double StateInitialrhoTime;
extern double StateInitialTime;
extern double ComputeTendenciesXfrhoTime;
extern double ComputeTendenciesXdrhoTime;
extern double SetBoundaryZrhoTime;
extern double SetBoundaryZuTime;
extern double SetBoundaryZwTime;
extern double SetBoundaryZptTime;
extern double ComputeTendenciesZfrhoTime;
extern double ComputeTendenciesZdrhoTime;
extern double UpdateStaterhoTime;
extern double OutputPreparedensTime;
extern double ncDefDimTime;
extern double ncCloseFileTime;

void ProfilingOutPut(int id);

void PhysicalVariableInit();
void PhysicalVariableFinish();

void MeshInit_cp(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh);
#endif
