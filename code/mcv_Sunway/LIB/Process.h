#ifndef PROCESS_H_INCLUDED
#define PROCESS_H_INCLUDED 1

#include"mpi.h"
#include<stdbool.h>

#include"ParaParam.h"
#include"Async.h"

//邻居
typedef struct{
  int id;     //进程号
  int orient; //方位 以接收方halo位置为准
  int side;   //本面的方位,以保留发送halo在本面内的位置
  int pos;    //CubedSphere的位置 0,1 中,边角
  int rotate; //CubedSphere的数据接收后是否翻转 0 不 1 翻
  int type;   //收发
  int tag;    //消息tag tag设置方式为接收方进程号+接收方halo方位编号
  int lon_beg,lon_end,lat_beg,lat_end; //邻居和本进程交互数据范围，若type=send则为发送范围，反之为接收范围
}NgbType; 

//各halo发送时的subarray信息
typedef struct{
  int host_id,ngb_id;
  int orient;
  MPI_Datatype dtype;
  int tag;
  int req; //CubedSphere的request号
  MPI_Datatype mpi_type_2d[2][2];
  MPI_Datatype mpi_type_3d[2][2][2];

}HaloType;

//本进程
typedef struct {
  MPI_Comm comm;
  int size;
  int id;
  int nlev,full_nlev,half_nlev;
  //LonLat
  int nlon,nlat;
  int lon_beg,lon_end,lat_beg,lat_end; //进程范围
  int lon_num,lat_num; //经纬长度
  int full_nlon,half_nlon,full_nlat,half_nlat; //在网格不同位置的三维长度
  bool at_north_pole,at_south_pole; //是否在极点处
  int lon_hw,lat_hw,lev_hw;
  //CubedSphere
  int panel,panelrank,xpid,ypid; //面号[1-6],面内进程号;xy两方向上的序号,用于方便定位邻面进程号
  int ncell;
  int x_beg,x_end,y_beg,y_end; //进程范围,水平方向为x
  int x_num,y_num;
  int full_nx,half_nx,full_ny,half_ny;
  int x_hw,y_hw,p_hw;
  int at_south,at_north,at_east,at_west;//是否在面边界处
  int re_west,re_east,re_north,re_south; //边界处halo是否翻转
  int NgbSendNum, NgbRecvNum ; //收发邻居数
  int NgbSSP, NgbSNP, NgbRSP, NgbRNP ; //send//recv south_pos north_pos , 收发的南北邻居的起始位置 LonLat
  int NgbWS,NgbES,NgbSS,NgbNS,NgbWR,NgbER,NgbSR,NgbNR; //东南西北个数 CubedShpere 代表各方向存在halo时需要做的收发操作
  //LonLat 下标从1开始 即 [1,NgbSendNum] 存放顺序为 西 东 recv(南 北) send(北 南) 
  NgbType  NgbSend[ngb_max];
  NgbType  NgbRecv[ngb_max];
  HaloType SendHalo[ngb_max];
  HaloType RecvHalo[ngb_max];
  //CubedSphere 下标从0开 东南西北全部拆开 方向同样为接收方halo方位
  NgbType  NgbSendW[ngb_cs],NgbSendE[ngb_cs],NgbSendS[ngb_cs],NgbSendN[ngb_cs];
  NgbType  NgbRecvW[ngb_cs],NgbRecvE[ngb_cs],NgbRecvS[ngb_cs],NgbRecvN[ngb_cs];
  HaloType SendHaloW[ngb_cs],SendHaloE[ngb_cs],SendHaloS[ngb_cs],SendHaloN[ngb_cs];
  HaloType RecvHaloW[ngb_cs],RecvHaloE[ngb_cs],RecvHaloS[ngb_cs],RecvHaloN[ngb_cs];
  AsyncType *FieldReq; 
  int *CSWait;
  //纬向通信子
  MPI_Comm zonal_comm;

}ProcType;


ProcType Proc;


//预处理进程基本信息
void ProcInit(ProcType *proc, MPI_Comm comm);

//预处理进程基本信息，经纬网格，传入划分方式，返回各进程计算区域
void ProcInit_LonLat_Domain(ProcType *Proc, MPI_Comm comm,int* ProcLatNum ,int* ProcLon ,int* ProcLat ,int NumLon, int NumLat, int NumLev);
void ProcInit_LonLat_Ngb(ProcType *Proc, MPI_Comm comm,int* ProcLatNum ,int* ProcLon ,int* ProcLat ,int NumLon, int NumLat, int NumLev, int HaloSizeLon, int HaloSizeLat, int HaloSizeLev, int HaloSizeLonSelf,  int CornerHalo, int dtype, int FieldNum);

//进程相关功能函数
bool Has_South_Pole();
bool Has_North_Pole();

//Cubed-Sphere
void ProcInit_CubedSphere_Domain(ProcType *Proc, MPI_Comm comm,int size, int rank, int* ProcLatNum ,int* ProcX ,int* ProcY ,int NumCell, int NumLev);
void ProcInit_CubedSphere_Ngb(ProcType *Proc, MPI_Comm comm,int size, int id, int ProcX ,int ProcY ,int NumCell, int NumLev, int HaloX, int HaloY, int HalosizeLev, int HaloP, int dtype, int FieldNum);
#endif