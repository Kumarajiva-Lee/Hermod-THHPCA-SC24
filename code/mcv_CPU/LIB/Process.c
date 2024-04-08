#include<stdio.h>
#include <stdlib.h>

#include"mpi.h"
#include "Process.h"
#include "ParaParam.h"

void ProcInit(ProcType *Proc, MPI_Comm comm){
    MPI_Comm_rank(MPI_COMM_WORLD, &Proc->id);
}


//LonLat
void NgbTypeInitLL(NgbType *Ngb,int id, int type, int orient, int tag, int lonbeg, int lonend, int latbeg, int latend){
    Ngb->id = id;
    Ngb->type = type;
    Ngb->orient = orient;
    Ngb->tag = tag;
    Ngb->lon_beg = lonbeg;
    Ngb->lon_end = lonend;
    Ngb->lat_beg = latbeg;
    Ngb->lat_end = latend;
}

//LonLat
void HaloTypeInitLL(HaloType *Halo, ProcType *Proc, int HaloType, int orient, int dtype, int host_id, int ngb_id, int tag, \
                  int lon_beg, int lon_end, int lat_beg, int lat_end, int NumLon, int NumLat, int NumLev, int HaloSizeLon, int HaloSizeLat){
    
    int full_lon_ibeg, full_lon_iend;
    int full_lat_ibeg, full_lat_iend;
    int half_lon_ibeg, half_lon_iend;
    int half_lat_ibeg, half_lat_iend;
    int num_full_lon, num_half_lon, num_full_lat, num_half_lat;
    int array_size[2][2][3];
    int subarray_size[2][2][3];
    int subarray_start[2][2][3];
    int array_size_2d[2];
    int subarray_size_2d[2];
    int subarray_start_2d[2];
    int num_lev[2];
    int i,j,k;

    if (dtype == DTYPE_INT){
        #ifdef halodtype
            #undef halodtype
        #endif  
        #define halodtype MPI_INT
    }
    else if(dtype == DTYPE_SINGLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_FLOAT
    }
    else if(dtype == DTYPE_DOUBLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_DOUBLE   
    }


    Halo->dtype = halodtype;
    Halo->host_id = host_id;
    Halo->ngb_id = ngb_id;
    Halo->tag = tag;
    Halo->orient = orient;
    
    //ToCheck 为了处理简单,full half都加上了halo
    num_lev[0] = NumLev + 2 * Proc->lev_hw;     //full   
    num_lev[1] = NumLev + 1 + 2 * Proc->lev_hw; //half

    num_full_lon = Proc->full_nlon;
    num_half_lon = Proc->half_nlon;
    num_full_lat = Proc->full_nlat;
    num_half_lat = Proc->half_nlat;

    full_lon_ibeg = lon_beg;
    full_lon_iend = lon_end;
    half_lon_ibeg = lon_beg;
    half_lon_iend = lon_end;
    full_lat_ibeg = lat_beg;
    full_lat_iend = lat_end;
    half_lat_ibeg = lat_beg;
    half_lat_iend = lat_end;
    //半网格在北极点少1圈，故北halo范围集体前移1圈
    if (Proc->at_north_pole && orient == NORTH){
        half_lat_ibeg --;
        half_lat_iend --;
    }
    else if (Proc->at_north_pole && (orient == WEST || orient == EAST)){
        half_lat_iend --;
    }



    for (k = 0; k < 2; k ++){
        array_size[0][0][0] = num_lev[k]; array_size[0][0][1] = num_full_lat+2*HaloSizeLat; array_size[0][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[0][1][0] = num_lev[k]; array_size[0][1][1] = num_full_lat+2*HaloSizeLat; array_size[0][1][2] = num_half_lon+2*HaloSizeLon;
        array_size[1][0][0] = num_lev[k]; array_size[1][0][1] = num_half_lat+2*HaloSizeLat; array_size[1][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[1][1][0] = num_lev[k]; array_size[1][1][1] = num_half_lat+2*HaloSizeLat; array_size[1][1][2] = num_half_lon+2*HaloSizeLon;

        subarray_size[0][0][0] = num_lev[k]; subarray_size[0][0][1] = full_lat_iend-full_lat_ibeg+1; subarray_size[0][0][2] = full_lon_iend-full_lon_ibeg+1;
        subarray_size[0][1][0] = num_lev[k]; subarray_size[0][1][1] = full_lat_iend-full_lat_ibeg+1; subarray_size[0][1][2] = half_lon_iend-half_lon_ibeg+1;
        subarray_size[1][0][0] = num_lev[k]; subarray_size[1][0][1] = half_lat_iend-half_lat_ibeg+1; subarray_size[1][0][2] = full_lon_iend-full_lon_ibeg+1;
        subarray_size[1][1][0] = num_lev[k]; subarray_size[1][1][1] = half_lat_iend-half_lat_ibeg+1; subarray_size[1][1][2] = half_lon_iend-half_lon_ibeg+1;

        subarray_start[0][0][0] = 0; subarray_start[0][0][1] = full_lat_ibeg; subarray_start[0][0][2] = full_lon_ibeg;
        subarray_start[0][1][0] = 0; subarray_start[0][1][1] = full_lat_ibeg; subarray_start[0][1][2] = half_lon_ibeg;
        subarray_start[1][0][0] = 0; subarray_start[1][0][1] = half_lat_ibeg; subarray_start[1][0][2] = full_lon_ibeg;
        subarray_start[1][1][0] = 0; subarray_start[1][1][1] = half_lat_ibeg; subarray_start[1][1][2] = half_lon_ibeg;

        for (j = 0; j < 2; j++)
            for (i = 0; i < 2; i++){
                //单层模式,不存在half_lat
                if (full_lat_iend-full_lat_ibeg+1 == 1 && j == 1) continue;
                if (dtype == DTYPE_INT){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_SINGLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_DOUBLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_3d[k][j][i]));
                }
                MPI_Type_commit(&(Halo->mpi_type_3d[k][j][i]));
            } 
    }

    for (j = 0; j < 2; j++)
        for (i = 0; i < 2; i++){
            //单层模式,不存在half_lat
            if (full_lat_iend-full_lat_ibeg+1 == 1 && j == 1) continue;

            array_size_2d[0] = array_size[j][i][1];
            array_size_2d[1] = array_size[j][i][2];

            subarray_size_2d[0] = subarray_size[j][i][1];
            subarray_size_2d[1] = subarray_size[j][i][2];

            subarray_start_2d[0] = subarray_start[j][i][1];
            subarray_start_2d[1] = subarray_start[j][i][2];

            if (dtype == DTYPE_INT){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_SINGLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_DOUBLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_2d[j][i]));
            }
            MPI_Type_commit(&(Halo->mpi_type_2d[j][i]));
        }
}

void ProcInit_LonLat_Domain(ProcType *Proc, MPI_Comm comm,int* ProcLatNum ,int* ProcLon ,int* ProcLat ,int NumLon, int NumLat, int NumLev){

    int Proccurr,latcurr,loncurr;
    int i,p,plon,plat,len,lencurr;
    int lonlen,latlen;
    int id,size;


    p = 0;
    Proccurr = 0;
    latcurr = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &Proc->size);
    Proc->comm = comm;
    Proc->id = id;
    Proc->nlon = NumLon;
    Proc->nlat = NumLat;
    Proc->nlev = NumLev;

    // //求进程负责范围
    // while (Proc->id < Proccurr){
    //     Proccurr += ProcLon[p] * ProcLat[p];
    //     p ++;
    // }

    // for (i = 0 ; i < p ; i++) latcurr += ProcLatNum[i];
    // Proccurr = Proc->id - Proccurr;

    plat = Proc->id % ProcLat[p];
    plon = Proc->id / ProcLat[p];



    loncurr = 0;
    lencurr = NumLon;

    for (i = 0 ; i < plon ; i ++){
        len = lencurr / (ProcLon[p] - i);
        if (lencurr % (ProcLon[p] - i)) len ++;
        loncurr += len;
        lencurr -= len;
    }

    len = lencurr / (ProcLon[p] - plon);
    if (lencurr % (ProcLon[p] - plon)) len ++;
    Proc->lon_beg = loncurr;
    Proc->lon_end = Proc->lon_beg + len - 1;
    lonlen = len;

    lencurr = ProcLatNum[p];
    for (i = 0 ; i < plat ; i ++){
        len = lencurr / (ProcLat[p] - i);
        if (lencurr % (ProcLat[p] - i)) len ++;
        latcurr += len;
        lencurr -= len;
    }

    len = lencurr / (ProcLat[p] - plat);
    if (lencurr % (ProcLat[p] - plat)) len ++;
    Proc->lat_beg = latcurr;
    Proc->lat_end = Proc->lat_beg + len - 1;
    latlen = len;  

    Proc->lon_num = lonlen;
    Proc->lat_num = latlen;

    Proc->full_nlon = Proc->lon_num;
    Proc->full_nlat = Proc->lat_num;
    Proc->full_nlev = NumLev;

    Proc->half_nlon = Proc->full_nlon;
    Proc->half_nlat = Proc->full_nlat - Has_North_Pole();
    Proc->half_nlev = Proc->full_nlev + 1;
}

void ProcInit_LonLat_Ngb(ProcType *Proc, MPI_Comm comm,int* ProcLatNum ,int* ProcLon ,int* ProcLat ,int NumLon, int NumLat, int NumLev, int HaloSizeLon, int HaloSizeLat, int HaloSizeLev, int HaloSizeLonSelf, int CornerHalo, int dtype, int FieldNum){

    int Proccurr,latcurr,loncurr;
    int i,p,plon,plat,len,lencurr;
    int lonlen,latlen;
    int wid,eid,nid,sid,id,pid;
    NgbType ngb;
    
    if (HaloSizeLon > HaloSizeLonSelf) Proc->lon_hw = HaloSizeLon;
    else Proc->lon_hw = HaloSizeLonSelf;
    Proc->lat_hw = HaloSizeLat;
    Proc->lev_hw = HaloSizeLev;

    if (HaloSizeLon < HaloSizeLonSelf) HaloSizeLonSelf = HaloSizeLon;

    p = 0;
    Proccurr = 0;
    latcurr = 0;

    id = Proc->id;


    Proc->FieldReq = (AsyncType*)malloc(FieldNum * sizeof(AsyncType));

    //求进程负责范围
    while (Proc->id < Proccurr){
        Proccurr += ProcLon[p] * ProcLat[p];
        p ++;
    }

    for (i = 0 ; i < p ; i++) latcurr += ProcLatNum[i];
    Proccurr = Proc->id - Proccurr;

    plat = Proccurr % ProcLat[p];
    plon = Proccurr / ProcLat[p];

    loncurr = 0;
    lencurr = NumLon;

    for (i = 0 ; i < plon ; i ++){
        len = lencurr / (ProcLon[p] - i);
        if (lencurr % (ProcLon[p] - i)) len ++;
        loncurr += len;
        lencurr -= len;
    }

    len = lencurr / (ProcLon[p] - plon);
    if (lencurr % (ProcLon[p] - plon)) len ++;
    lonlen = len;

    lencurr = ProcLatNum[p];
    for (i = 0 ; i < plat ; i ++){
        len = lencurr / (ProcLat[p] - i);
        if (lencurr % (ProcLat[p] - i)) len ++;
        latcurr += len;
        lencurr -= len;
    }

    len = lencurr / (ProcLat[p] - plat);
    if (lencurr % (ProcLat[p] - plat)) len ++;
    latlen = len;  


    //求邻居关系
    Proc->at_north_pole = false;
    Proc->at_south_pole = false;
    
    if (HaloSizeLon > 0){
        //west
        if (plon == 0){
            //west边界
            wid = id + ProcLat[p] * (ProcLon[p]-1);
        }
        else{
            //正常区域
            wid = id - ProcLat[p];
        };
        NgbTypeInitLL(&Proc->NgbRecv[WEST], wid, recv_type, WEST,  id+WEST, 0       , HaloSizeLon-1  , HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbSend[EAST], wid, send_type, EAST, wid+EAST, HaloSizeLon, 2*HaloSizeLon-1, HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbRecv[WEST_S], wid, recv_type, WEST,  id+WEST, HaloSizeLon-HaloSizeLonSelf, HaloSizeLon-1  , HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbSend[EAST_S], wid, send_type, EAST, wid+EAST, HaloSizeLon, HaloSizeLon+HaloSizeLonSelf-1, HaloSizeLat, HaloSizeLat+latlen-1);

        //east
        if (plon == ProcLon[p]-1) {
            //east边界
            eid = id - ProcLat[p] * (ProcLon[p]-1);
        }
        else{
            //正常区域
            eid = id + ProcLat[p];
        }
        NgbTypeInitLL(&Proc->NgbRecv[EAST], eid, recv_type, EAST,  id+EAST, HaloSizeLon+lonlen, 2*HaloSizeLon+lonlen-1, HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbSend[WEST], eid, send_type, WEST, eid+WEST, lonlen         , HaloSizeLon+lonlen-1  , HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbRecv[EAST_S], eid, recv_type, EAST,  id+EAST, HaloSizeLon+lonlen, HaloSizeLon+lonlen+HaloSizeLonSelf-1, HaloSizeLat, HaloSizeLat+latlen-1);
        NgbTypeInitLL(&Proc->NgbSend[WEST_S], eid, send_type, WEST, eid+WEST, lonlen+HaloSizeLon-HaloSizeLonSelf, HaloSizeLon+lonlen-1  , HaloSizeLat, HaloSizeLat+latlen-1);

        //南北邻居个数不固定，预处理位置关系
        Proc->NgbSendNum = 2+2;
        Proc->NgbRecvNum = 2+2;
        Proc->NgbSNP = 3+2;
        Proc->NgbRSP = 3+2;
    }


    //south 
    if (HaloSizeLat > 0){
        if (Proc->lat_beg == 0 && Proc->full_nlat >=3){
            //南极
            Proc->at_south_pole = true;
            if (ProcLon[p] == 1){
                pid = Proc->id;
            }
            else{
                pid = Proc->id - ProcLon[p] / 2 * ProcLat[p];
                if (pid < 0) pid = Proc->id + ProcLon[p] / 2 * ProcLat[p];
            }
            NgbTypeInitLL(&Proc->NgbRecv[OPPOSITE], pid, recv_type, SOUTH, id+SOUTH, HaloSizeLon, HaloSizeLon+lonlen-1, 0, HaloSizeLat-1);
            NgbTypeInitLL(&Proc->NgbSend[OPPOSITE], pid, send_type, SOUTH, id+SOUTH, HaloSizeLon, HaloSizeLon+lonlen-1, HaloSizeLat + 1, 2*HaloSizeLat);
        }
        else if(plat == 0){
            //在划分边界处，需处理非规则邻居
        }
        else{
            //在规则区域内，按8邻居处理
            //西南
            pid = wid - 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], pid, recv_type, SOUTH, id+SOUTH+WEST , HaloSizeLon-1 - CornerHalo + 1, HaloSizeLon-1, 0, HaloSizeLat-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], pid, send_type, NORTH, pid+NORTH+EAST, HaloSizeLon, HaloSizeLon + CornerHalo - 1, HaloSizeLat, 2*HaloSizeLat-1);
            //正南
            sid = id - 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], sid, recv_type, SOUTH, id+SOUTH , HaloSizeLon, HaloSizeLon+lonlen-1, 0, HaloSizeLat-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], sid, send_type, NORTH, sid+NORTH, HaloSizeLon, HaloSizeLon+lonlen-1, HaloSizeLat, 2*HaloSizeLat-1);
            //东南
            pid = eid - 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], pid, recv_type, SOUTH, id+SOUTH+EAST , HaloSizeLon+lonlen, HaloSizeLon+lonlen + CornerHalo -1, 0, HaloSizeLat-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], pid, send_type, NORTH, pid+NORTH+WEST, HaloSizeLon+lonlen - CornerHalo        , HaloSizeLon+lonlen-1, HaloSizeLat, 2*HaloSizeLat-1);
        
        }
        
        Proc->NgbSSP = Proc->NgbSendNum + 1;
        Proc->NgbRNP = Proc->NgbRecvNum + 1;
    }

    //north
    if (HaloSizeLat > 0){
        if (Proc->lat_end == NumLat - 1 && Proc->full_nlat >=3){
            //北极
            Proc->at_north_pole = true;
            if (ProcLon[p] == 1){
                pid = Proc->id;
            }
            else{
                pid = Proc->id + ProcLon[p] / 2 * ProcLat[p];
                if (pid >= Proc->size) pid = Proc->id - ProcLon[p] / 2 * ProcLat[p];
            }
            NgbTypeInitLL(&Proc->NgbRecv[OPPOSITE], pid, recv_type, NORTH, id+NORTH, HaloSizeLon, HaloSizeLon+lonlen-1, HaloSizeLat+latlen, 2*HaloSizeLat+latlen-1);
            NgbTypeInitLL(&Proc->NgbSend[OPPOSITE], pid, send_type, NORTH, id+NORTH, HaloSizeLon, HaloSizeLon+lonlen-1, latlen-1, latlen+HaloSizeLat-2);

        }
        else if(plat == ProcLat[p] - 1){
            //在划分边界处，需处理非规则邻居
        }
        else{
            //在规则区域内，按8邻居处理
            //西北
            pid = wid + 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], pid, recv_type, NORTH, id+NORTH+WEST , HaloSizeLon-1 - CornerHalo + 1, HaloSizeLon-1, HaloSizeLat+latlen, 2*HaloSizeLat+latlen-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], pid, send_type, SOUTH, pid+SOUTH+EAST, HaloSizeLon, HaloSizeLon + CornerHalo - 1, latlen, latlen+HaloSizeLat-1);
            //正北
            nid = id + 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], nid, recv_type, NORTH, id+NORTH , HaloSizeLon, HaloSizeLon+lonlen-1, HaloSizeLat+latlen, 2*HaloSizeLat+latlen-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], nid, send_type, SOUTH, nid+SOUTH, HaloSizeLon, HaloSizeLon+lonlen-1, latlen, latlen+HaloSizeLat-1);
            //东北
            pid = eid + 1;
            NgbTypeInitLL(&Proc->NgbRecv[++Proc->NgbRecvNum], pid, recv_type, NORTH, id+NORTH+EAST , HaloSizeLon+lonlen, HaloSizeLon+lonlen + CornerHalo -1, HaloSizeLat+latlen, 2*HaloSizeLat+latlen-1);
            NgbTypeInitLL(&Proc->NgbSend[++Proc->NgbSendNum], pid, send_type, SOUTH, pid+SOUTH+WEST, HaloSizeLon+lonlen - CornerHalo         , HaloSizeLon+lonlen-1, latlen, latlen+HaloSizeLat-1); 
        }
    }

    //设置Halo的Subarray信息
    for (i = 0; i <= Proc->NgbRecvNum; i++){
        ngb = Proc->NgbRecv[i];
        HaloTypeInitLL(&Proc->RecvHalo[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, NumLon, NumLat, NumLev,HaloSizeLon, HaloSizeLat);
    }
    for (i = 0; i <= Proc->NgbSendNum; i++){
        ngb = Proc->NgbSend[i];
        HaloTypeInitLL(&Proc->SendHalo[i], Proc, send_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, NumLon, NumLat, NumLev,HaloSizeLon, HaloSizeLat);
    }

    //预处理async_req
    for (i = 0 ; i < FieldNum ; i++){
        for (p = 0 ; p <= Proc->NgbSendNum ; p++) 
            Proc->FieldReq[i].send_req[p] = MPI_REQUEST_NULL;
        for (p = 0 ; p <= Proc->NgbRecvNum ; p++) 
            Proc->FieldReq[i].recv_req[p] = MPI_REQUEST_NULL;
    }

    //设置纬向通信子

    int color;
    
    color = Proc->lat_beg;

    MPI_Comm_split(MPI_COMM_WORLD, color, comm, &Proc->zonal_comm);

}

bool Has_South_Pole(){
    return Proc.lat_beg == 0;
}

bool Has_North_Pole(){
    return Proc.lat_end == (Proc.nlat - 1);
}

bool Has_Boundary_X(){
    return Proc.lon_end == (Proc.x_num - 1);
}

bool Has_Boundary_Y(){
    return Proc.lat_end == (Proc.y_num - 1);
}

//CubedShpere

void NgbTypeInitCS(NgbType *Ngb,int id, int type, int orient, int side, int rotate, int tag, int lonbeg, int lonend, int latbeg, int latend){
    Ngb->id = id;
    Ngb->type = type;
    Ngb->orient = orient;
    Ngb->side   = side;
    Ngb->rotate = rotate;
    Ngb->tag = tag;
    Ngb->lon_beg = lonbeg;
    Ngb->lon_end = lonend;
    Ngb->lat_beg = latbeg;
    Ngb->lat_end = latend;
}


//pos = 0 中心 pos = 1 边角
void NgbTypeInitCS_Z(NgbType *Ngb,int id, int type, int orient, int pos, int rotate, int tag, int lonbeg, int lonend, int latbeg, int latend){
    Ngb->id = id;
    Ngb->type = type;
    Ngb->orient = orient;
    Ngb->pos = pos;
    Ngb->rotate = rotate;
    Ngb->tag = tag;
    Ngb->lon_beg = lonbeg;
    Ngb->lon_end = lonend;
    Ngb->lat_beg = latbeg;
    Ngb->lat_end = latend;
}

void HaloTypeInitCS(HaloType *Halo, ProcType *Proc, int HaloType, int orient, int side, int dtype, int host_id, int ngb_id, int tag, \
                  int lon_beg, int lon_end, int lat_beg, int lat_end, int NumLon, int NumLat, int NumLev, int HaloSizeLon, int HaloSizeLat, int req){
    
    int full_lon_beg, full_lon_end;
    int full_lat_beg, full_lat_end;
    int half_lon_beg, half_lon_end;
    int half_lat_beg, half_lat_end;
    int num_full_lon, num_half_lon, num_full_lat, num_half_lat;
    int array_size[2][2][3];
    int subarray_size[2][2][3];
    int subarray_start[2][2][3];
    int array_size_2d[2];
    int subarray_size_2d[2];
    int subarray_start_2d[2];
    int num_lev[2];
    int i,j,k;

    if (dtype == DTYPE_INT){
        #ifdef halodtype
            #undef halodtype
        #endif  
        #define halodtype MPI_INT
    }
    else if(dtype == DTYPE_SINGLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_FLOAT
    }
    else if(dtype == DTYPE_DOUBLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_DOUBLE   
    }


    Halo->dtype = halodtype;
    Halo->host_id = host_id;
    Halo->ngb_id = ngb_id;
    Halo->tag = tag;
    Halo->orient = orient;
    Halo->req = req;
    
    //ToCheck 为了处理简单,full half都加上了halo
    num_lev[0] = NumLev + 2 * Proc->lev_hw;     //full   
    num_lev[1] = NumLev + 1 + 2 * Proc->lev_hw; //half

    num_full_lon = NumLon;
    num_half_lon = NumLon;
    num_full_lat = NumLat;
    num_half_lat = NumLat;

    if (Proc->at_north) num_half_lat --;
    if (Proc->at_east) num_half_lon --;

    full_lon_beg = lon_beg;
    full_lon_end = lon_end;
    half_lon_beg = lon_beg;
    half_lon_end = lon_end;
    full_lat_beg = lat_beg;
    full_lat_end = lat_end;
    half_lat_beg = lat_beg;
    half_lat_end = lat_end;

    if (Proc->at_west){
        //左侧,全网格发送+1
        if (HaloType == send_type  && side == WEST){ //&& orient == EAST
            full_lon_beg ++;
            full_lon_end ++;
        }
    }
    if (Proc->at_east){
        //右侧,全网格发送-1
        if (HaloType == send_type  && side == EAST){ //&& orient == WEST
            full_lon_beg --;
            full_lon_end --;
        }
        //半网格 -1
        if (side == EAST){
            half_lon_beg --;
            half_lon_end --;
        }
        else if (side == SOUTH || side == NORTH){
            half_lon_end --;
        }
    }
    if (Proc->at_south){
        //下侧,全网格发送+1
        if (HaloType == send_type && side == SOUTH){ //&& orient == NORTH 
            full_lat_beg ++;
            full_lat_end ++;
        }
    }
    if (Proc->at_north){
        //上侧,全网格发送-1
        if (HaloType == send_type && side == NORTH){ //&& orient == SOUTH 
            full_lat_beg --;
            full_lat_end --;
        }
        //半网格-1
        if (side == NORTH){
            half_lat_beg --;
            half_lat_end --;
        }
        else if (side == WEST || side == EAST){
            half_lat_end --;
        }
    }

    // if (Proc->at_east){
    //     num_half_lon++;
    //     half_lon_iend++;
    // }
    // if (Proc->at_north){
    //     num_half_lat++;
    //     half_lat_iend++;
    // }

    for (k = 0; k < 2; k ++){
        array_size[0][0][0] = num_lev[k]; array_size[0][0][1] = num_full_lat+2*HaloSizeLat; array_size[0][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[0][1][0] = num_lev[k]; array_size[0][1][1] = num_full_lat+2*HaloSizeLat; array_size[0][1][2] = num_half_lon+2*HaloSizeLon;
        array_size[1][0][0] = num_lev[k]; array_size[1][0][1] = num_half_lat+2*HaloSizeLat; array_size[1][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[1][1][0] = num_lev[k]; array_size[1][1][1] = num_half_lat+2*HaloSizeLat; array_size[1][1][2] = num_half_lon+2*HaloSizeLon;
        subarray_size[0][0][0] = num_lev[k]; subarray_size[0][0][1] = full_lat_end-full_lat_beg+1; subarray_size[0][0][2] = full_lon_end-full_lon_beg+1;
        subarray_size[0][1][0] = num_lev[k]; subarray_size[0][1][1] = full_lat_end-full_lat_beg+1; subarray_size[0][1][2] = half_lon_end-half_lon_beg+1;
        subarray_size[1][0][0] = num_lev[k]; subarray_size[1][0][1] = half_lat_end-half_lat_beg+1; subarray_size[1][0][2] = full_lon_end-full_lon_beg+1;
        subarray_size[1][1][0] = num_lev[k]; subarray_size[1][1][1] = half_lat_end-half_lat_beg+1; subarray_size[1][1][2] = half_lon_end-half_lon_beg+1;

        subarray_start[0][0][0] = 0; subarray_start[0][0][1] = full_lat_beg; subarray_start[0][0][2] = full_lon_beg;
        subarray_start[0][1][0] = 0; subarray_start[0][1][1] = full_lat_beg; subarray_start[0][1][2] = half_lon_beg;
        subarray_start[1][0][0] = 0; subarray_start[1][0][1] = half_lat_beg; subarray_start[1][0][2] = full_lon_beg;
        subarray_start[1][1][0] = 0; subarray_start[1][1][1] = half_lat_beg; subarray_start[1][1][2] = half_lon_beg;

        for (j = 0; j < 2; j++)
            for (i = 0; i < 2; i++){
                //if (Proc->id == 1) //printf("%d %d %d start [%d %d %d]\n",k,j,i,subarray_start[j][i][0],subarray_start[j][i][1],subarray_start[j][i][2]);
                if (dtype == DTYPE_INT){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_SINGLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_DOUBLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_3d[k][j][i]));
                }
                MPI_Type_commit(&(Halo->mpi_type_3d[k][j][i]));
            } 
    }

    for (j = 0; j < 2; j++)
        for (i = 0; i < 2; i++){
         
            array_size_2d[0] = array_size[j][i][1];
            array_size_2d[1] = array_size[j][i][2];

            subarray_size_2d[0] = subarray_size[j][i][1];
            subarray_size_2d[1] = subarray_size[j][i][2];

            subarray_start_2d[0] = subarray_start[j][i][1];
            subarray_start_2d[1] = subarray_start[j][i][2];

            if (dtype == DTYPE_INT){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_SINGLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_DOUBLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_2d[j][i]));
            }
            MPI_Type_commit(&(Halo->mpi_type_2d[j][i]));
        }
    

    return;

}

void HaloTypeInitCS_Z(HaloType *Halo, ProcType *Proc, int HaloType, int orient, int dtype, int host_id, int ngb_id, int tag, \
                  int lon_beg, int lon_end, int lat_beg, int lat_end, int NumLon, int NumLat, int NumLev, int HaloSizeLon, int HaloSizeLat, int req){
    
    int full_lon_ibeg, full_lon_iend;
    int full_lat_ibeg, full_lat_iend;
    int half_lon_ibeg, half_lon_iend;
    int half_lat_ibeg, half_lat_iend;
    int num_full_lon, num_half_lon, num_full_lat, num_half_lat;
    int array_size[2][2][3];
    int subarray_size[2][2][3];
    int subarray_start[2][2][3];
    int array_size_2d[2];
    int subarray_size_2d[2];
    int subarray_start_2d[2];
    int num_lev[2];
    int i,j,k;

    if (dtype == DTYPE_INT){
        #ifdef halodtype
            #undef halodtype
        #endif  
        #define halodtype MPI_INT
    }
    else if(dtype == DTYPE_SINGLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_FLOAT
    }
    else if(dtype == DTYPE_DOUBLE){
        #ifdef halodtype
            #undef halodtype
        #endif
        #define halodtype MPI_DOUBLE   
    }


    Halo->dtype = halodtype;
    Halo->host_id = host_id;
    Halo->ngb_id = ngb_id;
    Halo->tag = tag;
    Halo->orient = orient;
    Halo->req = req;
    
    //ToCheck 为了处理简单,full half都加上了halo
    num_lev[0] = NumLev + 2 * Proc->lev_hw;     //full   
    num_lev[1] = NumLev + 1 + 2 * Proc->lev_hw; //half

    num_full_lon = NumLon;
    num_half_lon = NumLon;
    num_full_lat = NumLat;
    num_half_lat = NumLat;

    full_lon_ibeg = lon_beg;
    full_lon_iend = lon_end;
    half_lon_ibeg = lon_beg;
    half_lon_iend = lon_end;
    full_lat_ibeg = lat_beg;
    full_lat_iend = lat_end;
    half_lat_ibeg = lat_beg;
    half_lat_iend = lat_end;

    // if (Proc->at_east){
    //     num_half_lon++;
    //     half_lon_iend++;
    // }
    // if (Proc->at_north){
    //     num_half_lat++;
    //     half_lat_iend++;
    // }

    for (k = 0; k < 2; k ++){
        array_size[0][0][0] = num_lev[k]; array_size[0][0][1] = num_full_lat+2*HaloSizeLat; array_size[0][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[0][1][0] = num_lev[k]; array_size[0][1][1] = num_full_lat+2*HaloSizeLat; array_size[0][1][2] = num_half_lon+2*HaloSizeLon;
        array_size[1][0][0] = num_lev[k]; array_size[1][0][1] = num_half_lat+2*HaloSizeLat; array_size[1][0][2] = num_full_lon+2*HaloSizeLon;
        array_size[1][1][0] = num_lev[k]; array_size[1][1][1] = num_half_lat+2*HaloSizeLat; array_size[1][1][2] = num_half_lon+2*HaloSizeLon;
        subarray_size[0][0][0] = num_lev[k]; subarray_size[0][0][1] = full_lat_iend-full_lat_ibeg+1; subarray_size[0][0][2] = full_lon_iend-full_lon_ibeg+1;
        subarray_size[0][1][0] = num_lev[k]; subarray_size[0][1][1] = full_lat_iend-full_lat_ibeg+1; subarray_size[0][1][2] = half_lon_iend-half_lon_ibeg+1;
        subarray_size[1][0][0] = num_lev[k]; subarray_size[1][0][1] = half_lat_iend-half_lat_ibeg+1; subarray_size[1][0][2] = full_lon_iend-full_lon_ibeg+1;
        subarray_size[1][1][0] = num_lev[k]; subarray_size[1][1][1] = half_lat_iend-half_lat_ibeg+1; subarray_size[1][1][2] = half_lon_iend-half_lon_ibeg+1;

        subarray_start[0][0][0] = 0; subarray_start[0][0][1] = full_lat_ibeg; subarray_start[0][0][2] = full_lon_ibeg;
        subarray_start[0][1][0] = 0; subarray_start[0][1][1] = full_lat_ibeg; subarray_start[0][1][2] = half_lon_ibeg;
        subarray_start[1][0][0] = 0; subarray_start[1][0][1] = half_lat_ibeg; subarray_start[1][0][2] = full_lon_ibeg;
        subarray_start[1][1][0] = 0; subarray_start[1][1][1] = half_lat_ibeg; subarray_start[1][1][2] = half_lon_ibeg;

        for (j = 0; j < 2; j++)
            for (i = 0; i < 2; i++){
                if (dtype == DTYPE_INT){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_SINGLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_3d[k][j][i]));
                }
                else if (dtype == DTYPE_DOUBLE){
                    MPI_Type_create_subarray(3, array_size[j][i], subarray_size[j][i], subarray_start[j][i], MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_3d[k][j][i]));
                }
                MPI_Type_commit(&(Halo->mpi_type_3d[k][j][i]));
            } 
    }

    for (j = 0; j < 2; j++)
        for (i = 0; i < 2; i++){
            array_size_2d[0] = array_size[j][i][1];
            array_size_2d[1] = array_size[j][i][2];

            subarray_size_2d[0] = subarray_size[j][i][1];
            subarray_size_2d[1] = subarray_size[j][i][2];

            subarray_start_2d[0] = subarray_start[j][i][1];
            subarray_start_2d[1] = subarray_start[j][i][2];

            if (dtype == DTYPE_INT){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_INT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_SINGLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_FLOAT, &(Halo->mpi_type_2d[j][i]));
            }
            else if (dtype == DTYPE_DOUBLE){
                MPI_Type_create_subarray(2, array_size_2d, subarray_size_2d, subarray_start_2d, MPI_ORDER_C, MPI_DOUBLE, &(Halo->mpi_type_2d[j][i]));
            }
            MPI_Type_commit(&(Halo->mpi_type_2d[j][i]));
        }
    

    return;

}

//默认cell可以整除划分数
void ProcInit_CubedSphere_Domain(ProcType *Proc, MPI_Comm comm,int size, int id, int* ProcLatNum ,int* ProcX ,int* ProcY ,int NumCell, int NumLev){
    int panelproc;

    Proc->comm = comm;
    Proc->id = id;
    Proc->size = size;
    Proc->ncell = NumCell;
    Proc->nlev = NumLev;

    if (size != 1 && size != ProcX[0] * ProcY[0] * 6){
        //printf("Wrong Proc Size! size = %d , Proc In = %d * %d * 6\n",size,ProcX[0],ProcY[0]);
        return;
    }
    else if (ProcLatNum[0] != NumCell)
    {
        //printf("Wrong Cell Size!\n");
        return;
    }


    panelproc = ProcX[0] * ProcY[0];
    Proc->panel = id / panelproc + 1;
    Proc->panelrank = id % panelproc;

    Proc->xpid = Proc->panelrank / ProcY[0];
    Proc->ypid = Proc->panelrank % ProcY[0];

    Proc->x_num = NumCell / ProcX[0];
    Proc->y_num = NumCell / ProcY[0];

    Proc->x_beg = Proc->x_num * Proc->xpid;
    Proc->x_end = Proc->x_num * (Proc->xpid + 1) - 1;

    Proc->y_beg = Proc->y_num * Proc->ypid;
    Proc->y_end = Proc->y_num * (Proc->ypid + 1) - 1;

    //ToCheck 与经纬保持一致,半网格在最后-1
    Proc->full_nx = Proc->x_num;
    if (Proc->xpid != ProcX[0] - 1) Proc->half_nx = Proc->full_nx;
    else Proc->half_nx = Proc->full_nx - 1;

    Proc->full_ny = Proc->y_num;
    if (Proc->ypid != ProcY[0] - 1) Proc->half_ny = Proc->full_ny;
    else Proc->half_ny = Proc->full_ny - 1;

    Proc->full_nlev = NumLev;
    Proc->half_nlev = NumLev + 1;

    Proc->at_south = (Proc->ypid == 0)?true:false;
    Proc->at_north = (Proc->ypid == ProcY[0] - 1)?true:false;
    Proc->at_west  = (Proc->xpid == 0)?true:false;
    Proc->at_east  = (Proc->xpid == ProcX[0] - 1)?true:false;

    //简化生成
    Proc->lon_beg = Proc->x_beg;
    Proc->lon_end = Proc->x_end;
    Proc->lat_beg = Proc->y_beg;
    Proc->lat_end = Proc->y_end;
    Proc->full_nlon = Proc->full_nx;
    Proc->half_nlon = Proc->half_nx;
    Proc->full_nlat = Proc->full_ny;
    Proc->half_nlat = Proc->half_ny;


}

//暂时支持6*N*N
void ProcInit_CubedSphere_Ngb(ProcType *Proc, MPI_Comm comm,int size, int id, int ProcX ,int ProcY ,int NumCell, int NumLev, int HaloX, int HaloY, int HalosizeLev, int HaloP, int dtype, int FieldNum){
    int i,j,k,p;
    int wid,eid,nid,sid,pid;
    int nx,ny;
    int pn; //panelNgb
    NgbType ngb;

    Proc->x_hw = HaloX;
    Proc->y_hw = HaloY;
    Proc->lon_hw = HaloX;
    Proc->lat_hw = HaloY;

    nx = Proc->x_num;
    ny = Proc->y_num;

    Proc->NgbSendNum = 0;
    Proc->NgbRecvNum = 0;

    Proc->NgbWS = 0;
    Proc->NgbES = 0;
    Proc->NgbSS = 0;
    Proc->NgbNS = 0;
    Proc->NgbWR = 0;
    Proc->NgbER = 0;
    Proc->NgbSR = 0;
    Proc->NgbNR = 0;

    if (Proc->panel == 1){
        Proc->re_west  = 0;
        Proc->re_east  = 0;
        Proc->re_north = 0;
        Proc->re_south = 0;
    }
    else if (Proc->panel == 2){
        Proc->re_west  = 0;
        Proc->re_east  = 0;
        Proc->re_north = 0;
        Proc->re_south = 1;
    }
    else if (Proc->panel == 3){
        Proc->re_west  = 0;
        Proc->re_east  = 0;
        Proc->re_north = 1;
        Proc->re_south = 1;
    }
    else if (Proc->panel == 4){
        Proc->re_west  = 0;
        Proc->re_east  = 0;
        Proc->re_north = 1;
        Proc->re_south = 0;
    }
    else if (Proc->panel == 5){
        Proc->re_west  = 1;
        Proc->re_east  = 0;
        Proc->re_north = 1;
        Proc->re_south = 0;
    }   
    else{//panel == 6
        Proc->re_west  = 0;
        Proc->re_east  = 1;
        Proc->re_north = 0;
        Proc->re_south = 1;
    }

    //WEST_Ngb 处理西边邻居关系
    if (! Proc->at_west){
        //不在边界
        wid = id - ProcY;
        NgbTypeInitCS(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, WEST, 0, id+WEST, 0, HaloX-1, HaloY, HaloY+ny-1);
        NgbTypeInitCS(&Proc->NgbSendE[Proc->NgbES++], wid, send_type, EAST, WEST, 0, wid+EAST, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);
    }
    else{
        //在不同面
        if (Proc->panel <= 4){
            if (Proc->panel > 1) wid = id - ProcY;
            else wid = id + 4 * ProcX * ProcY - ProcY;
            NgbTypeInitCS(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, WEST, 0, id+WEST, 0, HaloX-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendE[Proc->NgbES++], wid, send_type, EAST, WEST, 0, wid+EAST, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);
        }
        else if (Proc->panel == 5){
            wid = 3 * ProcX * ProcY -1 + (ProcY - Proc->ypid) * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, WEST, 1, id+WEST, 0, HaloX-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], wid, send_type, NORTH, WEST, 0, wid+NORTH, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);
        } 
        else{//panel == 6
            wid = 3 * ProcX * ProcY -1 + Proc->ypid * ProcY + 1;
            NgbTypeInitCS(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, WEST, 0, id+WEST, 0, HaloX-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], wid, send_type, SOUTH, WEST, 0, wid+SOUTH, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);
        }
    }

    //EAST_Ngb 处理东边邻居关系
    if (!Proc->at_east){
        //不在边界
        eid = id + ProcY;
        NgbTypeInitCS(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, EAST, 0, id+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY, HaloY+ny-1);
        NgbTypeInitCS(&Proc->NgbSendW[Proc->NgbWS++], eid, send_type, WEST, EAST, 0, eid+WEST, nx, nx+HaloX-1, HaloY, HaloY+ny-1);
    }
    else{
        if (Proc->panel <= 4){
            if (Proc->panel <= 3) eid = id + ProcY;
            else eid = Proc->ypid;
            NgbTypeInitCS(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, EAST, 0, id+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendW[Proc->NgbWS++], eid, send_type, WEST, EAST, 0, eid+WEST, nx, nx+HaloX-1, HaloY, HaloY+ny-1);
        }
        else if (Proc->panel == 5){
            eid = ProcX * ProcY - 1 + (Proc->ypid + 1) * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, EAST, 0, id+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], eid, send_type, NORTH, EAST, 0, eid+NORTH, nx, nx+HaloX-1, HaloY, HaloY+ny-1);
        }
        else{//panel == 6
            eid = ProcX * ProcY - 1 + (ProcY - Proc->ypid - 1) * ProcY + 1;
            NgbTypeInitCS(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, EAST, 0, id+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY, HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], eid, send_type, SOUTH, EAST, 1, eid+SOUTH, nx, nx+HaloX-1, HaloY, HaloY+ny-1);
        }
    }
    
    //SOUTH_Ngb 处理南边邻居关系
    if (!Proc->at_south){
        //不在边界
        sid = id - 1;
        NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
        NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], sid, send_type, NORTH, SOUTH, 0, sid+NORTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
    }
    else{
        if (Proc->panel == 1){
            sid = 5 * ProcX * ProcY - 1 + (Proc->xpid + 1) * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], sid, send_type, NORTH, SOUTH,  0, sid+NORTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
        else if (Proc->panel == 2){
            sid = 6 * ProcX * ProcY - 1 - Proc->xpid;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendE[Proc->NgbES++], sid, send_type, EAST, SOUTH, 1, sid+EAST, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
        else if (Proc->panel == 3){
            sid = 5 * ProcX * ProcY + (ProcX - Proc->xpid - 1)*ProcY;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], sid, send_type, SOUTH, SOUTH, 1, sid+SOUTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
        else if (Proc->panel == 4){
            sid = 5 * ProcX * ProcY + Proc->xpid;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendW[Proc->NgbWS++], sid, send_type, WEST, SOUTH, 0, sid+WEST, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
        else if (Proc->panel == 5){
            sid = (Proc->xpid + 1) * ProcY - 1;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], sid, send_type, NORTH, SOUTH, 0, sid+NORTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
        else{//panel == 6
            sid = 2 * ProcX * ProcY + (ProcX - Proc->xpid - 1) * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, SOUTH, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], sid, send_type, SOUTH, SOUTH, 1, sid+SOUTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
        }
    }

    //NORTH_Ngb 处理南边邻居关系
    if (!Proc->at_north){
        nid = id + 1;
        NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
        NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], nid, send_type, SOUTH, NORTH,  0, nid+SOUTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
    }
    else{
        if (Proc->panel == 1){
            nid = 4 * ProcX * ProcY + Proc->xpid * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], nid, send_type, SOUTH, NORTH,  0, nid+SOUTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
        else if (Proc->panel == 2){
            nid = 4 * ProcX * ProcY + (ProcX-1)*ProcY + Proc->xpid;
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendE[Proc->NgbES++], nid, send_type, EAST, NORTH,  0, nid+EAST, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
        else if (Proc->panel == 3){
            nid = 4 * ProcX * ProcY - 1 + (ProcX - Proc->xpid)*ProcY;
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], nid, send_type, NORTH, NORTH,  1, nid+NORTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
        else if (Proc->panel == 4){
            nid = 4 * ProcX * ProcY - 1 + (ProcX - Proc->xpid);
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendW[Proc->NgbWS++], nid, send_type, WEST, NORTH,  1, nid+WEST, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
        else if (Proc->panel == 5){
            nid = 2 * ProcX * ProcY - 1 + (ProcX - Proc->xpid) * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendN[Proc->NgbNS++], nid, send_type, NORTH, NORTH,  1, nid+NORTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
        else{//panel == 6
            nid = Proc->xpid * ProcY;
            NgbTypeInitCS(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, NORTH,  0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
            NgbTypeInitCS(&Proc->NgbSendS[Proc->NgbSS++], nid, send_type, SOUTH, NORTH,  0, nid+SOUTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
        }
    }

  

    //Ngb转Halo
    for (k = 0 ; k < size ; k++){
        if (id == k){
            //printf("=============================\n");
            //printf("ID %d Panel %d Pid %d [%d %d %d %d]\n", Proc->id, Proc->panel, Proc->panelrank, Proc->at_west, Proc->at_east, Proc->at_south, Proc->at_north);
            //Send
            //printf("Send\n");
            for (i = 0; i < Proc->NgbWS; i++){
                ngb = Proc->NgbSendW[i];
                HaloTypeInitCS(&Proc->SendHaloW[i], Proc, send_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
                //printf("West ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            } 
            for (i = 0; i < Proc->NgbES; i++){
                ngb = Proc->NgbSendE[i];
                HaloTypeInitCS(&Proc->SendHaloE[i], Proc, send_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
                //printf("East ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            for (i = 0; i < Proc->NgbSS; i++){
                ngb = Proc->NgbSendS[i];
                HaloTypeInitCS(&Proc->SendHaloS[i], Proc, send_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
                //printf("South ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            for (i = 0; i < Proc->NgbNS; i++){
                ngb = Proc->NgbSendN[i];
                HaloTypeInitCS(&Proc->SendHaloN[i], Proc, send_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                     ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
                //printf("North ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            
            //Recv
            //printf("Recv\n");
            for (i = 0; i < Proc->NgbWR; i++){
                ngb = Proc->NgbRecvW[i];
                HaloTypeInitCS(&Proc->RecvHaloW[i], Proc, recv_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                    ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
                //printf("West ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            for (i = 0; i < Proc->NgbER; i++){
                ngb = Proc->NgbRecvE[i];
                HaloTypeInitCS(&Proc->RecvHaloE[i], Proc, recv_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                    ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
                //printf("East ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            for (i = 0; i < Proc->NgbSR; i++){
                ngb = Proc->NgbRecvS[i];
                HaloTypeInitCS(&Proc->RecvHaloS[i], Proc, recv_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                    ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
                //printf("South ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
            for (i = 0; i < Proc->NgbNR; i++){
                ngb = Proc->NgbRecvN[i];
                HaloTypeInitCS(&Proc->RecvHaloN[i], Proc, recv_type, ngb.orient, ngb.side, dtype, id, ngb.id, ngb.tag, \
                    ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
                //printf("North ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
            }
        }
        // MPI_Barrier(comm);
        // sleep(2);
        MPI_Barrier(comm);
    }

    //预处理async_req
    Proc->FieldReq = (AsyncType*)malloc(FieldNum * sizeof(AsyncType));

    for (i = 0 ; i < FieldNum ; i++){
        for (p = 1 ; p <= Proc->NgbSendNum ; p++) 
            Proc->FieldReq[i].send_req[p] = MPI_REQUEST_NULL;
        for (p = 1 ; p <= Proc->NgbRecvNum ; p++) 
            Proc->FieldReq[i].recv_req[p] = MPI_REQUEST_NULL;
    }

    
    //预处理CSWait
    Proc->CSWait = (int*)malloc(FieldNum * sizeof(int));

}

//Z字展开
// void ProcInit_CubedSphere_Ngb_Z(ProcType *Proc, MPI_Comm comm,int size, int id, int ProcX ,int ProcY ,int NumCell, int NumLev, int HaloX, int HaloY, int HalosizeLev, int HaloP, int dtype, int FieldNum){
//     int i,j,k,p;
//     int wid,eid,nid,sid,pid;
//     int nx,ny;
//     int pn; //panelNgb
//     NgbType ngb;

//     Proc->x_hw = HaloX;
//     Proc->y_hw = HaloY;
//     Proc->p_hw = HaloP;
//     Proc->lon_hw = HaloX;
//     Proc->lat_hw = HaloY;

//     nx = Proc->x_num;
//     ny = Proc->y_num;

//     Proc->NgbSendNum = 0;
//     Proc->NgbRecvNum = 0;

//     Proc->NgbWS = 0;
//     Proc->NgbES = 0;
//     Proc->NgbSS = 0;
//     Proc->NgbSS = 0;
//     Proc->NgbWR = 0;
//     Proc->NgbER = 0;
//     Proc->NgbSR = 0;
//     Proc->NgbSR = 0;

//     //WEST_Ngb 处理西边邻居关系
//     if (! Proc->at_west){
//         //不在边界
//         wid = id - ProcY;
//         NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, 0, 0, id+WEST, 0, HaloX-1, HaloY, HaloY+ny-1);
//         NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], wid, send_type, EAST, 0, 0, wid+EAST, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);
//     }
//     else{
//         //在边,分奇偶
//         if (Proc->panel & 1){
//             //奇

//             pn = (Proc->panel - 1 + 4) % 6;
//             wid = (ProcY-1-Proc->ypid+1) * ProcY - 1 + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, 1, 1, id+WEST, HaloP, HaloP+ny-1, 0, HaloX-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], wid, send_type, NORTH, 1, 0, wid+NORTH, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);

//             if (Proc->at_south){
//                 //在左下角               
//                 pid = wid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 1, id+WEST+NORTH, 0, HaloP-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+EAST, HaloX, 2*HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//             }
//             else if (Proc->at_north){
//                 //在左上角
//                 pid = wid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 1, id+WEST+SOUTH, HaloP+ny, 2*HaloP+ny-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+WEST, HaloX, 2*HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//             else{
//                 //普通边
//                 pid = wid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 1, id+WEST+NORTH, 0, HaloP-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+EAST, HaloX, 2*HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//                 pid = wid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 1, id+WEST+SOUTH, HaloP+ny, 2*HaloP+ny-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+WEST, HaloX, 2*HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//         }
//         else{
//             //偶

//             pn = (Proc->panel - 1 + 5) % 6;
//             wid = (ProcX - 1)*ProcY + Proc->ypid + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], wid, recv_type, WEST, 1, 0, id+WEST, 0, HaloX-1, HaloP, HaloP+ny-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], wid, send_type, EAST, 1, 0, wid+EAST, HaloX, 2*HaloX-1, HaloY, HaloY+ny-1);

//             if (Proc->at_south){
//                 //在左下角
//                 pid = wid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 0, id+WEST+NORTH, 0, HaloX-1, HaloP+ny, 2*HaloP+ny-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+SOUTH, HaloX, 2*HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//             }
//             else if (Proc->at_north){
//                 //在左上角
//                 pid = wid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 0, id+WEST+SOUTH, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+NORTH, HaloX, 2*HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//             else{
//                 //普通边
//                 pid = wid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 0, id+WEST+NORTH, 0, HaloX-1, HaloP+ny, 2*HaloP+ny-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+SOUTH, HaloX, 2*HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//                 pid = wid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvW[Proc->NgbWR++], pid, recv_type, WEST, 1, 0, id+WEST+SOUTH, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+NORTH, HaloX, 2*HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//         }
//     }

//     //EAST_Ngb
//     if (! Proc->at_east){
//         //不在边界
//         eid = id + ProcY;
//         NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, 0, 0, id+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY, HaloY+ny-1);
//         NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], eid, send_type, WEST, 0, 0, eid+WEST, nx, nx+HaloX-1, HaloY, HaloY+ny-1);
//     }
//     else{
//         //在边,分奇偶
//         if (Proc->panel & 1){
//             //奇

//             pn = (Proc->panel - 1 + 1) % 6;
//             eid = Proc->ypid + pn*ProcX*ProcY;
//             NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, 1, 0, id+EAST, 0, HaloX-1, HaloP, HaloP+ny-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], eid, send_type, WEST, 1, 0, eid+WEST, nx, nx+HaloX-1, HaloY, HaloY+ny-1);

//             if (Proc->at_south){
//                 //右下角
//                 pid = eid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 0, id+EAST+NORTH, 0, HaloX-1, HaloP+ny, 2*HaloP+ny-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+SOUTH, nx, nx+HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//             }
//             else if (Proc->at_north){
//                 //右上角
//                 pid = eid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 0, id+EAST+SOUTH, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid ,send_type, WEST, 1, 0, pid+WEST+NORTH, nx, nx+HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//             else{
//                 //普通边
//                 pid = eid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 0, id+EAST+NORTH, 0, HaloX-1, HaloP+ny, 2*HaloP+ny-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+SOUTH, nx, nx+HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//                 pid = eid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 0, id+EAST+SOUTH, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid ,send_type, WEST, 1, 0, pid+WEST+NORTH, nx, nx+HaloX-1, HaloY, HaloY+HaloP-1);

//             }
//         }
//         else{
//             //偶
//             pn = (Proc->panel - 1 + 2) % 6;
//             eid = (ProcY-Proc->ypid-1) * ProcY + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], eid, recv_type, EAST, 1, 1, id+EAST, HaloP, HaloP+ny-1, 0, HaloX-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], eid, send_type, SOUTH, 1, 0, eid+SOUTH, nx, nx+HaloX-1, HaloY, HaloY+ny-1);

//             if (Proc->at_south){
//                 //右下角
//                 pid = eid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 1, id+EAST+NORTH, 0, HaloP-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+EAST,nx, nx+HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);

//             }
//             else if (Proc->at_north){
//                 //右上角
//                 pid = eid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 1, id+EAST+SOUTH, HaloP+ny, 2*HaloP+ny-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+WEST, nx, nx+HaloX-1, HaloY, HaloY+HaloP-1);

//             }
//             else{
//                 //普通边
//                 pid = eid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 1, id+EAST+NORTH, 0, HaloP-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+EAST,nx, nx+HaloX-1, HaloY+ny-HaloP, HaloY+ny-1);
//                 pid = eid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvE[Proc->NgbER++], pid, recv_type, EAST, 1, 1, id+EAST+SOUTH, HaloP+ny, 2*HaloP+ny-1, 0, HaloX-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+WEST, nx, nx+HaloX-1, HaloY, HaloY+HaloP-1);
//             }
//         }
//     }

//     //SOUTH_Ngb
//     if (!Proc->at_south){
//         //不在边界
//         sid = id - 1;
//         //西南
//         if (!Proc->at_west){
//             pid = sid - ProcY;
//             NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 0, 0, id+SOUTH+WEST, 0, HaloX-1, 0, HaloY-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 0, 0, pid+NORTH+EAST, HaloX, HaloX+HaloX-1, HaloY, 2*HaloY-1);
//         }
//         //正南
//         pid = sid;
//         NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 0, 0, id+SOUTH, HaloX, HaloX+nx-1, 0, HaloY-1);  
//         NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 0, 0, pid+NORTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);
//         //东南
//         if (!Proc->at_east){
//             pid = sid + ProcY;
//             NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 0, 0, id+SOUTH+EAST, HaloX+nx, 2*HaloX+nx-1, 0, HaloY-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 0, 0, pid+NORTH+WEST, nx, nx+HaloX-1, HaloY, 2*HaloY-1);
//         }
//     }
//     else{
//         //在边,分奇偶
//         if (Proc->panel & 1){
//             //奇
//             pn = (Proc->panel - 1 + 5) % 6;
//             sid = (Proc->xpid+1)*ProcY-1 + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, 1, 0, id+SOUTH, HaloP, HaloP+nx-1, 0, HaloY-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], sid, send_type, NORTH, 1, 0, sid+NORTH, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);

//             if (Proc->at_west){
//                 //左下角
//                 pid = sid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 0, id+SOUTH+EAST, HaloP+nx, 2*HaloP+nx-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+WEST, nx, nx+HaloX-1, HaloY, 2*HaloY-1);

//             }
//             else if (Proc->at_east){
//                 //右下角
//                 pid = sid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 0, id+SOUTH+WEST, 0, HaloP-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+EAST,HaloX, HaloX+HaloX-1, HaloY, 2*HaloY-1);
//             }
//             else{
//                 //普通边
//                 pid = sid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 0, id+SOUTH+EAST, HaloP+nx, 2*HaloP+nx-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+WEST, nx, nx+HaloX-1, HaloY, 2*HaloY-1);
//                 pid = sid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 0, id+SOUTH+WEST, 0, HaloP-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendN[Proc->NgbNS++], pid, send_type, NORTH, 1, 0, pid+NORTH+EAST,HaloX, HaloX+HaloX-1, HaloY, 2*HaloY-1);
//             }
//         }
//         else{
//             //偶
//             pn = (Proc->panel - 1 + 4) % 6;
//             sid = (ProcX-1)*ProcY + ProcY-Proc->xpid-1 + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], sid, recv_type, SOUTH, 1, 1, id+SOUTH, 0, HaloX-1, HaloP, HaloP+nx-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], sid, send_type, EAST, 1, 0, sid+EAST, HaloX, HaloX+nx-1, HaloY, 2*HaloY-1);

//             if (Proc->at_west){
//                 //左下角
//                 pid = sid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 1, id+SOUTH+EAST, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+NORTH, nx, nx+HaloX-1, HaloY, 2*HaloY-1);
                
//             }
//             else if (Proc->at_east){
//                 //右下角
//                 pid = sid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 1, id+SOUTH+WEST,0, HaloX-1, HaloP+nx, 2*HaloP+nx-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+SOUTH, HaloX, HaloX+HaloX-1, HaloY, 2*HaloY-1);
//             }
//             else{
//                 //普通边
//                 pid = sid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 1, id+SOUTH+EAST, 0, HaloX-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+NORTH, nx, nx+HaloX-1, HaloY, 2*HaloY-1);
//                 pid = sid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvS[Proc->NgbSR++], pid, recv_type, SOUTH, 1, 1, id+SOUTH+WEST,0, HaloX-1, HaloP+nx, 2*HaloP+nx-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendE[Proc->NgbES++], pid, send_type, EAST, 1, 0, pid+EAST+SOUTH, HaloX, HaloX+HaloX-1, HaloY, 2*HaloY-1);
//             }
//         }
//     }

//     //NORTH_Ngb
//     if (!Proc->at_north){
//         //不在边界
//         nid = id + 1;
//         //西北
//         if (!Proc->at_west){
//             pid = nid - ProcY;
//             NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 0, 0, id+NORTH+WEST, 0, HaloX-1, HaloY+ny, 2*HaloY+ny-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 0, 0, pid+SOUTH+EAST, HaloX, HaloX+HaloX-1, ny, ny+HaloY-1);
//         }
//         //正北
//         pid = nid;
//         NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 0, 0, id+NORTH, HaloX, HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1);
//         NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 0, 0, pid+SOUTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);
//         //东北
//         if (!Proc->at_east){
//             pid = nid + ProcY;
//             NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 0, 0, id+NORTH+EAST, HaloX+nx, 2*HaloX+nx-1, HaloY+ny, 2*HaloY+ny-1); 
//             NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 0, 0, pid+SOUTH+WEST, nx, HaloX+nx-1, ny, ny+HaloY-1);
//         } 
//     }
//     else{
//         //在边,分奇偶
//         if (Proc->panel & 1){
//             //奇
//             pn = (Proc->panel - 1 + 2) % 6;
//             nid = ProcY - Proc->xpid - 1 + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, 1, 1, id+NORTH, 0, HaloY-1, HaloP, HaloP+nx-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], nid, send_type, WEST, 1, 0, nid+WEST, HaloX, HaloX+nx-1, ny, ny+HaloY-1);

//             if (Proc->at_west){
//                 //左上角
//                 pid = nid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 1, id+NORTH+EAST, 0, HaloY-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+NORTH, nx, HaloX+nx-1, ny, ny+HaloY-1);
//             }
//             else if (Proc->at_east){
//                 //右上角
//                 pid = nid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 1, id+NORTH+WEST, 0, HaloY-1, HaloP+nx, 2*HaloP+nx-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+SOUTH, HaloX, HaloX+HaloX-1, ny, ny+HaloY-1);
//             }
//             else{
//                 //普通边
//                 pid = nid - 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 1, id+NORTH+EAST, 0, HaloY-1, 0, HaloP-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+NORTH, nx, HaloX+nx-1, ny, ny+HaloY-1);
//                 pid = nid + 1;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 1, id+NORTH+WEST, 0, HaloY-1, HaloP+nx, 2*HaloP+nx-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendW[Proc->NgbWS++], pid, send_type, WEST, 1, 0, pid+WEST+SOUTH, HaloX, HaloX+HaloX-1, ny, ny+HaloY-1);
//             }
//         }
//         else{
//             //偶
//             pn = (Proc->panel - 1 + 1) % 6;
//             nid = Proc->xpid * ProcY + pn*ProcX*ProcY;

//             NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], nid, recv_type, NORTH, 1, 0, id+NORTH, HaloP, HaloP+nx-1, 0, HaloY-1);
//             NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], nid, send_type, SOUTH, 1, 0, nid+SOUTH, HaloX, HaloX+nx-1, ny, ny+HaloY-1);

//             if (Proc->at_west){
//                 //左上角
//                 pid = nid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 0, id+NORTH+EAST, HaloP+nx, 2*HaloP+nx-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+WEST, nx, HaloX+nx-1, ny, ny+HaloY-1);
//             }
//             else if (Proc->at_east){
//                 //右上角
//                 pid = nid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 0, id+NORTH+WEST, 0, HaloP-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+EAST, HaloX, HaloX+HaloX-1, ny, ny+HaloY-1);
//             }
//             else{
//                 //普通边
//                 pid = nid + ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 0, id+NORTH+EAST, HaloP+nx, 2*HaloP+nx-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+WEST, nx, HaloX+nx-1, ny, ny+HaloY-1);
//                 pid = nid - ProcY;
//                 NgbTypeInitCS_Z(&Proc->NgbRecvN[Proc->NgbNR++], pid, recv_type, NORTH, 1, 0, id+NORTH+WEST, 0, HaloP-1, 0, HaloY-1);
//                 NgbTypeInitCS_Z(&Proc->NgbSendS[Proc->NgbSS++], pid, send_type, SOUTH, 1, 0, pid+SOUTH+EAST, HaloX, HaloX+HaloX-1, ny, ny+HaloY-1);
//             }
//         }
//     }

//     //Ngb转Halo
//     for (k = 0 ; k < size ; k++){
//         if (id == k){
//             //printf("=============================\n");
//             //printf("ID %d Panel %d Pid %d [%d %d %d %d]\n", Proc->id, Proc->panel, Proc->panelrank, Proc->at_west, Proc->at_east, Proc->at_south, Proc->at_north);
//             //Send
//             //printf("Send\n");
//             for (i = 0; i < Proc->NgbWS; i++){
//                 ngb = Proc->NgbSendW[i];
//                 HaloTypeInitCS(&Proc->SendHaloW[i], Proc, send_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                      ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
//                 //printf("West ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             } 
//             for (i = 0; i < Proc->NgbES; i++){
//                 ngb = Proc->NgbSendE[i];
//                 HaloTypeInitCS(&Proc->SendHaloE[i], Proc, send_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                      ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
//                 //printf("East ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//             for (i = 0; i < Proc->NgbSS; i++){
//                 ngb = Proc->NgbSendS[i];
//                 HaloTypeInitCS(&Proc->SendHaloS[i], Proc, send_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                      ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
//                 //printf("South ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//             for (i = 0; i < Proc->NgbNS; i++){
//                 ngb = Proc->NgbSendN[i];
//                 HaloTypeInitCS(&Proc->SendHaloN[i], Proc, send_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                      ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbSendNum);
//                 //printf("North ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
            
//             //Recv
//             //printf("Recv\n");
//             for (i = 0; i < Proc->NgbWR; i++){
//                 ngb = Proc->NgbRecvW[i];
//                 if (ngb.pos == 0)
//                     HaloTypeInitCS(&Proc->RecvHaloW[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                         ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
//                 else{
//                     if (ngb.rotate == 0)
//                         HaloTypeInitCS(&Proc->RecvHaloW[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, HaloX, ny+2*HaloP, NumLev,0, 0, ++Proc->NgbRecvNum);
//                     else
//                         HaloTypeInitCS(&Proc->RecvHaloW[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, ny+2*HaloP, HaloX, NumLev,0, 0, ++Proc->NgbRecvNum);
//                 }
//                 //printf("West ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//             for (i = 0; i < Proc->NgbER; i++){
//                 ngb = Proc->NgbRecvE[i];
//                 if (ngb.pos == 0)
//                     HaloTypeInitCS(&Proc->RecvHaloE[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                         ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
//                 else{
//                     if (ngb.rotate == 0)
//                         HaloTypeInitCS(&Proc->RecvHaloE[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, HaloX, ny+2*HaloP, NumLev,0, 0, ++Proc->NgbRecvNum);
//                     else
//                         HaloTypeInitCS(&Proc->RecvHaloE[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, ny+2*HaloP, HaloX, NumLev,0, 0, ++Proc->NgbRecvNum);
//                 }
//                 //printf("East ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//             for (i = 0; i < Proc->NgbSR; i++){
//                 ngb = Proc->NgbRecvS[i];
//                 if (ngb.pos == 0)
//                     HaloTypeInitCS(&Proc->RecvHaloS[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                         ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
//                 else{
//                     if (ngb.rotate == 0)
//                         HaloTypeInitCS(&Proc->RecvHaloS[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx+2*HaloP, HaloY, NumLev,0, 0, ++Proc->NgbRecvNum);
//                     else
//                         HaloTypeInitCS(&Proc->RecvHaloS[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, HaloY, nx+2*HaloP, NumLev,0, 0, ++Proc->NgbRecvNum);
//                 }
//                 //printf("South ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//             for (i = 0; i < Proc->NgbNR; i++){
//                 ngb = Proc->NgbRecvN[i];
//                 if (ngb.pos == 0)
//                     HaloTypeInitCS(&Proc->RecvHaloN[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                         ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx, ny, NumLev,HaloX, HaloY, ++Proc->NgbRecvNum);
//                 else{
//                     if (ngb.rotate == 0)
//                         HaloTypeInitCS(&Proc->RecvHaloN[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, nx+2*HaloP, HaloY, NumLev,0, 0, ++Proc->NgbRecvNum);
//                     else
//                         HaloTypeInitCS(&Proc->RecvHaloN[i], Proc, recv_type, ngb.orient, dtype, id, ngb.id, ngb.tag, \
//                             ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end, HaloY, nx+2*HaloP, NumLev,0, 0, ++Proc->NgbRecvNum);
//                 }
//                 //printf("North ngb = %d tag %d [%d %d %d %d]\n",ngb.id,ngb.tag,ngb.lon_beg, ngb.lon_end, ngb.lat_beg, ngb.lat_end);
//             }
//         }
//         MPI_Barrier(comm);
//     }

//     //预处理async_req
//     Proc->FieldReq = (AsyncType*)malloc(FieldNum * sizeof(AsyncType));

//     for (i = 0 ; i < FieldNum ; i++){
//         for (p = 1 ; p <= Proc->NgbSendNum ; p++) 
//             Proc->FieldReq[i].send_req[p] = MPI_REQUEST_NULL;
//         for (p = 1 ; p <= Proc->NgbRecvNum ; p++) 
//             Proc->FieldReq[i].recv_req[p] = MPI_REQUEST_NULL;
//     }

// }