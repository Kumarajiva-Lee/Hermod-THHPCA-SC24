#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include"mpi.h"

#include "Process.h"
#include "Diagnose.h"
#include "Time.h"
#include "../src/physical_variable.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void Diagnose(struct HybridStateField *state, int xstep){
    // double vmin,vmax,vmean,vvar,sum;
    // double send[2],recv[2]; //sum tot
    // int lx,ly,lz;
    // int i,j,k;
    // int ids,ide,jds,jde,kds,kde,tot;

    // if (Timeinfo.timestep % xstep != 0) return;


    // kds = Proc.lev_hw;
    // kde = 32+Proc.lev_hw;

    // ids = Proc.lon_hw;
    // ide = Proc.full_nlon+Proc.lon_hw;
    
    // //u

    // jds = Proc.lat_hw;
    // jde = Proc.full_nlat+Proc.lat_hw;
    // if (Proc.id == 0) jds ++;
    // if (Proc.id == 19) jde --;

    // vmin = 10000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->u_lon[k][j][i] < vmin){
    //         vmin = state->u_lon[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }
    
    // send[0] = vmin;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MIN,Proc.comm);
    // vmin = recv[0];

    // vmax = -10000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->u_lon[k][j][i] > vmax){
    //         vmax = state->u_lon[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }
    
    // send[0] = vmax;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MAX,Proc.comm);
    // vmax = recv[0];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += state->u_lon[k][j][i];
    //       tot++;
    //     }

    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vmean = recv[0] / recv[1];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += pow(vmean - state->u_lon[k][j][i],2);
    //       tot++;
    //     }
    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vvar = recv[0] / recv[1];

    // if (Proc.id == 0) printf("U %.16f %.16f %.16f %.16f\n",vmin,vmax,vmean,vvar);

    // //v
    // jds = Proc.lat_hw;
    // jde = Proc.full_nlat+Proc.lat_hw;
    // if (Proc.id == 19) jde --;

    //   //min
    // vmin = 10000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->v_lat[k][j][i] < vmin){
    //         vmin = state->v_lat[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }

    // send[0] = vmin;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MIN,Proc.comm);
    // vmin = recv[0];

    // vmax = -10000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->v_lat[k][j][i] > vmax){
    //         vmax = state->v_lat[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }

    // send[0] = vmax;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MAX,Proc.comm);
    // vmax = recv[0];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += state->v_lat[k][j][i];
    //       tot++;
    //     }
    
    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vmean = recv[0] / recv[1];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += pow(vmean - state->v_lat[k][j][i],2);
    //       tot++;
    //     }
    
    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vvar = recv[0] / recv[1];

    // if (Proc.id == 0) printf("V %.16f %.16f %.16f %.16f\n",vmin,vmax,vmean,vvar);

    // //pt
    // jds = Proc.lat_hw;
    // jde = Proc.full_nlat+Proc.lat_hw;
    //   //min
    // vmin = 100000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->pt[k][j][i] < vmin){
    //         vmin = state->pt[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }

    // send[0] = vmin;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MIN,Proc.comm);
    // vmin = recv[0];

    // vmax = -100000;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++)
    //       if (state->pt[k][j][i] > vmax){
    //         vmax = state->pt[k][j][i];
    //         lx = i - Proc.lon_hw + Proc.lon_beg;
    //         ly = j - Proc.lat_hw + Proc.lat_beg;
    //         lz = k - Proc.lev_hw;
    //       }

    // send[0] = vmax;
    // MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_MAX,Proc.comm);
    // vmax = recv[0];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += state->pt[k][j][i];
    //       tot++;
    //     }
    
    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vmean = recv[0] / recv[1];

    // sum = 0;
    // tot = 0;
    // for (k = kds; k < kde ; k++)
    //   for (j = jds; j < jde ; j++)
    //     for (i = ids; i< ide ; i++){
    //       sum += pow(vmean - state->pt[k][j][i],2);
    //       tot++;
    //     }
    
    // send[0] = sum;
    // send[1] = tot;
    // MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,Proc.comm);
    // vvar = recv[0] / recv[1];

    // if (Proc.id == 0) printf("PT %.16f %.16f %.16f %.16f\n",vmin,vmax,vmean,vvar);
    
}