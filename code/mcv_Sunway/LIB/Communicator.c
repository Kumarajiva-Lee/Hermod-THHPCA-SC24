#include<stdlib.h>
#include<stdbool.h>
#include"mpi.h"

#include"Process.h"
#include"ParaParam.h"
#include"Communicator.h"
#include"Memory.h"

void swap(double *a, double *b){
  double t;
  t = *a;
  *a = *b;
  *b = t;
}

void UpdateHalo_2d_I(ProcType Proc, int *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,p;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;

    if (west_halo){
        if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
        if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_2d[j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
        MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_2d[j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
    }

    if (east_halo){
        if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
        if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

        MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_2d[j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
        MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_2d[j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    }

    if (south_halo){
        if (!Proc.at_south_pole)
            for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }   

        if (!Proc.at_north_pole)
            for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (north_halo ){
        if (!Proc.at_north_pole)
            for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }
        if (!Proc.at_south_pole)
            for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (Proc.at_north_pole || Proc.at_south_pole){
        int nx,mx,ny,nz,hx,hy,hz;
        hx = Proc.lon_hw;
        hy = Proc.lat_hw;
        hz = Proc.lev_hw;
        nx = (i == 0)?Proc.full_nlon:Proc.half_nlon;
        nx += 2 * hx;
        mx = (nx + 1) / 2;
        ny = (j == 0)?Proc.full_nlat:Proc.half_nlat;
        ny += 2 * hy;
        nz = 1;
        nz += 2 * hz;

        double tmp[nz][ny][nx];

        if (south_halo && Proc.at_south_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_2d[j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_2d[j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = 0; jj < hy; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii]; 
            } 
        }
        if (north_halo && Proc.at_north_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_2d[j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_2d[j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = ny - hy; jj < ny; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj - ny + hy][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii];
            }
        }

    }
}

void UpdateHalo_3d_I(ProcType Proc, int *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,k,p;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;
    k = (full_lev)?0:1;

    if (west_halo){
        if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
        if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_3d[k][j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
        MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_3d[k][j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
    }

    if (east_halo){
        if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
        if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

        MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
        MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_3d[k][j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    }

    if (south_halo){
        if (!Proc.at_south_pole)
            for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }   
        if (!Proc.at_north_pole)    
            for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (north_halo){
        if (!Proc.at_north_pole)
            for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }
        if (!Proc.at_south_pole)
            for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (Proc.at_north_pole || Proc.at_south_pole){
        int nx,mx,ny,nz,hx,hy,hz;
        hx = Proc.lon_hw;
        hy = Proc.lat_hw;
        hz = Proc.lev_hw;
        nx = (i == 0)?Proc.full_nlon:Proc.half_nlon;
        nx += 2 * hx;
        mx = (nx + 1) / 2;
        ny = (j == 0)?Proc.full_nlat:Proc.half_nlat;
        ny += 2 * hy;
        nz = (k == 0)?Proc.full_nlev:Proc.half_nlev;
        nz += 2 * hz;

        double tmp[nz][ny][nx];

        if (south_halo && Proc.at_south_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = 0; jj < hy; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii]; 
            } 
        }
        if (north_halo && Proc.at_north_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_2d[j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_2d[j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = ny - hy; jj < ny; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj - ny + hy][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii];
            }
        }

    }
}

void UpdateHalo_2d_S(ProcType Proc, float *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,p;
    int ierr;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;

    if (west_halo){
        if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
        if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_2d[j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
        MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_2d[j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
    }

    if (east_halo){
        if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
        if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

        MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_2d[j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
        MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_2d[j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    }

    if (south_halo){
        for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
            if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
            MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
        }       
        for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
            if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
            MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
        }
    }

    if (north_halo){
        for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
            if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
            MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
        }
        for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
            if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
            MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
        }
    }
}

void UpdateHalo_3d_S(ProcType Proc, float *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,k,p;
    int ierr;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;
    k = (full_lev)?0:1;

    if (west_halo){
        if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
        if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_3d[k][j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
        MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_3d[k][j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
    }

    if (east_halo){
        if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
        if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

        MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
        MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_3d[k][j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    }

    if (south_halo){
        for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
            if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
            MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
        }       
        for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
            if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
            MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
        }
    }

    if (north_halo){
        for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
            if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
            MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
        }
        for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
            if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
            MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
        }
    }
}

//ToRemember: Add west_s to other
void UpdateHalo_2d_D(ProcType Proc, double *array, AsyncType *async, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo, bool small){
    int i,j,p;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;

    if (small){
      if (west_halo){
          if (async->recv_req[WEST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST_S], MPI_STATUSES_IGNORE);
          if (async->send_req[WEST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST_S], MPI_STATUSES_IGNORE);
          MPI_Irecv(array, 1, Proc.RecvHalo[WEST_S].mpi_type_2d[j][i], Proc.RecvHalo[WEST_S].ngb_id, Proc.RecvHalo[WEST_S].tag, Proc.comm, &async->recv_req[WEST_S]);
          MPI_Isend(array, 1, Proc.SendHalo[WEST_S].mpi_type_2d[j][i], Proc.SendHalo[WEST_S].ngb_id, Proc.SendHalo[WEST_S].tag, Proc.comm, &async->send_req[WEST_S]); 
      }

      if (east_halo){
          if (async->recv_req[EAST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST_S], MPI_STATUSES_IGNORE);
          if (async->send_req[EAST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST_S], MPI_STATUSES_IGNORE);

          MPI_Irecv(array, 1, Proc.RecvHalo[EAST_S].mpi_type_2d[j][i], Proc.RecvHalo[EAST_S].ngb_id, Proc.RecvHalo[EAST_S].tag, Proc.comm, &async->recv_req[EAST_S]);
          MPI_Isend(array, 1, Proc.SendHalo[EAST_S].mpi_type_2d[j][i], Proc.SendHalo[EAST_S].ngb_id, Proc.SendHalo[EAST_S].tag, Proc.comm, &async->send_req[EAST_S]);
      }
    }
    else{
      if (west_halo){
          if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
          if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
          MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_2d[j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
          MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_2d[j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
      }

      if (east_halo){
          if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
          if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

          MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_2d[j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
          MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_2d[j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
      }
    }

    if (south_halo){
        if (!Proc.at_south_pole)
            for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }       
        if (!Proc.at_north_pole)
            for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (north_halo){
        if (!Proc.at_north_pole)
            for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_2d[j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }
        if (!Proc.at_south_pole)
            for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_2d[j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if ((Proc.at_north_pole || Proc.at_south_pole) && Proc.full_nlat >=3){
        int nx,mx,ny,nz,hx,hy,hz;
        hx = Proc.lon_hw;
        hy = Proc.lat_hw;
        hz = Proc.lev_hw;
        nx = (i == 0)?Proc.full_nlon:Proc.half_nlon;
        nx += 2 * hx;
        mx = (nx + 1) / 2;
        ny = (j == 0)?Proc.full_nlat:Proc.half_nlat;
        ny += 2 * hy;
        nz = 1;
        //nz += 2 * hz;

        double tmp[nz][ny][nx];

        if (south_halo && Proc.at_south_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_2d[j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_2d[j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = 0; jj < hy; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii]; 
            } 
        }
        if (north_halo && Proc.at_north_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_2d[j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_2d[j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = ny - hy; jj < ny; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj - ny + hy][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii];
            }
        }

    }
}

void UpdateHalo_3d_D(ProcType Proc, double *array, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo, bool small){
    int i,j,k,p;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;
    int count;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;
    k = (full_lev)?0:1;

    if (small){
      if (west_halo){
          if (async->recv_req[WEST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST_S], MPI_STATUSES_IGNORE);
          if (async->send_req[WEST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST_S], MPI_STATUSES_IGNORE);
          MPI_Irecv(array, 1, Proc.RecvHalo[WEST_S].mpi_type_3d[k][j][i], Proc.RecvHalo[WEST_S].ngb_id, Proc.RecvHalo[WEST_S].tag, Proc.comm, &async->recv_req[WEST_S]);
          MPI_Isend(array, 1, Proc.SendHalo[WEST_S].mpi_type_3d[k][j][i], Proc.SendHalo[WEST_S].ngb_id, Proc.SendHalo[WEST_S].tag, Proc.comm, &async->send_req[WEST_S]); 
      }

      if (east_halo){
          if (async->recv_req[EAST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST_S], MPI_STATUSES_IGNORE);
          if (async->send_req[EAST_S] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST_S], MPI_STATUSES_IGNORE);

          MPI_Irecv(array, 1, Proc.RecvHalo[EAST_S].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST_S].ngb_id, Proc.RecvHalo[EAST_S].tag, Proc.comm, &async->recv_req[EAST_S]);
          MPI_Isend(array, 1, Proc.SendHalo[EAST_S].mpi_type_3d[k][j][i], Proc.SendHalo[EAST_S].ngb_id, Proc.SendHalo[EAST_S].tag, Proc.comm, &async->send_req[EAST_S]);
      }
    }
    else{
      if (west_halo){
          if (async->recv_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[WEST], MPI_STATUSES_IGNORE);
          if (async->send_req[WEST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[WEST], MPI_STATUSES_IGNORE);
          MPI_Irecv(array, 1, Proc.RecvHalo[WEST].mpi_type_3d[k][j][i], Proc.RecvHalo[WEST].ngb_id, Proc.RecvHalo[WEST].tag, Proc.comm, &async->recv_req[WEST]);
          MPI_Isend(array, 1, Proc.SendHalo[WEST].mpi_type_3d[k][j][i], Proc.SendHalo[WEST].ngb_id, Proc.SendHalo[WEST].tag, Proc.comm, &async->send_req[WEST]); 
      }

      if (east_halo){
          if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
          if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

          MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
          MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_3d[k][j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
      }
    }

    if (south_halo){
        if (!Proc.at_south_pole)
            for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }     
        if (!Proc.at_north_pole)  
            for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if (north_halo){
        if (!Proc.at_north_pole)
            for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
                if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
                MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
            }
        if (!Proc.at_south_pole)
            for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
                if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
                MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
            }
    }

    if ((Proc.at_north_pole || Proc.at_south_pole) && Proc.full_nlat >=3){
        int nx,mx,ny,nz,hx,hy,hz;
        hx = Proc.lon_hw;
        hy = Proc.lat_hw;
        hz = Proc.lev_hw;
        nx = (i == 0)?Proc.full_nlon:Proc.half_nlon;
        nx += 2 * hx;
        mx = (nx + 1) / 2;
        ny = (j == 0)?Proc.full_nlat:Proc.half_nlat;
        ny += 2 * hy;
        nz = (k == 0)?Proc.full_nlev:Proc.half_nlev;
        nz += 2 * hz;

        double tmp[nz][ny][nx];

        if (south_halo && Proc.at_south_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = 0; jj < hy; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (hy - 1 - jj) * nx + ii) = tmp[kk][jj][ii]; 
            } 
        }
        if (north_halo && Proc.at_north_pole){
            MPI_Sendrecv(array, 1, Proc.SendHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.SendHalo[OPPOSITE].ngb_id, 0,
                         array, 1, Proc.RecvHalo[OPPOSITE].mpi_type_3d[k][j][i], Proc.RecvHalo[OPPOSITE].ngb_id, 0,
                         Proc.comm, &status);
            for (kk = 0 ; kk < nz; kk++)
              for (jj = ny - hy; jj < ny; jj++)
                for (ii = 0 ; ii < nx ; ii++)
                  tmp[kk][jj - ny + hy][ii] = *(array + kk * nx * ny + jj * nx + ii);
            if (Proc.SendHalo[OPPOSITE].ngb_id == Proc.id){
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = hx; ii < mx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii + mx - hx];
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = mx; ii < nx - hx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii - mx + hx];
            }
            else{
              for (kk = 0 ; kk < nz; kk++)
                for (jj = 0; jj < hy; jj++)
                  for (ii = 0; ii < nx; ii++)
                    *(array + kk * nx * ny + (ny - jj - 1) * nx + ii) = tmp[kk][jj][ii];
            }
        }

    }
}

void UpdateHaloCS_2d_D(ProcType Proc, double *array, AsyncType *async, int varid, bool full_lon, bool full_lat, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,k,p,req;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;
    int count;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;

    Proc.CSWait[varid] = 1;

    if (west_halo){
        for (p = 0; p < Proc.NgbWS; p++){
          req = Proc.SendHaloW[p].req;

          if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
          MPI_Isend(array, 1, Proc.SendHaloW[p].mpi_type_2d[j][i], Proc.SendHaloW[p].ngb_id, Proc.SendHaloW[p].tag, Proc.comm, &async->send_req[req]); 
        }
        for (p = 0; p < Proc.NgbWR; p++){
          req = Proc.RecvHaloW[p].req;
          if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
          MPI_Irecv(array, 1, Proc.RecvHaloW[p].mpi_type_2d[j][i], Proc.RecvHaloW[p].ngb_id, Proc.RecvHaloW[p].tag, Proc.comm, &async->recv_req[req]);
        }
    }

    if (east_halo){
      for (p = 0; p < Proc.NgbES; p++){
        req = Proc.SendHaloE[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloE[p].mpi_type_2d[j][i], Proc.SendHaloE[p].ngb_id, Proc.SendHaloE[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbER; p++){
        req = Proc.RecvHaloE[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHaloE[p].mpi_type_2d[j][i], Proc.RecvHaloE[p].ngb_id, Proc.RecvHaloE[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    if (south_halo){
      for (p = 0; p < Proc.NgbSS; p++){
        req = Proc.SendHaloS[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloS[p].mpi_type_2d[j][i], Proc.SendHaloS[p].ngb_id, Proc.SendHaloS[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbSR; p++){
        req = Proc.RecvHaloS[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHaloS[p].mpi_type_2d[j][i], Proc.RecvHaloS[p].ngb_id, Proc.RecvHaloS[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    if (north_halo){
      for (p = 0; p < Proc.NgbNS; p++){
        req = Proc.SendHaloN[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloN[p].mpi_type_2d[j][i], Proc.SendHaloN[p].ngb_id, Proc.SendHaloN[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbNR; p++){
        req = Proc.RecvHaloN[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        MPI_Irecv(array, 1, Proc.RecvHaloN[p].mpi_type_2d[j][i], Proc.RecvHaloN[p].ngb_id, Proc.RecvHaloN[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    // if (east_halo){
    //     if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
    //     if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

    //     MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
    //     MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_3d[k][j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    // }

    // if (south_halo){
    //     if (!Proc.at_south_pole)
    //         for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
    //             if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
    //         }     
    //     if (!Proc.at_north_pole)  
    //         for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
    //             if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
    //         }
    // }

    // if (north_halo){
    //     if (!Proc.at_north_pole)
    //         for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
    //             if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
    //         }
    //     if (!Proc.at_south_pole)
    //         for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
    //             if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
    //         }
    // }    
}

//arrayX arrayY代表x y方向上接收外Halo的临时数组
void UpdateHaloCS_3d_D(ProcType Proc, double *array, double *arrayX, double *arrayY, AsyncType *async, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int i,j,k,p,req;
    int ierr;
    MPI_Status status;
    int ii,jj,kk;
    int count;

    i = (full_lon)?0:1;
    j = (full_lat)?0:1;
    k = (full_lev)?0:1;



    if (west_halo){
        for (p = 0; p < Proc.NgbWS; p++){
          req = Proc.SendHaloW[p].req;
          // printf("%d %d\n",Proc.id,req);
          // if (req >= ngb_max) printf("ERROR!!!!!\n");
          if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
          MPI_Isend(array, 1, Proc.SendHaloW[p].mpi_type_3d[k][j][i], Proc.SendHaloW[p].ngb_id, Proc.SendHaloW[p].tag, Proc.comm, &async->send_req[req]); 
        }
        for (p = 0; p < Proc.NgbWR; p++){
          req = Proc.RecvHaloW[p].req;
          if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
          if (Proc.at_west)
            MPI_Irecv(arrayX, 1, Proc.RecvHaloW[p].mpi_type_3d[k][j][i], Proc.RecvHaloW[p].ngb_id, Proc.RecvHaloW[p].tag, Proc.comm, &async->recv_req[req]);
          else
            MPI_Irecv(array, 1, Proc.RecvHaloW[p].mpi_type_3d[k][j][i], Proc.RecvHaloW[p].ngb_id, Proc.RecvHaloW[p].tag, Proc.comm, &async->recv_req[req]);
        }
    }

    if (east_halo){
      for (p = 0; p < Proc.NgbES; p++){
        req = Proc.SendHaloE[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloE[p].mpi_type_3d[k][j][i], Proc.SendHaloE[p].ngb_id, Proc.SendHaloE[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbER; p++){
        req = Proc.RecvHaloE[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        if (Proc.at_east)
          MPI_Irecv(arrayX, 1, Proc.RecvHaloE[p].mpi_type_3d[k][j][i], Proc.RecvHaloE[p].ngb_id, Proc.RecvHaloE[p].tag, Proc.comm, &async->recv_req[req]);
        else
          MPI_Irecv(array, 1, Proc.RecvHaloE[p].mpi_type_3d[k][j][i], Proc.RecvHaloE[p].ngb_id, Proc.RecvHaloE[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    if (south_halo){
      for (p = 0; p < Proc.NgbSS; p++){
        req = Proc.SendHaloS[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloS[p].mpi_type_3d[k][j][i], Proc.SendHaloS[p].ngb_id, Proc.SendHaloS[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbSR; p++){
        req = Proc.RecvHaloS[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        if (Proc.at_south)
          MPI_Irecv(arrayY, 1, Proc.RecvHaloS[p].mpi_type_3d[k][j][i], Proc.RecvHaloS[p].ngb_id, Proc.RecvHaloS[p].tag, Proc.comm, &async->recv_req[req]);
        else
          MPI_Irecv(array, 1, Proc.RecvHaloS[p].mpi_type_3d[k][j][i], Proc.RecvHaloS[p].ngb_id, Proc.RecvHaloS[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    if (south_halo){
      for (p = 0; p < Proc.NgbNS; p++){
        req = Proc.SendHaloN[p].req;
        if (async->send_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[req], MPI_STATUSES_IGNORE);
        MPI_Isend(array, 1, Proc.SendHaloN[p].mpi_type_3d[k][j][i], Proc.SendHaloN[p].ngb_id, Proc.SendHaloN[p].tag, Proc.comm, &async->send_req[req]); 
      }
      for (p = 0; p < Proc.NgbNR; p++){
        req = Proc.RecvHaloN[p].req;
        if (async->recv_req[req] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[req], MPI_STATUSES_IGNORE);
        if (Proc.at_north)
          MPI_Irecv(arrayY, 1, Proc.RecvHaloN[p].mpi_type_3d[k][j][i], Proc.RecvHaloN[p].ngb_id, Proc.RecvHaloN[p].tag, Proc.comm, &async->recv_req[req]);
        else
          MPI_Irecv(array, 1, Proc.RecvHaloN[p].mpi_type_3d[k][j][i], Proc.RecvHaloN[p].ngb_id, Proc.RecvHaloN[p].tag, Proc.comm, &async->recv_req[req]);
      }
    }

    // if (east_halo){
    //     if (async->recv_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[EAST], MPI_STATUSES_IGNORE);
    //     if (async->send_req[EAST] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[EAST], MPI_STATUSES_IGNORE);

    //     MPI_Irecv(array, 1, Proc.RecvHalo[EAST].mpi_type_3d[k][j][i], Proc.RecvHalo[EAST].ngb_id, Proc.RecvHalo[EAST].tag, Proc.comm, &async->recv_req[EAST]);
    //     MPI_Isend(array, 1, Proc.SendHalo[EAST].mpi_type_3d[k][j][i], Proc.SendHalo[EAST].ngb_id, Proc.SendHalo[EAST].tag, Proc.comm, &async->send_req[EAST]);
    // }

    // if (south_halo){
    //     if (!Proc.at_south_pole)
    //         for (p = Proc.NgbRSP; p < Proc.NgbRNP; p++){
    //             if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
    //         }     
    //     if (!Proc.at_north_pole)  
    //         for (p = Proc.NgbSSP; p <= Proc.NgbSendNum; p++){
    //             if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
    //         }
    // }

    // if (north_halo){
    //     if (!Proc.at_north_pole)
    //         for (p = Proc.NgbRNP; p <= Proc.NgbRecvNum; p++){
    //             if (async->recv_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->recv_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Irecv(array, 1, Proc.RecvHalo[p].mpi_type_3d[k][j][i], Proc.RecvHalo[p].ngb_id, Proc.RecvHalo[p].tag, Proc.comm, &async->recv_req[p]);
    //         }
    //     if (!Proc.at_south_pole)
    //         for (p = Proc.NgbSNP; p < Proc.NgbSSP; p++){
    //             if (async->send_req[p] != MPI_REQUEST_NULL) MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    //             MPI_Isend(array, 1, Proc.SendHalo[p].mpi_type_3d[k][j][i], Proc.SendHalo[p].ngb_id, Proc.SendHalo[p].tag, Proc.comm, &async->send_req[p]);
    //         }
    // }    
}

void HaloWait(ProcType Proc, AsyncType *async){
    int p;
    MPI_Status status;
    for (p = 1; p <= Proc.NgbSendNum; p++)
        MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    for (p = 1; p <= Proc.NgbRecvNum; p++)
        MPI_Wait(&async->recv_req[p], &status);

}

void HaloWaitCS(ProcType Proc, AsyncType *async, int varid, double ***array, bool full_lon, bool full_lat, bool full_lev, bool west_halo, bool east_halo, bool south_halo, bool north_halo){
    int p;
    int lon_beg,lon_end,lat_beg,lat_end;
    int nlon,nlat;
    int i,j;
    MPI_Status status;
    for (p = 1; p <= Proc.NgbSendNum; p++)
        MPI_Wait(&async->send_req[p], MPI_STATUSES_IGNORE);
    for (p = 1; p <= Proc.NgbRecvNum; p++)
        MPI_Wait(&async->recv_req[p], &status);

    //Todo:将各种维度的halo区数据翻转,目前只支持3d空间的1维halo
    if (Proc.CSWait[varid] == 1){
      if (west_halo && Proc.re_west && Proc.at_west){
        if (full_lat) nlat = Proc.full_nlat;
        else nlat = Proc.half_nlat;
        for (j = Proc.lat_hw; j <= (nlat + Proc.lat_hw*2-1) / 2; j++) swap(&array[0][j][0],&array[0][nlat + Proc.lat_hw*2-1-j][0]);
      }
      if (east_halo && Proc.re_east && Proc.at_east){
        if (full_lon) nlon = Proc.full_nlon;
        else nlon = Proc.half_nlon;
        if (full_lat) nlat = Proc.full_nlat;
        else nlat = Proc.half_nlat;
        for (j = Proc.lat_hw; j <= (nlat + Proc.lat_hw*2-1) / 2; j++) swap(&array[0][j][nlon+Proc.lon_hw*2-1],&array[0][nlat + Proc.lat_hw*2-1-j][nlon+Proc.lon_hw*2-1]);
      }
      if (south_halo && Proc.re_south && Proc.at_south){
        if (full_lon) nlon = Proc.full_nlon;
        else nlon = Proc.half_nlon;
        for (i = Proc.lon_hw; i <= (nlon + Proc.lon_hw*2-1) / 2; i++) swap(&array[0][0][i],&array[0][0][nlon + Proc.lon_hw*2-1-i]);
      }
      if (north_halo && Proc.re_north && Proc.at_north){
        if (full_lon) nlon = Proc.full_nlon;
        else nlon = Proc.half_nlon;
        if (full_lat) nlat = Proc.full_nlat;
        else nlat = Proc.half_nlat;
        for (i = Proc.lon_hw; i <= (nlon + Proc.lon_hw*2-1) / 2; i++) swap(&array[0][nlat + Proc.lat_hw*2-1][i],&array[0][nlat + Proc.lat_hw*2-1][nlon + Proc.lon_hw*2-1-i]);
      }
    }

    Proc.CSWait[varid] = 0;
}

void Zonal_Sum_1d_D(ProcType Proc, double** sum,int size){
    int ierr;
    double *res = allocate_1d_array_D(size);
    MPI_Allreduce(*sum, res, size, MPI_DOUBLE, MPI_SUM, Proc.zonal_comm);
    *sum = res;

}
