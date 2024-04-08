#ifndef ASYNC_H_INCLUDED
#define ASYNC_H_INCLUDED

#include"mpi.h"
#include"ParaParam.h"

typedef struct{
    MPI_Request send_req[ngb_max];
    MPI_Request recv_req[ngb_max];
}AsyncType;

#endif