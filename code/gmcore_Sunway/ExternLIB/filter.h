#ifndef FILTER_H_INCLUDED
#define FILTER_H_INCLUDED 1

#include"../LIB/Process.h"
#include"../src/physical_variable.h"

typedef struct{
    double *width_lon;
    int *ngrid_lon;
    double **wgt_lon;
    double *width_lat;
    int *ngrid_lat;
    double **wgt_lat;
}FilterType; 

// FilterType BigFilter;
FilterType Filter[4]; //0 = Big; 1 = phs; 2 = pt; 3 = uv



void Filter_Init(struct HybridMeshField* mesh, double dt);
void Filter_On_Cell(int id, bool is3d, double *x);
void Filter_on_lon_edge(int id, double *x);
void Filter_on_lat_edge(int id, double *x);
#endif
