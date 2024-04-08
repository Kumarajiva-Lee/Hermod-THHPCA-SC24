from frontend.pytoy.lang import space_op, Float64
from frontend.pytoy.lang.dtype_ext import LonLatField


import gmcore_smallball.mesh as gm
import gmcore_smallball.global_config as gc
import gmcore_smallball.mesh as gm

#简化版
@space_op
def interp_lev_edge_to_lev_lon_edge(mesh:gm.HybridMeshField,x_lev:LonLatField[Float64],x_lev_lon:LonLatField[Float64]):
    # for k in range(gm.half_kds,gm.half_kde+1):
    #     for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
    #         for i in range(gm.half_ids,gm.half_ide+1):
    #             x_lev_lon[i,j,k] = mf.upwind3(mf.sign(1.0, u[i,j,k]), beta, x_lev[i-1:i+2,j,k])
    
    for k in range(gm.half_kds,gm.half_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                x_lev_lon[i,j,k] = (mesh.area_lon_west[0,j,0] * x_lev[i,j,k] + mesh.area_lon_east[0,j,0] * x_lev[i+1,j,k]) / mesh.area_lon[0,j,0]

@space_op
def interp_lev_edge_to_lev_lat_edge(mesh:gm.HybridMeshField,x_lev:LonLatField[Float64],x_lev_lat:LonLatField[Float64]):
    for k in range(gm.half_kds,gm.half_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                x_lev_lat[i,j,k] = (mesh.area_lat_north[0,j,0] * x_lev[i,j+1,k] + mesh.area_lat_south[0,j,0] * x_lev[i,j,k]) / mesh.area_lat[0,j,0]

@space_op
def averageCellToLonEdge(x:LonLatField[Float64], x_lon:LonLatField[Float64]):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                x_lon[i,j,k] = (x[i,j,k] + x[i+1,j,k]) * 0.5

@space_op
def averageCellToLatEdge(x:LonLatField[Float64], x_lat:LonLatField[Float64]):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                x_lat[i,j,k] = (x[i,j,k] + x[i,j+1,k]) * 0.5

@space_op
def interpCellToVtx(x:LonLatField[Float64], x_vtx:LonLatField[Float64], mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds, gm.half_jde+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                x_vtx[i,j,k] = (
                    (x[i,j,k] + x[i+1,j,k])   * mesh.area_subcell_1[0,j,0] + 
                    (x[i,j+1,k]+x[i+1,j+1,k]) * mesh.area_subcell_0[0,j+1,0] 
                ) / mesh.area_vtx[0,j,0]

# @space_op
# def levEdgeToLevLonEdge(x_lev:LonLatField[Float64], x_lev_lon:LonLatField[Float64], mesh:gm.HybridMeshField):
#     for k in range(0,gm.num_half_lev):
#         for j in range(gm.full_lat_ibeg_no_pole,gm.full_lat_ibeg_no_pole+1):
#             for i in range(0,gm.num_half_lon):
#                 x_lev_lon[i,j,k] = (mesh.area_lon_west[j] * x_lev[i,j,k] + 
#                                     mesh.area_lon_east[j] * x_lev[i+1,j,k]) / mesh.area_lon[j]

# @space_op
# def levEdgeToLevLatEdge(x_lev:LonLatField[Float64], x_lev_lat:LonLatField[Float64], mesh:gm.HybridMeshField):
#     for k in range(0,gm.num_half_lev):
#         for j in range(gm.half_lat_ibeg_no_pole,gm.half_lat_ibeg_no_pole+1):
#             for i in range(0,gm.num_full_lon):
#                 x_lev_lat[i,j,k] = (mesh.area_lon_north[j] * x_lev[i,j+1,k] + 
#                                     mesh.area_lon_south[j] * x_lev[i, j,k]) / mesh.area_lat[j]

