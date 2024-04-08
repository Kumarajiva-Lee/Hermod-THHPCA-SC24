import math

from frontend.pytoy.lang import space_op, Float64, Int32, extern_func
from frontend.pytoy.lang.dtype_ext import LonLatField

import gmcore_smallball.global_config as gc
import gmcore_smallball.physical_variable as phv
import gmcore_smallball.mesh as gm
import gmcore_smallball.const_value as cv
import gmcore_smallball.sphere_geometry as sg
import gmcore_smallball.interp as interp
import gmcore_smallball.math_formual as mf
import gmcore_smallball.ffsl as fs
import gmcore_smallball.filter as filter
from gmcore_smallball.utils import Vector4


@space_op
def prepareStatic(staticv:phv.HybridStaticField,mesh:gm.HybridMeshField):
    for j in range(gm.full_jds_no_pole, gm.full_jde_no_pole+1):
        for i in range(gm.half_ids,gm.half_ide+1):
            staticv.dzsdlon[i,j,0] = (staticv.gzs[i+1,j,0] - staticv.gzs[i,j,0]) / cv.g / mesh.de_lon[0,j,0]
    
    for j in range(gm.half_jds,gm.half_jde+1):
        for i in range(gm.full_ids,gm.full_ide+1):
            staticv.dzsdlat[i,j,0] = (staticv.gzs[i,j+1,0] - staticv.gzs[i,j,0]) / cv.g / mesh.de_lat[0,j,0]

@space_op
def preparegzlev(state:phv.HybridStateField,staticv:phv.HybridStaticField):
    for k in range(gm.half_kde,gm.half_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.gz_lev[i,j,k] = staticv.gzs[i,j,0]

@space_op
def c2a(state:phv.HybridStateField):
    for k in range(gm.full_kds, gm.full_kde+1):
        for j in range(gm.full_jds_no_pole, gm.full_jde_no_pole+1):
            for i in range(gm.full_ids, gm.full_ide):
                state.u[i,j,k] = 0.5 * (state.u_lon[i,j,k] + state.u_lon[i-1,j,k])
                state.v[i,j,k] = 0.5 * (state.v_lat[i,j,k] + state.v_lat[i,j-1,k])

@space_op 
def calc_ph(state:phv.HybridStateField,mesh:gm.HybridMeshField):

    for k in range(gm.half_kds,gm.half_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.ph_lev[i,j,k] = sg.hybrid_coord_calc_ph_lev(mesh.hyai[0,0,k], mesh.hybi[0,0,k], state.phs[i,j,0])
                #state.ph_exn_lev[i,j,k] = state.ph_lev[i,j,k]**cv.Rd_o_cpd
    
    for k in range(gm.half_kds,gm.half_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                #state.ph_lev[i,j,k] = sg.hybrid_coord_calc_ph_lev(mesh.hyai[0,0,k], mesh.hybi[0,0,k], state.phs[i,j,0])
                state.ph_exn_lev[i,j,k] = state.ph_lev[i,j,k]**cv.Rd_o_cpd
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.ph[i,j,k] = 0.5 * (state.ph_lev[i,j,k] + state.ph_lev[i,j,k+1])

@space_op
def calc_m(state:phv.HybridStateField,mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.m[i,j,k] = state.ph_lev[i,j,k+1] - state.ph_lev[i,j,k]
    
    for k in range(gm.half_kds+1,gm.half_kde-1+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.m_lev[i,j,k] = state.ph[i,j,k] - state.ph[i,j,k-1]

    #Top boundary
    for k in range(gm.half_kds,gm.half_kds+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.m_lev[i,j,k] = state.ph[i,j,k] - state.ph_lev[i,j,k]
    
    #Bottom boundary
    for k in range(gm.half_kde,gm.half_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.m_lev[i,j,k] = state.ph_lev[i,j,k] - state.ph[i,j,k-1]
 
    interp.averageCellToLonEdge(state.m,state.m_lon)
    interp.averageCellToLatEdge(state.m,state.m_lat)
    interp.interpCellToVtx(state.m,state.m_vtx,mesh)

@space_op
def calc_t(state:phv.HybridStateField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.t[i,j,k] = mf.temperature(state.pt[i,j,k],state.ph[i,j,k],state.qv[i,j,k])

@space_op
def calc_mf(state:phv.HybridStateField, advpt:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    
    accum_uv_cell(state,advpt,advPara,mesh,dt)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.mfx_lon[i,j,k] = state.m_lon[i,j,k] * state.u_lon[i,j,k]
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.mfy_lat[i,j,k] = state.m_lat[i,j,k] * state.v_lat[i,j,k]

    accum_mf_cell(state,advpt,advPara,mesh)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.mfx_lat[i,j,k] = mesh.half_tangent_wgt_0[0,j,0] * (state.mfx_lon[i-1,j  ,k] + state.mfx_lon[i,j  ,k]) + \
                                       mesh.half_tangent_wgt_1[0,j,0] * (state.mfx_lon[i-1,j+1,k] + state.mfx_lon[i,j+1,k])
                state.u_lat[i,j,k]   = state.mfx_lat[i,j,k] / state.m_lat[i,j,k]
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.mfy_lon[i,j,k] = mesh.full_tangent_wgt_0[0,j,0] * (state.mfy_lat[i,j-1,k] + state.mfy_lat[i+1,j-1,k]) + \
                                       mesh.full_tangent_wgt_1[0,j,0] * (state.mfy_lat[i,j  ,k] + state.mfy_lat[i+1,j  ,k])
                state.v_lon[i,j,k]   = state.mfy_lon[i,j,k] / state.m_lon[i,j,k]

@space_op
def calc_ke(state:phv.HybridStateField,mesh:gm.HybridMeshField):
    ke_vtx = Vector4(0.0, 0.0, 0.0, 0.0)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.ke[i,j,k] = (mesh.area_lon_west[0,j,0] * state.u_lon[i-1,j,k]**2 + 
                                   mesh.area_lon_east[0,j,0] * state.u_lon[i,j,k]**2 +
                                   mesh.area_lat_north[0,j-1,0]*state.v_lat[i,j-1,k]**2 +
                                   mesh.area_lat_south[0,j,0] * state.v_lat[i,j,k]**2) / mesh.area_cell[0,j,0]

    if gc.ke_scheme == 2:
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    # ke_vtx.e0 = (mesh.area_lat_east[0,j,0] * state.v_lat[i-1,j,k]**2 + 
                    #              mesh.area_lat_west[0,j,0] * state.v_lat[i,j,k]**2 +
                    #              mesh.area_lon_north[0,j,0]* state.u_lon[i-1,j,k]**2 +
                    #              mesh.area_lon_south[0,j+1,0]*state.u_lon[i-1,j+1,k]**2) / mesh.area_vtx[0,j,0]
                    # ke_vtx.e1 = (mesh.area_lat_east[0,j-1,0] * state.v_lat[i-1,j-1,k]**2 + 
                    #              mesh.area_lat_west[0,j-1,0] * state.v_lat[i,j-1,k]**2 +
                    #              mesh.area_lon_north[0,j-1,0]*state.u_lon[i-1,j-1,k]**2 +
                    #              mesh.area_lon_south[0,j,0]*state.u_lon[i-1,j,k]**2) / mesh.area_vtx[0,j-1,0]
                    # ke_vtx.e2 = (mesh.area_lat_east[0,j-1,0] * state.v_lat[i,j-1,k]**2 + 
                    #              mesh.area_lat_west[0,j-1,0] * state.v_lat[i+1,j-1,k]**2 +
                    #              mesh.area_lon_north[0,j-1,0]*state.u_lon[i,j-1,k]**2 +
                    #              mesh.area_lon_south[0,j,0]*state.u_lon[i,j,k]**2) / mesh.area_vtx[0,j-1,0]
                    # ke_vtx.e3 = (mesh.area_lat_east[0,j,0] * state.v_lat[i,j,k]**2 + 
                    #              mesh.area_lat_west[0,j,0] * state.v_lat[i+1,j,k]**2 +
                    #              mesh.area_lon_north[0,j,0]*state.u_lon[i,j,k]**2 +
                    #              mesh.area_lon_south[j+1]*state.u_lon[i,j+1,k]**2) / mesh.area_vtx[0,j,0]
                    # state.ke[i,j,k] = (1.0-gc.ke_cell_wgt)*(
                    #                     (ke_vtx.e0+ke_vtx.e3) * mesh.area_subcell_1[0,j,0] +
                    #                     (ke_vtx.e1+ke_vtx.e2) * mesh.area_subcell_0[0,j,0]
                    #                   ) / mesh.area_cell[0,j,0] + gc.ke_cell_wgt * state.ke[i,j,k]
                    state.ke[i,j,k] = (1.0-gc.ke_cell_wgt)*(
                                        (((mesh.area_lat_east[0,j,0] * state.v_lat[i-1,j,k]**2 + 
                                            mesh.area_lat_west[0,j,0] * state.v_lat[i,j,k]**2 +
                                            mesh.area_lon_north[0,j,0]* state.u_lon[i-1,j,k]**2 +
                                            mesh.area_lon_south[0,j+1,0]*state.u_lon[i-1,j+1,k]**2) / mesh.area_vtx[0,j,0]) + 
                                         ((mesh.area_lat_east[0,j,0] * state.v_lat[i,j,k]**2 + 
                                            mesh.area_lat_west[0,j,0] * state.v_lat[i+1,j,k]**2 +
                                            mesh.area_lon_north[0,j,0]*state.u_lon[i,j,k]**2 +
                                            mesh.area_lon_south[j+1]*state.u_lon[i,j+1,k]**2) / mesh.area_vtx[0,j,0])) * mesh.area_subcell_1[0,j,0] +
                                        (((mesh.area_lat_east[0,j-1,0] * state.v_lat[i-1,j-1,k]**2 + 
                                            mesh.area_lat_west[0,j-1,0] * state.v_lat[i,j-1,k]**2 +
                                            mesh.area_lon_north[0,j-1,0]*state.u_lon[i-1,j-1,k]**2 +
                                            mesh.area_lon_south[0,j,0]*state.u_lon[i-1,j,k]**2) / mesh.area_vtx[0,j-1,0]) + 
                                        ((mesh.area_lat_east[0,j-1,0] * state.v_lat[i,j-1,k]**2 + 
                                            mesh.area_lat_west[0,j-1,0] * state.v_lat[i+1,j-1,k]**2 +
                                            mesh.area_lon_north[0,j-1,0]*state.u_lon[i,j-1,k]**2 +
                                            mesh.area_lon_south[0,j,0]*state.u_lon[i,j,k]**2) / mesh.area_vtx[0,j-1,0])) * mesh.area_subcell_0[0,j,0]
                                      ) / mesh.area_cell[0,j,0] + gc.ke_cell_wgt * state.ke[i,j,k]

    #south reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = state.v_lat[i,j,k] ** 2

    state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.ke[i,j,k] = state.tmpsum[i,j,k] / gc.nlon

    #north reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = state.v_lat[i,j-1,k] ** 2
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.ke[i,j,k] = state.tmpsum[i,j,k] / gc.nlon
    
@space_op
def calc_vor(state:phv.HybridStateField, mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.vor[i,j,k] = (
                    (state.u_lon[i,j,k] * mesh.de_lon[0,j,0] - state.u_lon[i,j+1,k] * mesh.de_lon[0,j+1,0]) +
                    (state.v_lat[i+1,j,k]*mesh.de_lat[0,j,0] - state.v_lat[i,j,k]   * mesh.de_lat[0,j,0])
                    ) / mesh.area_vtx[0,j,0]
    
    if gc.pv_pole_stokes:
        #south reduce
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jds,gm.half_jds+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    state.tmpsum[i,j,k] = - state.u_lat[i,j,k] * mesh.le_lat[0,j,0]
    
        state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)
        
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jds,gm.half_jds+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    state.vor[i,j,k] = state.tmpsum[i,j,k] / gc.nlon / mesh.area_cell[0,j,0]

        #north reduce
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jde,gm.half_jde+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    state.tmpsum[i,j,k] = state.u_lat[i,j,k] * mesh.le_lat[0,j,0]

        state.tmpsum.sum(False,0,gm.full_nlev,False,gm.half_jde,gm.half_jde+1,True,0,gm.full_nlon)

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jde,gm.half_jde+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    state.vor[i,j,k] = state.tmpsum[i,j,k] / gc.nlon / mesh.area_cell[0,j+1,0]

@space_op
def calc_pv(state:phv.HybridStateField, mesh:gm.HybridMeshField):
    calc_vor(state,mesh)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.pv[i,j,k] = (state.vor[i,j,k] + mesh.half_f[0,j,0]) / state.m_vtx[i,j,k]

@space_op 
def interp_pv(state:phv.HybridStateField):
    if gc.pv_scheme == 2:
        interp_pv_upwind(state)

@space_op
def interp_pv_upwind(state:phv.HybridStateField):
    if gc.upwind_order_pv==3:
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    b = math.fabs(state.v_lon[i,j,k]) / (math.sqrt(state.u_lon[i,j,k]**2 + state.v_lon[i,j,k]**2) + cv.eps)
                    state.pv_lon[i,j,k] = b*mf.upwind3(mf.sign(1.0,state.v_lon[i,j,k]),gc.upwind_wgt_pv,state.pv[i,j-2,k],state.pv[i,j-1,k],state.pv[i,j,k],state.pv[i,j+1,k]) + \
                                          (1-b) * 0.5 * (state.pv[i,j-1,k] + state.pv[i,j,k])  

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jds,gm.half_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    b = math.fabs(state.u_lat[i,j,k]) / (math.sqrt(state.u_lat[i,j,k]**2 + state.v_lat[i,j,k]**2) + cv.eps)
                    state.pv_lat[i,j,k] = b*mf.upwind3(mf.sign(1.0,state.u_lat[i,j,k]),gc.upwind_wgt_pv,state.pv[i-2,j,k],state.pv[i-1,j,k],state.pv[i,j,k],state.pv[i+1,j,k]) + \
                                          (1-b) * 0.5 * (state.pv[i-1,j,k] + state.pv[i,j,k]) 

@space_op
def calc_div(state:phv.HybridStateField, mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.div[i,j,k] = (
                    (state.u_lon[i,j,k] * mesh.le_lon[0,j,0] - state.u_lon[i-1,j,k] * mesh.le_lon[0,j,0]) +
                    (state.v_lat[i,j,k] * mesh.le_lat[0,j,0] - state.v_lat[i,j-1,k] * mesh.le_lat[0,j-1,0])
                ) / mesh.area_cell[0,j,0]
    
    #south reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = state.v_lat[i,j,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.div[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0]
    
    #north reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = - state.v_lat[i,j-1,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.div[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0]

@space_op
def calc_gz_lev(state:phv.HybridStateField,staticv:phv.HybridStaticField):
    for k in range(gm.half_kds,gm.half_kde-1+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                dgz= Float64(0.0)
                #调整为左闭右开
                for l in range(k,gm.full_nlev-1+1):
                    dgz += cv.Rd * state.t[i,j,l] * math.log(state.ph_lev[i,j,l+1] / state.ph_lev[i,j,l])
                state.gz_lev[i,j,k] = staticv.gzs[i,j,0] + dgz       

@space_op
def calc_grad_mf(state:phv.HybridStateField,tend:phv.HybridTendField,mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dmfdlon[i,j,k] = (state.mfx_lon[i,j,k] - state.mfx_lon[i-1,j,k]) * mesh.le_lon[0,j,0] / mesh.area_cell[0,j,0]
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dmfdlat[i,j,k] = (state.mfy_lat[i,j,k] * mesh.le_lat[0,j,0] - state.mfy_lat[i,j-1,k] * mesh.le_lat[0,j-1,0]) / mesh.area_cell[0,j,0]

    #south reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = state.mfy_lat[i,j,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dmfdlat[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0]

    #north reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = - state.mfy_lat[i,j-1,k]

    state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dmfdlat[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0]

@space_op
def calc_dphsdt(tend:phv.HybridTendField):
    for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dphs[i,j,0] = 0.0
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dphs[i,j,0] = tend.dphs[i,j,0] - tend.dmfdlon[i,j,k] - tend.dmfdlat[i,j,k]

@space_op
def calc_we_lev(state:phv.HybridStateField, tend:phv.HybridTendField, advpt:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    
    #ToRemember 目前无kernel概念，同一变量的连续多个for循环只有最后一个会保留update,要把三维留在最后
    #Set vertical boundary conditions.
    for k in range(gm.half_kds,gm.half_kds+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.we_lev[i,j,k] = 0.0
    for k in range(gm.half_kde,gm.half_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.we_lev[i,j,k] = 0.0
    
    #ToRemember 目前循环内的坐标运算无法自动添加halo，只会自动添加在循环变量的范围上
    for k in range(gm.half_kds+1,gm.half_kde-1+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                mfs = Float64(0.0)
                for l in range(0, k - 1 + 1):
                    mfs = mfs + tend.dmfdlon[i,j,l] + tend.dmfdlat[i,j,l]
                state.we_lev[i,j,k] = - mesh.hybi[0,0,k] * tend.dphs[i,j,0] - mfs
                # state.we_lev[i,j,k] = - sg.hybrid_coord_calc_dphdt_lev(mesh.hybi[0,0,k],tend.dphs[i,j,0]) - mfs


    # for j in range(gm.full_jds,gm.full_jde+1):
    #     for i in range(gm.full_ids,gm.full_ide+1):
    #         mfs = Float64(0.0)
    #         for k in range(gm.half_kds+1,gm.half_kde-1+1):                  
    #             for l in range(k - 1, k - 1 + 1):
    #                 mfs = mfs + tend.dmfdlon[i,j,l] + tend.dmfdlat[i,j,l]
    #             state.we_lev[i,j,k] = - mesh.hybi[0,0,k] * tend.dphs[i,j,0] - mfs
                # state.we_lev[i,j,k] = - sg.hybrid_coord_calc_dphdt_lev(mesh.hybi[0,0,k],tend.dphs[i,j,0]) - mfs



    accum_we_lev(state,advpt,advPara,dt)
    
    interp.interp_lev_edge_to_lev_lon_edge(mesh,state.we_lev,state.we_lev_lon)
    interp.interp_lev_edge_to_lev_lat_edge(mesh,state.we_lev,state.we_lev_lat)

@space_op
def calc_wedudlev_wedvdlev(state:phv.HybridStateField,tend:phv.HybridTendField):
    for k in range(gm.full_kds+1,gm.full_kde-1 +1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.wedudlev[i,j,k] = (
                    state.we_lev_lon[i,j,k+1] * (state.u_lon[i,j,k+1] - state.u_lon[i,j,k  ]) +
                    state.we_lev_lon[i,j,k  ] * (state.u_lon[i,j,k  ] - state.u_lon[i,j,k-1])
                ) / state.m_lon[i,j,k] / 2.0
    
    for k in range(gm.full_kds,gm.full_kds+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.wedudlev[i,j,k] = (state.we_lev_lon[i,j,k+1] *
                        (state.u_lon[i,j,k+1] - state.u_lon[i,j,k])) / state.m_lon[i,j,k] / 2.0

    for k in range(gm.full_kde,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.wedudlev[i,j,k] = (state.we_lev_lon[i,j,k] *
                        (state.u_lon[i,j,k] - state.u_lon[i,j,k-1])) / state.m_lon[i,j,k] / 2.0

    for k in range(gm.full_kds+1,gm.full_kde-1 +1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.wedvdlev[i,j,k] = (
                    state.we_lev_lat[i,j,k+1] * (state.v_lat[i,j,k+1] - state.v_lat[i,j,k  ]) +
                    state.we_lev_lat[i,j,k  ] * (state.v_lat[i,j,k  ] - state.v_lat[i,j,k-1])
                ) / state.m_lat[i,j,k] / 2.0

    for k in range(gm.full_kds,gm.full_kds +1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.wedvdlev[i,j,k] = (state.we_lev_lat[i,j,k+1] * 
                    (state.v_lat[i,j,k+1] - state.v_lat[i,j,k]) / state.m_lat[i,j,k] / 2.0)
    
    for k in range(gm.full_kde,gm.full_kde +1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.wedvdlev[i,j,k] = (state.we_lev_lat[i,j,k] * 
                    (state.v_lat[i,j,k] - state.v_lat[i,j,k-1]) / state.m_lat[i,j,k] / 2.0)

@space_op
def calc_grad_ptf(state:phv.HybridStateField,tend:phv.HybridTendField,adv:phv.HybridAdvField, mesh:gm.HybridMeshField, dt:Float64):
    
    fs.ffsl_calc_tracer_hflx(state,adv,state.pt,state.ptf_lon,state.ptf_lat,mesh,dt)
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dptfdlon[i,j,k] = (
                    state.ptf_lon[i,j,k] - state.ptf_lon[i-1,j,k]    
                ) * mesh.le_lon[0,j,0] / mesh.area_cell[0,j,0]
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dptfdlat[i,j,k] = (
                    state.ptf_lat[i,j,k]   * mesh.le_lat[0,j,0] - 
                    state.ptf_lat[i,j-1,k] * mesh.le_lat[0,j-1,0]
                ) / mesh.area_cell[0,j,0]

    #south reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = state.ptf_lat[i,j,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dptfdlat[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0]
    
    #north reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = - state.ptf_lat[i,j-1,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dptfdlat[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0]

    #---------------FFSL---------------
    #Set upper and lower boundary conditions.
    for k in range(gm.full_kds-1,gm.full_kds):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.pt[i,j,k] = 2 * state.pt[i,j,k+1] - state.pt[i,j,k+2]

    for k in range(gm.full_kds-2,gm.full_kds-1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.pt[i,j,k] = 2 * state.pt[i,j,k+1] - state.pt[i,j,k+2]

    for k in range(gm.full_kde+1,gm.full_kde+2+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.pt[i,j,k] = 2 * state.pt[i,j,k-1] - state.pt[i,j,k-2]

    fs.ffsl_calc_tracer_vflx(adv,state.pt,state.ptf_lev)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dptfdlev[i,j,k] = state.ptf_lev[i,j,k+1] - state.ptf_lev[i,j,k]

@space_op
def calc_coriolis(state:phv.HybridStateField,tend:phv.HybridTendField,mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                if gc.coriolis_scheme == 1:
                    tend.qhu[i,j,k] = (
                        mesh.half_tangent_wgt_0[0,j,0] * (
                            state.mfx_lon[i-1,j,k] * (state.pv_lat[i,j,k] + state.pv_lon[i-1,j,k]) +
                            state.mfx_lon[i  ,j,k] * (state.pv_lat[i,j,k] + state.pv_lon[i  ,j,k])
                        )  +
                        mesh.half_tangent_wgt_1[0,j,0] * (
                            state.mfx_lon[i-1,j+1,k] * (state.pv_lat[i,j,k] + state.pv_lon[i-1,j+1,k]) +
                            state.mfx_lon[i  ,j+1,k] * (state.pv_lat[i,j,k] + state.pv_lon[i  ,j+1,k])
                        )
                    ) * 0.5
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                if gc.coriolis_scheme == 1:
                    tend.qhv[i,j,k] = (
                        mesh.full_tangent_wgt_0[0,j,0] * (
                            state.mfy_lat[i  ,j-1,k] * (state.pv_lon[i,j,k] + state.pv_lat[i  ,j-1,k]) +
                            state.mfy_lat[i+1,j-1,k] * (state.pv_lon[i,j,k] + state.pv_lat[i+1,j-1,k])
                        )  +
                        mesh.full_tangent_wgt_1[0,j,0] * (
                            state.mfy_lat[i  ,j  ,k] * (state.pv_lon[i,j,k] + state.pv_lat[i  ,j  ,k]) +
                            state.mfy_lat[i+1,j  ,k] * (state.pv_lon[i,j,k] + state.pv_lat[i+1,j  ,k])
                        )
                    ) * 0.5

@space_op
def calc_grad_ke(state:phv.HybridStateField,tend:phv.HybridTendField,mesh:gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.dkedlon[i,j,k] = (state.ke[i+1,j,k] - state.ke[i,j,k]) / mesh.de_lon[0,j,0]
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dkedlat[i,j,k] = (state.ke[i,j+1,k] - state.ke[i,j,k]) / mesh.de_lat[0,j,0]

@space_op
def pgf_lin97(state: phv.HybridStateField, tend: phv.HybridTendField, mesh: gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tl = 1 + 0.5 * (state.qm[i,j,k] + state.qm[i+1,j,k])
                dph1 = state.ph_exn_lev[i+1,j,k+1] - state.ph_exn_lev[i,j,k]
                dph2 = state.ph_exn_lev[i,j,k+1] - state.ph_exn_lev[i+1,j,k]
                dgz1 = state.gz_lev[i,j,k+1] - state.gz_lev[i+1,j,k]
                dgz2 = state.gz_lev[i,j,k] - state.gz_lev[i+1,j,k+1]
                tend.pgf_lon[i,j,k] = -(dph1 * dgz1 + dph2 * dgz2) / mesh.de_lon[0,j,0] / (dph1 + dph2) / tl

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tl = 1 + 0.5 * (state.qm[i,j,k] + state.qm[i,j+1,k])
                dph1 = state.ph_exn_lev[i,j+1,k+1] - state.ph_exn_lev[i,j,k]
                dph2 = state.ph_exn_lev[i,j,k+1] - state.ph_exn_lev[i,j+1,k]
                dgz1 = state.gz_lev[i,j,k+1] - state.gz_lev[i,j+1,k]
                dgz2 = state.gz_lev[i,j,k] - state.gz_lev[i,j+1,k+1]
                tend.pgf_lat[i,j,k] = -(dph1 * dgz1 + dph2 * dgz2) / mesh.de_lat[0,j,0] / (dph1 + dph2) / tl

@space_op
def calc_tend_forward(tend: phv.HybridTendField, tendPara: phv.HybridTendPara):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.du[i,j,k] = tend.qhv[i,j,k] - tend.dkedlon[i,j,k] - tend.wedudlev[i,j,k]
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dv[i,j,k] = - tend.qhu[i,j,k] - tend.dkedlat[i,j,k] - tend.wedvdlev[i,j,k]

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.dpt[i,j,k] = -tend.dptfdlon[i,j,k] - tend.dptfdlat[i,j,k] - tend.dptfdlev[i,j,k]
    
    tendPara.phs = True
    tendPara.pt = True

@space_op
def calc_tend_backward(tend1: phv.HybridTendField, tend2: phv.HybridTendField, tendPara: phv.HybridTendPara):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend1.du[i,j,k] = tend2.du[i,j,k] - tend1.pgf_lon[i,j,k]
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend1.dv[i,j,k] = tend2.dv[i,j,k] - tend1.pgf_lat[i,j,k]

    tendPara.u = True
    tendPara.v = True

@space_op
def update_state(old_state: phv.HybridStateField, new_state: phv.HybridStateField, tend: phv.HybridTendField, tendPara: phv.HybridTendPara, mesh: gm.HybridMeshField, dt: Float64):
    if tendPara.phs:
        filter.Filter_On_Cell(0,0,tend.dphs)
        #filter_on_cell(block%big_filter, dtend%dphs)
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                new_state.phs[i,j,0] = old_state.phs[i,j,0] + dt * tend.dphs[i,j,0]
        
        calc_ph(new_state, mesh)
        calc_m(new_state, mesh)
    
    if tendPara.pt:
        filter.Filter_On_Cell(0,1,tend.dpt)
        #filter_on_cell(block%big_filter, dtend%dpt)
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    new_state.pt[i,j,k] = (old_state.pt[i,j,k] * old_state.m[i,j,k] + dt * tend.dpt[i,j,k]) / new_state.m[i,j,k]
    
    if tendPara.gz:
        #filter_on_cell(block%big_filter, dtend%dgz)
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    new_state.gz[i,j,k] = old_state.gz[i,j,k] + dt * tend.dgz[i,j,k]
    
        calc_m(new_state, mesh)

    if tendPara.u and tendPara.v:
        filter.Filter_on_lon_edge(0,tend.du)
        filter.Filter_on_lat_edge(0,tend.dv)
        #filter_on_lon_edge(block%big_filter, dtend%du)
        #filter_on_lat_edge(block%big_filter, dtend%dv)
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    new_state.u_lon[i,j,k] = old_state.u_lon[i,j,k] + dt * tend.du[i,j,k]

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jds,gm.half_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    new_state.v_lat[i,j,k] = old_state.v_lat[i,j,k] + dt * tend.dv[i,j,k]

@space_op
def div_damp_run(state: phv.HybridStateField, mesh: gm.HybridMeshField):
    for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
                for i in range(gm.half_ids,gm.half_ide+1):
                    state.u_lon[i,j,k] = state.u_lon[i,j,k] + mesh.c_lon[0,j,k] * (state.div[i+1,j,k] - state.div[i,j,k]) / mesh.de_lon[0,j,0]
    
    for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.half_jds,gm.half_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state.v_lat[i,j,k] = state.v_lat[i,j,k] + mesh.c_lat[0,j,k] * (state.div[i,j+1,k] - state.div[i,j,k]) / mesh.de_lat[0,j,0]

@space_op
def smag_damp_run(state: phv.HybridStateField, tend: phv.HybridTendField, mesh: gm.HybridMeshField, dt:Float64):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.smag_t[i,j,k] = (
                    state.u_lon[i,j,k] - state.u_lon[i-1,j,k]
                ) / mesh.de_lon[0,j,0] - (
                    state.v_lat[i,j,k] * mesh.half_cos_lat[0,j,0] - 
                    state.v_lat[i,j-1,k] * mesh.half_cos_lat[0,j-1,0]
                ) / mesh.le_lon[0,j,0] / mesh.full_cos_lat[0,j,0]
    
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.smag_s[i,j,k] = (
                    state.v_lat[i+1,j,k] - state.v_lat[i,j,k] 
                ) / mesh.le_lat[0,j,0] + (
                    state.u_lon[i,j+1,k] * mesh.full_cos_lat[0,j+1,0] - 
                    state.u_lon[i,j,k] * mesh.full_cos_lat[0,j,0]
                ) / mesh.de_lat[0,j,0] / mesh.half_cos_lat[0,j,0]

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.kmh_lon[i,j,k] = (
                    gc.smag_damp_coef / (1.0 / mesh.de_lon[0,j,0]**2 + 1 / mesh.le_lon[0,j,0]**2)
                ) * math.sqrt(
                    0.5 * (state.smag_t[i,j,k]**2 + state.smag_t[i+1,j,k]**2) + 
                    0.5 * (state.smag_s[i,j,k]**2 + state.smag_s[i,j-1,k]**2)
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.kmh_lat[i,j,k] = (
                    gc.smag_damp_coef / (1.0 / mesh.le_lat[0,j,0]**2 + 1 / mesh.de_lat[0,j,0]**2)
                ) * math.sqrt(
                    0.5 * (state.smag_t[i,j,k]**2 + state.smag_t[i,j+1,k]**2) + 
                    0.5 * (state.smag_s[i,j,k]**2 + state.smag_s[i-1,j,k]**2)
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                tend.smag_dudt[i,j,k] = state.kmh_lon[i,j,k] * (
                    (state.u_lon[i-1,j,k] - 2*state.u_lon[i,j,k] + state.u_lon[i+1,j,k]) / mesh.de_lon[0,j,0]**2 +
                    ((state.u_lon[i,j+1,k] - state.u_lon[i,j,k]) / mesh.de_lat[0,j,0] * mesh.half_cos_lat[0,j,0] - 
                     (state.u_lon[i,j,k] - state.u_lon[i,j-1,k]) / mesh.de_lat[0,j-1,0] * mesh.half_cos_lat[0,j-1,0]
                    ) / mesh.le_lon[0,j,0] / mesh.full_cos_lat[0,j,0]
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                state.u_lon[i,j,k] = state.u_lon[i,j,k] + dt * tend.smag_dudt[i,j,k]

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.smag_dvdt[i,j,k] = state.kmh_lat[i,j,k] * (
                    (state.v_lat[i-1,j,k] - 2*state.v_lat[i,j,k] + state.v_lat[i+1,j,k]) / mesh.le_lat[0,j,0]**2
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds+1,gm.half_jde-1+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.smag_dvdt[i,j,k] = state.kmh_lat[i,j,k] * (
                    (state.v_lat[i-1,j,k] - 2*state.v_lat[i,j,k] + state.v_lat[i+1,j,k]) / mesh.le_lat[0,j,0]**2 +
                    ((state.v_lat[i,j+1,k] - state.v_lat[i,j,k]) / mesh.le_lon[0,j+1,0] * mesh.full_cos_lat[0,j+1,0] -
                     (state.v_lat[i,j,k] - state.v_lat[i,j-1,k]) / mesh.le_lon[0,j,0] * mesh.full_cos_lat[0,j,0]
                     ) / mesh.de_lat[0,j,0] / mesh.half_cos_lat[0,j,0]
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jde,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                tend.smag_dvdt[i,j,k] = state.kmh_lat[i,j,k] * (
                    (state.v_lat[i-1,j,k] - 2*state.v_lat[i,j,k] + state.v_lat[i+1,j,k]) / mesh.le_lat[0,j,0]**2
                )

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.v_lat[i,j,k] = state.v_lat[i,j,k] + dt * tend.smag_dvdt[i,j,k]

@space_op
def trickypt(state: phv.HybridStateField):
    for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state.pt[i,j,k] = state.pt[i,j,k] * state.m[i,j,k]

@space_op
def pole_damp_run(state: phv.HybridStateField, mesh: gm.HybridMeshField):

    #ToDo 去除reduce sum造成的多次update,此处先简单拆开
    trickypt(state)
    # for k in range(gm.full_kds,gm.full_kde+1):
    #         for j in range(gm.full_jds,gm.full_jde+1):
    #             for i in range(gm.full_ids,gm.full_ide+1):
    #                 state.pt[i,j,k] = state.pt[i,j,k] * state.m[i,j,k]
    
    filter.Filter_On_Cell(1,0,state.phs)
    #call filter_on_cell(block%small_filter_phs, dstate%phs)

    calc_ph(state, mesh)
    calc_m(state,mesh)

    filter.Filter_On_Cell(2,1,state.pt)
    #call filter_on_cell(block%small_filter_pt, dstate%pt)

    for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state.pt[i,j,k] = state.pt[i,j,k] / state.m[i,j,k]

    filter.Filter_on_lon_edge(3,state.u_lon)
    filter.Filter_on_lat_edge(3,state.v_lat)
    #call filter_on_lon_edge(block%small_filter_uv, dstate%u_lon)
    #call filter_on_lat_edge(block%small_filter_uv, dstate%v_lat)

@space_op
def accum_uv_cell(state:phv.HybridStateField, adv:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField, dt: Float64):
    if advPara.uv_step == -1:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.uu[i,j,k] = adv.u0[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.vv[i,j,k] = adv.v0[i,j,k]
        advPara.uv_step = 1
    
    if advPara.uv_step == 0:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.uu[i,j,k] = state.u_lon[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.vv[i,j,k] = state.v_lat[i,j,k]
    elif advPara.uv_step == advPara.nstep:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.uu[i,j,k] = (adv.uu[i,j,k] + state.u_lon[i,j,k]) / (advPara.nstep + 1)
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.vv[i,j,k] = (adv.vv[i,j,k] + state.v_lat[i,j,k]) / (advPara.nstep + 1)
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.u0[i,j,k] = adv.uu[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.v0[i,j,k] = adv.vv[i,j,k]
    else:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.uu[i,j,k] = adv.uu[i,j,k] + state.u_lon[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.vv[i,j,k] = adv.vv[i,j,k] + state.v_lat[i,j,k]

    if advPara.dynamic:
        advPara.uv_step = 0
    else:
        advPara.uv_step += 1

    if advPara.dynamic or advPara.uv_step > advPara.nstep:
        if not advPara.dynamic:
            advPara.uv_step = -1
        
        #eul
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds_no_pole, gm.full_jde_no_pole + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.cflx[i,j,k] = adv.uu[i,j,k] * dt / mesh.de_lon[0,j,0]
        
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids , gm.full_ide + 1):
                    adv.cfly[i,j,k] = adv.vv[i,j,k] * dt / mesh.de_lat[0,j,0]
    
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds_no_pole, gm.full_jde_no_pole + 1):
                for i in range(gm.full_ids , gm.full_ide + 1):
                    adv.divx[i,j,k] = (adv.uu[i,j,k] - adv.uu[i-1,j,k]) * mesh.le_lon[0,j,0] / mesh.area_cell[0,j,0]
                    adv.divy[i,j,k] = (adv.vv[i,j,k] * mesh.le_lat[0,j,0] - adv.vv[i,j-1,k] * mesh.le_lat[0,j-1,0]) / mesh.area_cell[0,j,0]

        #south reduce
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jds + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state.tmpsum[i, j, k] = adv.vv[i, j, k]
        
        state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jds + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.divy[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0]

        #north reduce
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jde, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state.tmpsum[i,j,k] = - adv.vv[i,j-1,k]
        
        state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jde, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.divy[i,j,k] = state.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0]

@space_op
def accum_mf_cell(state:phv.HybridStateField, adv:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField):
    if advPara.mf_step == -1:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.mfx[i,j,k] = 0.0

        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.mfy[i,j,k] = 0.0

        advPara.mf_step = 1
    
    if advPara.mf_step == 0:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.mfx[i,j,k] = state.mfx_lon[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.mfy[i,j,k] = state.mfy_lat[i,j,k]
    elif advPara.mf_step == advPara.nstep:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.mfx[i,j,k] = (adv.mfx[i,j,k] + state.mfx_lon[i,j,k]) / advPara.nstep
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.mfy[i,j,k] = (adv.mfy[i,j,k] + state.mfy_lat[i,j,k]) / advPara.nstep
    else:
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.half_ids, gm.half_ide + 1):
                    adv.mfx[i,j,k] = adv.mfx[i,j,k] + state.mfx_lon[i,j,k]
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.half_jds, gm.half_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.mfy[i,j,k] = adv.mfy[i,j,k] + state.mfy_lat[i,j,k]

    if advPara.dynamic:
        advPara.mf_step = 0
    else:
        advPara.mf_step += 1

    if not advPara.dynamic and advPara.mf_step > advPara.nstep:
        advPara.mf_step = -1

@space_op
def accum_we_lev(state:phv.HybridStateField, adv:phv.HybridAdvField, advPara:phv.HybridAdvPara, dt: Float64):
    if advPara.we_step == -1:
        for k in range(gm.half_kds, gm.half_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.we[i,j,k] = adv.we0[i,j,k]
                    adv.mm[i,j,k] = adv.m0[i,j,k]
        
        advPara.we_step = 1
    
    if advPara.we_step == 0:
        for k in range(gm.half_kds, gm.half_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.we[i,j,k] = state.we_lev[i,j,k]
                    adv.mm[i,j,k] = state.m_lev[i,j,k]
    elif advPara.we_step == advPara.nstep:
        for k in range(gm.half_kds, gm.half_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.we[i,j,k] = (adv.we[i,j,k] + state.we_lev[i,j,k]) / (advPara.nstep + 1)
                    adv.mm[i,j,k] = (adv.mm[i,j,k] + state.m_lev[i,j,k]) / (advPara.nstep + 1)
                    adv.we0[i,j,k] = state.we_lev[i,j,k]
                    adv.m0[i,j,k] = state.m_lev[i,j,k]
    else:
        for k in range(gm.half_kds, gm.half_kde + 1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.we[i,j,k] = adv.we[i,j,k] + state.we_lev[i,j,k]
                    adv.mm[i,j,k] = adv.mm[i,j,k] + state.m_lev[i,j,k]
    
    if advPara.dynamic:
        advPara.we_step = 0
    else:
        advPara.we_step = advPara.we_step + 1
    
    if (advPara.dynamic or advPara.we_step > advPara.nstep):
        if advPara.dynamic == False:
            advPara.we_step = -1
        for k in range(gm.half_kds+1, gm.half_kde -1+1):
            for j in range(gm.full_jds, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    adv.cflz[i,j,k] = adv.we[i,j,k] / adv.mm[i,j,k] * dt

@space_op
def copy_old_m(state:phv.HybridStateField, adv:phv.HybridAdvField):
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                adv.old_m[i,j,k] = state.m[i,j,k]

#ToCheck accum not in field,but in hybridfield
@space_op
def adv_accum_wind(state:phv.HybridStateField, adv:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField, dt: Float64):
    accum_uv_cell(state, adv, advPara, mesh, dt)
    accum_mf_cell(state, adv, advPara, mesh)
    accum_we_lev(state, adv, advPara, dt)

@space_op
def adv_run(state_old:phv.HybridStateField, state_new:phv.HybridStateField, adv:phv.HybridAdvField, advPara:phv.HybridAdvPara, mesh:gm.HybridMeshField, dt: Float64):
    adv_accum_wind(state_old, adv, advPara, mesh, dt)
    fs.ffsl_calc_tracer_hflx(state_old, adv, state_old.qv, adv.qmf_lon, adv.qmf_lat, mesh, dt)

    if gc.dt_adv == 0:
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = adv.old_m[i,j,k] * state_old.qv[i,j,k] / state_old.m[i,j,k]
        
        fs.ffsl_calc_tracer_vflx(adv, state_new.qv, adv.qmf_lev)

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = state_new.qv[i,j,k] * state_old.m[i,j,k]
        
        filter.Filter_On_Cell(1,1,state_new.qv)
        #call filter_on_cell(block%small_filter_phs, q(:,:,:,:,l,new))

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = state_new.qv[i,j,k] / state_old.m[i,j,k]
    else:
        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds+1,gm.full_jde):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = adv.old_m[i,j,k] * state_old.qv[i,j,k] - (
                        adv.qmf_lon[i,j,k] - adv.qmf_lon[i-1,j,k]
                    ) * mesh.le_lon[0,j,0] + (
                        adv.qmf_lat[i,j,k] * mesh.le_lat[0,j,0] -
                        adv.qmf_lat[i,j-1,k] * mesh.le_lat[0,j-1,0]
                    ) / mesh.area_cell[0,j,0] * gc.dt_adv

        #south reduce
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jds + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state_new.tmpsum[i, j, k] = adv.qmf_lat[i,j,k]

        state_new.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jds, gm.full_jds + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state_new.qv[i,j,k] = adv.old_m[i,j,k] * state_old.qv[i,j,k] - (
                        state_new.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0] * gc.dt_adv
                    )

        #north reduce
        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jde, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state_new.tmpsum[i, j, k] = adv.qmf_lat[i,j-1,k]
        
        state_new.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)

        for k in range(gm.full_kds, gm.full_kde + 1):
            for j in range(gm.full_jde, gm.full_jde + 1):
                for i in range(gm.full_ids, gm.full_ide + 1):
                    state_new.qv[i,j,k] = adv.old_m[i,j,k] * state_old.qv[i,j,k] + (
                        state_new.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0] * gc.dt_adv
                    )

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = state_new.qv[i,j,k] / state_old.m[i,j,k]

        fs.ffsl_calc_tracer_vflx(adv, state_new.qv, adv.qmf_lev)

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = state_new.qv[i,j,k] * state_old.m[i,j,k] - (adv.qmf_lev[i,j,k+1] - adv.qmf_lev[i,j,k]) * gc.dt_adv

        filter.Filter_On_Cell(1,1,state_new.qv)

        for k in range(gm.full_kds,gm.full_kde+1):
            for j in range(gm.full_jds,gm.full_jde+1):
                for i in range(gm.full_ids,gm.full_ide+1):
                    state_new.qv[i,j,k] = state_new.qv[i,j,k] / state_old.m[i,j,k]

    copy_old_m(state_old, adv)

@space_op
def reset_flag(tendPara: phv.HybridTendPara):
    tendPara.u = False
    tendPara.v = False
    tendPara.pt = False
    tendPara.gz = False
    tendPara.phs = False

@space_op
def timeadvance():
    Time_Advance()

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 0, parallel = True)
def Time_Advance():
    pass

@space_op
def Diagnose_State(state:phv.HybridStateField):
    Diagnose(state,15)

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 0, parallel = True)
def Diagnose(state:phv.HybridStateField, xstep:Int32):
    pass
