import math

from frontend.pytoy.lang import space_op, Float64, Int32
from frontend.pytoy.lang.dtype_ext import LonLatField
import frontend.pytoy.lang.builtin_func

import gmcore_smallball.global_config as gc
import gmcore_smallball.physical_variable as phv
import gmcore_smallball.mesh as gm
import gmcore_smallball.const_value as cv
import gmcore_smallball.sphere_geometry as sg
import gmcore_smallball.interp as interp
import gmcore_smallball.math_formual as mf
from gmcore_smallball.utils import Vector4, Vector3

@space_op
def ffsl_calc_tracer_hflx(state: phv.HybridStateField, adv: phv.HybridAdvField, q: LonLatField[Float64], qmfx: LonLatField[Float64], qmfy: LonLatField[Float64],mesh: gm.HybridMeshField, dt: Float64):
    
    hflx_ppm_inner(adv, adv.uu, adv.vv, q, q, qmfx, qmfy)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                adv.qx[i,j,k] = q[i,j,k] - 0.5 * (
                    (qmfx[i,j,k] - qmfx[i-1,j,k]) * 
                    mesh.le_lon[0,j,0] / mesh.area_cell[0,j,0] - 
                    adv.divx[i,j,k] * q[i,j,k]
                ) * dt
                adv.qy[i,j,k] = q[i,j,k] - 0.5 * (
                    (qmfy[i,j,k]   * mesh.le_lat[0,j,0] - 
                     qmfy[i,j-1,k] * mesh.le_lat[0,j-1,0]) / mesh.area_cell[0,j,0] - 
                     adv.divy[i,j,k] * q[i,j,k]
                ) * dt
    
    #south reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = qmfy[i,j,k]
    
    state.tmpsum.sum(False,0,gm.full_nlev,False,0,1,True,0,gm.full_nlon)

    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jds+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                adv.qx[i,j,k] = q[i,j,k]
                adv.qy[i,j,k] = q[i,j,k] + 0.5 * (
                    state.tmpsum[i,j,k] * mesh.le_lat[0,j,0] / gc.nlon / mesh.area_cell[0,j,0] -
                    adv.divy[i,j,k] * q[i,j,k]) * dt

    #north reduce
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                state.tmpsum[i,j,k] = qmfy[i,j-1,k]

    state.tmpsum.sum(False,0,gm.full_nlev,False,gm.full_jde,gm.full_jde+1,True,0,gm.full_nlon)


    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jde,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                adv.qx[i,j,k] = q[i,j,k]
                adv.qy[i,j,k] = q[i,j,k] + 0.5 * (
                                state.tmpsum[i,j,k] * mesh.le_lat[0,j-1,0] / gc.nlon / mesh.area_cell[0,j,0] -
                                adv.divy[i,j,k] * q[i,j,k]) * dt

    hflx_ppm_outer(adv, adv.mfx, adv.mfy, adv.qy, adv.qx, qmfx, qmfy)

@space_op
def ffsl_calc_tracer_vflx(adv: phv.HybridAdvField, q: LonLatField[Float64], qmfz: LonLatField[Float64]):
    vflx_ppm(adv,adv.we,q,qmfz)

#Todo 更改减少updatehalo的逻辑,目前逻辑为一个kernel内对同一变量只update一次,但一kernel在一段流程前后分别调用时，前一次的updata会被吃掉
@space_op
def hflx_ppm_inner(adv: phv.HybridAdvField, u: LonLatField[Float64], v: LonLatField[Float64], mx: LonLatField[Float64], my: LonLatField[Float64], mfx: LonLatField[Float64], mfy: LonLatField[Float64]):
    buf = Vector3(0.0, 0.0, 0.0)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                buf = mf.ppm(mx[i-2,j,k], mx[i-1,j,k], mx[i,j,k], mx[i+1,j,k], mx[i+2,j,k])
                adv.qlx[i,j,k] = buf.x
                adv.dqx[i,j,k] = buf.y
                adv.q6x[i,j,k] = buf.z
                buf = mf.ppm(my[i,j-2,k], my[i,j-1,k], my[i,j,k], my[i,j+1,k], my[i,j+2,k])
                adv.qly[i,j,k] = buf.x
                adv.dqy[i,j,k] = buf.y
                adv.q6y[i,j,k] = buf.z
    
    cursum = Float64(0.0)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                ci = math.trunc(adv.cflx[i,j,k])
                cf = adv.cflx[i,j,k] - ci
                if (adv.cflx[i,j,k] > 0.0):
                    s1 = 1.0 - cf
                    s2 = 1.0
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    cursum = 0.0
                    for l in range(i+1-ci,i+1):
                        cursum = cursum + mx[l,j,k]
                    mfx[i,j,k] = u[i,j,k] * (cursum + adv.qlx[i - ci,j,k] * ds1 + 0.5 * adv.dqx[i - ci,j,k] * ds2 + adv.q6x[i - ci,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflx[i,j,k]
                elif (adv.cflx[i,j,k] < 0.0):
                    s1 = 0.0
                    s2 = -cf
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    cursum = 0.0
                    for l in range(i+1,i-ci+1):
                        cursum = cursum + mx[l,j,k]
                    mfx[i,j,k] = -u[i,j,k] * (cursum + adv.qlx[i - ci + 1,j,k] * ds1 + 0.5 * adv.dqx[i - ci + 1,j,k] * ds2 + adv.q6x[i - ci + 1,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflx[i,j,k]
                else:
                    mfx[i,j,k] = 0.0
                
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                if (adv.cfly[i,j,k] > 0):
                    s1 = 1 - adv.cfly[i,j,k]
                    s2 = 1
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    mfy[i,j,k] =  v[i,j,k] * (adv.qly[i,j,k] * ds1 + 0.5 * adv.dqy[i,j,k] * ds2 + adv.q6y[i,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cfly[i,j,k]
                elif (adv.cfly[i,j,k] < 0):
                    s1 = 0
                    s2 = -adv.cfly[i,j,k]
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    mfy[i,j,k] = -v[i,j,k] * (adv.qly[i,j + 1,k] * ds1 + 0.5 * adv.dqy[i,j + 1,k] * ds2 + adv.q6y[i,j + 1,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cfly[i,j,k]
                else:
                    mfy[i,j,k] = 0

@space_op
def hflx_ppm_outer(adv: phv.HybridAdvField, u: LonLatField[Float64], v: LonLatField[Float64], mx: LonLatField[Float64], my: LonLatField[Float64], mfx: LonLatField[Float64], mfy: LonLatField[Float64]):
    buf = Vector3(0.0, 0.0, 0.0)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                buf = mf.ppm(mx[i-2,j,k], mx[i-1,j,k], mx[i,j,k], mx[i+1,j,k], mx[i+2,j,k])
                adv.qlx[i,j,k] = buf.x
                adv.dqx[i,j,k] = buf.y
                adv.q6x[i,j,k] = buf.z
                buf = mf.ppm(my[i,j-2,k], my[i,j-1,k], my[i,j,k], my[i,j+1,k], my[i,j+2,k])
                adv.qly[i,j,k] = buf.x
                adv.dqy[i,j,k] = buf.y
                adv.q6y[i,j,k] = buf.z
    
    cursum = Float64(0.0)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds_no_pole,gm.full_jde_no_pole+1):
            for i in range(gm.half_ids,gm.half_ide+1):
                ci = math.trunc(adv.cflx[i,j,k])
                cf = adv.cflx[i,j,k] - ci
                if (adv.cflx[i,j,k] > 0.0):
                    s1 = 1.0 - cf
                    s2 = 1.0
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    cursum = 0.0
                    for l in range(i+1-ci,i+1):
                        cursum = cursum + mx[l,j,k]
                    mfx[i,j,k] = u[i,j,k] * (cursum + adv.qlx[i - ci,j,k] * ds1 + 0.5 * adv.dqx[i - ci,j,k] * ds2 + adv.q6x[i - ci,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflx[i,j,k]
                elif (adv.cflx[i,j,k] < 0.0):
                    s1 = 0.0
                    s2 = -cf
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    cursum = 0.0
                    for l in range(i+1,i-ci+1):
                        cursum = cursum + mx[l,j,k]
                    mfx[i,j,k] = -u[i,j,k] * (cursum + adv.qlx[i - ci + 1,j,k] * ds1 + 0.5 * adv.dqx[i - ci + 1,j,k] * ds2 + adv.q6x[i - ci + 1,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflx[i,j,k]
                else:
                    mfx[i,j,k] = 0.0
                
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.half_jds,gm.half_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                if (adv.cfly[i,j,k] > 0):
                    s1 = 1 - adv.cfly[i,j,k]
                    s2 = 1
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    mfy[i,j,k] =  v[i,j,k] * (adv.qly[i,j,k] * ds1 + 0.5 * adv.dqy[i,j,k] * ds2 + adv.q6y[i,j,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cfly[i,j,k]
                elif (adv.cfly[i,j,k] < 0):
                    s1 = 0
                    s2 = -adv.cfly[i,j,k]
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    mfy[i,j,k] = -v[i,j,k] * (adv.qly[i,j + 1,k] * ds1 + 0.5 * adv.dqy[i,j + 1,k] * ds2 + adv.q6y[i,j + 1,k] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cfly[i,j,k]
                else:
                    mfy[i,j,k] = 0
                    
@space_op
def vflx_ppm(adv: phv.HybridAdvField, w: LonLatField[Float64], m: LonLatField[Float64], mfz: LonLatField[Float64]):
    buf = Vector3(0.0, 0.0, 0.0)
    for k in range(gm.full_kds,gm.full_kde+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                buf = mf.ppm(m[i,j,k-2],m[i,j,k-1],m[i,j,k],m[i,j,k+1],m[i,j,k+2])
                adv.qlx[i,j,k] = buf.x
                adv.dqx[i,j,k] = buf.y
                adv.q6x[i,j,k] = buf.z
    
    cursum = Float64(0.0)
    for k in range(gm.half_kds+1,gm.half_kde-1+1):
        for j in range(gm.full_jds,gm.full_jde+1):
            for i in range(gm.full_ids,gm.full_ide+1):
                ci = math.trunc(adv.cflz[i,j,k])
                cf = adv.cflz[i,j,k] - ci
                if (adv.cflz[i,j,k] > 0):
                    s1 = 1.0 - cf
                    s2 = 1.0
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    cursum = 0.0
                    for l in range(k - ci,k - 1 + 1):
                        cursum = cursum + m[i,j,l]
                    mfz[i,j,k] =  w[i,j,k] * (cursum + adv.qlx[i,j,k - ci - 1] * ds1 + 0.5 * adv.dqx[i,j,k - ci - 1] * ds2 + adv.q6x[i,j,k - ci - 1] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflz[i,j,k]
                elif (adv.cflz[i,j,k] < 0):
                    s1 = 0.0
                    s2 = -cf
                    ds1 = s2    - s1
                    ds2 = s2**2 - s1**2
                    ds3 = s2**3 - s1**3
                    for l in range(k,k-ci-1+1):
                        cursum = cursum + m[i,j,l]
                    mfz[i,j,k] = -w[i,j,k] * (cursum + adv.qlx[i,j,k - ci] * ds1 + 0.5 * adv.dqx[i,j,k - ci] * ds2 + adv.q6x[i,j,k - ci] * (ds2 / 2.0 - ds3 / 3.0)) / adv.cflz[i,j,k]
                else:
                    mfz[i,j,k] = 0.0