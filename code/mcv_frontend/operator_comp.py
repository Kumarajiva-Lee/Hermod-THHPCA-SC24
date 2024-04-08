import math

from frontend.pytoy.lang import space_op, Float64, Int32, extern_func
from frontend.pytoy.lang.dtype_ext import LonLatField, CubedSphereField

import frontend.mcv_stroud.const_value as cv
import frontend.mcv_stroud.global_config as gc
import frontend.mcv_stroud.physical_variable as phv
import frontend.mcv_stroud.mesh as gm
import frontend.mcv_stroud.math_formual as mf
from frontend.mcv_stroud.utils import Vector2

@space_op
def StateInitial(state:phv.HybridStateField,  mesh:gm.HybridMeshField):
    pypanel = Int32(0)

    #x=1 y=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_1[i,0,0],mesh.y_1[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_1[i,0,0],mesh.y_1[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_11[i,j,0] = jab * h
            state.u_11[i,j,0] = uv.x
            state.v_11[i,j,0] = uv.y
            state.jab_11[i,j,0] = jab
            alpha=mesh.x_1[i,0,0]/cv.r
            beta=mesh.y_1[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_11[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_11[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_11[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_11[i,j,0],state.v_11[i,j,0],mesh.x_1[i,0,0],mesh.y_1[0,j,0])
            state.uc_11[i,j,0] = uvc.x
            state.vc_11[i,j,0] = uvc.y

    #x=2 y=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_2[i,0,0],mesh.y_1[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_2[i,0,0],mesh.y_1[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_21[i,j,0] = jab * h
            state.u_21[i,j,0] = uv.x
            state.v_21[i,j,0] = uv.y
            state.jab_21[i,j,0] = jab
            alpha=mesh.x_2[i,0,0]/cv.r
            beta=mesh.y_1[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_21[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_21[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_21[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_21[i,j,0],state.v_21[i,j,0],mesh.x_2[i,0,0],mesh.y_1[0,j,0])
            state.uc_21[i,j,0] = uvc.x
            state.vc_21[i,j,0] = uvc.y

    #x=3 y=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_3[i,0,0],mesh.y_1[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_3[i,0,0],mesh.y_1[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_31[i,j,0] = jab * h
            state.u_31[i,j,0] = uv.x
            state.v_31[i,j,0] = uv.y
            state.jab_31[i,j,0] = jab
            alpha=mesh.x_3[i,0,0]/cv.r
            beta=mesh.y_1[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_31[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_31[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_31[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_31[i,j,0],state.v_31[i,j,0],mesh.x_3[i,0,0],mesh.y_1[0,j,0])
            state.uc_31[i,j,0] = uvc.x
            state.vc_31[i,j,0] = uvc.y

    #x=1 y=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_1[i,0,0],mesh.y_2[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_1[i,0,0],mesh.y_2[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_12[i,j,0] = jab * h
            state.u_12[i,j,0] = uv.x
            state.v_12[i,j,0] = uv.y
            state.jab_12[i,j,0] = jab
            alpha=mesh.x_1[i,0,0]/cv.r
            beta=mesh.y_2[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_12[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_12[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_12[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_12[i,j,0],state.v_12[i,j,0],mesh.x_1[i,0,0],mesh.y_2[0,j,0])
            state.uc_12[i,j,0] = uvc.x
            state.vc_12[i,j,0] = uvc.y

    #x=2 y=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_2[i,0,0],mesh.y_2[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_2[i,0,0],mesh.y_2[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_22[i,j,0] = jab * h
            state.u_22[i,j,0] = uv.x
            state.v_22[i,j,0] = uv.y
            state.jab_22[i,j,0] = jab
            alpha=mesh.x_2[i,0,0]/cv.r
            beta=mesh.y_2[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_22[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_22[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_22[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_22[i,j,0],state.v_22[i,j,0],mesh.x_2[i,0,0],mesh.y_2[0,j,0])
            state.uc_22[i,j,0] = uvc.x
            state.vc_22[i,j,0] = uvc.y

    #x=3 y=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_3[i,0,0],mesh.y_2[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_3[i,0,0],mesh.y_2[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_32[i,j,0] = jab * h
            state.u_32[i,j,0] = uv.x
            state.v_32[i,j,0] = uv.y
            state.jab_32[i,j,0] = jab
            alpha=mesh.x_3[i,0,0]/cv.r
            beta=mesh.y_2[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_32[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_32[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_32[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_32[i,j,0],state.v_32[i,j,0],mesh.x_3[i,0,0],mesh.y_2[0,j,0])
            state.uc_32[i,j,0] = uvc.x
            state.vc_32[i,j,0] = uvc.y

    #x=1 y=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_1[i,0,0],mesh.y_3[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_1[i,0,0],mesh.y_3[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_13[i,j,0] = jab * h
            state.u_13[i,j,0] = uv.x
            state.v_13[i,j,0] = uv.y
            state.jab_13[i,j,0] = jab
            alpha=mesh.x_1[i,0,0]/cv.r
            beta=mesh.y_3[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_13[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_13[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_13[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_13[i,j,0],state.v_13[i,j,0],mesh.x_1[i,0,0],mesh.y_3[0,j,0])
            state.uc_13[i,j,0] = uvc.x
            state.vc_13[i,j,0] = uvc.y

    #x=2 y=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_2[i,0,0],mesh.y_3[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_2[i,0,0],mesh.y_3[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_23[i,j,0] = jab * h
            state.u_23[i,j,0] = uv.x
            state.v_23[i,j,0] = uv.y
            state.jab_23[i,j,0] = jab
            alpha=mesh.x_2[i,0,0]/cv.r
            beta=mesh.y_3[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_23[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_23[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_23[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_23[i,j,0],state.v_23[i,j,0],mesh.x_2[i,0,0],mesh.y_3[0,j,0])
            state.uc_23[i,j,0] = uvc.x
            state.vc_23[i,j,0] = uvc.y

    #x=3 y=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ltret = Vector2(0.0,0.0)
            ltret = mf.pprop2sp(pypanel,mesh.x_3[i,0,0],mesh.y_3[0,j,0])
            lamb = ltret.x
            theta = ltret.y
            h = mf.computeh(lamb, theta)
            jab = mf.computejab(mesh.x_3[i,0,0],mesh.y_3[0,j,0])
            uv = Vector2(0.0,0.0)
            uv = mf.computevel(pypanel,lamb, theta)
            state.h_33[i,j,0] = jab * h
            state.u_33[i,j,0] = uv.x
            state.v_33[i,j,0] = uv.y
            state.jab_33[i,j,0] = jab
            alpha=mesh.x_3[i,0,0]/cv.r
            beta=mesh.y_3[0,j,0]/cv.r
            rho = math.sqrt(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))
            state.jab5_33[i,j,0] = (1.0+math.tan(beta)**2)*(rho**2.0*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab6_33[i,j,0] = math.tan(alpha)*math.tan(beta)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            state.jab7_33[i,j,0] = (1.0+math.tan(alpha)**2.)*(rho**2*math.cos(alpha)**2*math.cos(beta)**2)
            uvc = Vector2(0.0,0.0)
            uvc = mf.cov2contrav(state.u_33[i,j,0],state.v_33[i,j,0],mesh.x_3[i,0,0],mesh.y_3[0,j,0])
            state.uc_33[i,j,0] = uvc.x
            state.vc_33[i,j,0] = uvc.y

@space_op
def updateState(state_old: phv.HybridStateField, state_new: phv.HybridStateField):
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_11[i,j,0] = state_old.h_11[i,j,0]
            state_new.u_11[i,j,0] = state_old.u_11[i,j,0]
            state_new.v_11[i,j,0] = state_old.v_11[i,j,0]
            state_new.jab_11[i,j,0] = state_old.jab_11[i,j,0]
            state_new.jab5_11[i,j,0] = state_old.jab5_11[i,j,0]
            state_new.jab6_11[i,j,0] = state_old.jab6_11[i,j,0]
            state_new.jab7_11[i,j,0] = state_old.jab7_11[i,j,0]
            state_new.uc_11[i,j,0] = state_old.uc_11[i,j,0]
            state_new.vc_11[i,j,0] = state_old.vc_11[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_12[i,j,0] = state_old.h_12[i,j,0]
            state_new.u_12[i,j,0] = state_old.u_12[i,j,0]
            state_new.v_12[i,j,0] = state_old.v_12[i,j,0]
            state_new.jab_12[i,j,0] = state_old.jab_12[i,j,0]
            state_new.jab5_12[i,j,0] = state_old.jab5_12[i,j,0]
            state_new.jab6_12[i,j,0] = state_old.jab6_12[i,j,0]
            state_new.jab7_12[i,j,0] = state_old.jab7_12[i,j,0]
            state_new.uc_12[i,j,0] = state_old.uc_12[i,j,0]
            state_new.vc_12[i,j,0] = state_old.vc_12[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_13[i,j,0] = state_old.h_13[i,j,0]
            state_new.u_13[i,j,0] = state_old.u_13[i,j,0]
            state_new.v_13[i,j,0] = state_old.v_13[i,j,0]
            state_new.jab_13[i,j,0] = state_old.jab_13[i,j,0]
            state_new.jab5_13[i,j,0] = state_old.jab5_13[i,j,0]
            state_new.jab6_13[i,j,0] = state_old.jab6_13[i,j,0]
            state_new.jab7_13[i,j,0] = state_old.jab7_13[i,j,0]
            state_new.uc_13[i,j,0] = state_old.uc_13[i,j,0]
            state_new.vc_13[i,j,0] = state_old.vc_13[i,j,0]

    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_21[i,j,0] = state_old.h_21[i,j,0]
            state_new.u_21[i,j,0] = state_old.u_21[i,j,0]
            state_new.v_21[i,j,0] = state_old.v_21[i,j,0]
            state_new.jab_21[i,j,0] = state_old.jab_21[i,j,0]
            state_new.jab5_21[i,j,0] = state_old.jab5_21[i,j,0]
            state_new.jab6_21[i,j,0] = state_old.jab6_21[i,j,0]
            state_new.jab7_21[i,j,0] = state_old.jab7_21[i,j,0]
            state_new.uc_21[i,j,0] = state_old.uc_21[i,j,0]
            state_new.vc_21[i,j,0] = state_old.vc_21[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_22[i,j,0] = state_old.h_22[i,j,0]
            state_new.u_22[i,j,0] = state_old.u_22[i,j,0]
            state_new.v_22[i,j,0] = state_old.v_22[i,j,0]
            state_new.jab_22[i,j,0] = state_old.jab_22[i,j,0]
            state_new.jab5_22[i,j,0] = state_old.jab5_22[i,j,0]
            state_new.jab6_22[i,j,0] = state_old.jab6_22[i,j,0]
            state_new.jab7_22[i,j,0] = state_old.jab7_22[i,j,0]
            state_new.uc_22[i,j,0] = state_old.uc_22[i,j,0]
            state_new.vc_22[i,j,0] = state_old.vc_22[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_23[i,j,0] = state_old.h_23[i,j,0]
            state_new.u_23[i,j,0] = state_old.u_23[i,j,0]
            state_new.v_23[i,j,0] = state_old.v_23[i,j,0]
            state_new.jab_23[i,j,0] = state_old.jab_23[i,j,0]
            state_new.jab5_23[i,j,0] = state_old.jab5_23[i,j,0]
            state_new.jab6_23[i,j,0] = state_old.jab6_23[i,j,0]
            state_new.jab7_23[i,j,0] = state_old.jab7_23[i,j,0]
            state_new.uc_23[i,j,0] = state_old.uc_23[i,j,0]
            state_new.vc_23[i,j,0] = state_old.vc_23[i,j,0]

    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_31[i,j,0] = state_old.h_31[i,j,0]
            state_new.u_31[i,j,0] = state_old.u_31[i,j,0]
            state_new.v_31[i,j,0] = state_old.v_31[i,j,0]
            state_new.jab_31[i,j,0] = state_old.jab_31[i,j,0]
            state_new.jab5_31[i,j,0] = state_old.jab5_31[i,j,0]
            state_new.jab6_31[i,j,0] = state_old.jab6_31[i,j,0]
            state_new.jab7_31[i,j,0] = state_old.jab7_31[i,j,0]
            state_new.uc_31[i,j,0] = state_old.uc_31[i,j,0]
            state_new.vc_31[i,j,0] = state_old.vc_31[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_32[i,j,0] = state_old.h_32[i,j,0]
            state_new.u_32[i,j,0] = state_old.u_32[i,j,0]
            state_new.v_32[i,j,0] = state_old.v_32[i,j,0]
            state_new.jab_32[i,j,0] = state_old.jab_32[i,j,0]
            state_new.jab5_32[i,j,0] = state_old.jab5_32[i,j,0]
            state_new.jab6_32[i,j,0] = state_old.jab6_32[i,j,0]
            state_new.jab7_32[i,j,0] = state_old.jab7_32[i,j,0]
            state_new.uc_32[i,j,0] = state_old.uc_32[i,j,0]
            state_new.vc_32[i,j,0] = state_old.vc_32[i,j,0]
    for j in range(1, gc.ny):
        for i in range(1,gc.nx):
            state_new.h_33[i,j,0] = state_old.h_33[i,j,0]
            state_new.u_33[i,j,0] = state_old.u_33[i,j,0]
            state_new.v_33[i,j,0] = state_old.v_33[i,j,0]
            state_new.jab_33[i,j,0] = state_old.jab_33[i,j,0]
            state_new.jab5_33[i,j,0] = state_old.jab5_33[i,j,0]
            state_new.jab6_33[i,j,0] = state_old.jab6_33[i,j,0]
            state_new.jab7_33[i,j,0] = state_old.jab7_33[i,j,0]
            state_new.uc_33[i,j,0] = state_old.uc_33[i,j,0]
            state_new.vc_33[i,j,0] = state_old.vc_33[i,j,0]

@space_op
def updateX(state: phv.HybridStateField, tend: phv.HybridTendField):
    #yp=1
    #计算f/q _1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_31[i-1,j,0])+ math.sqrt(state.jab5_31[i-1,j,0] *cv.gra * state.h_31[i-1,j,0] / state.jab_31[i-1,j,0])
            tend.fh_1[i,j,0] = (((state.uc_31[i-1,j,0] * state.h_31[i-1,j,0]) + (state.uc_11[i,j,0] * state.h_11[i,j,0])) + ulocal*(state.h_31[i-1,j,0]-state.h_11[i,j,0])) / 2
            tend.fu_1[i,j,0] = (((cv.gra*state.h_31[i-1,j,0] / state.jab_31[i-1,j,0] + 0.5*(state.u_31[i-1,j,0]*state.uc_31[i-1,j,0] + state.v_31[i-1,j,0]*state.vc_31[i-1,j,0])) + (cv.gra*state.h_11[i,j,0] / state.jab_11[i,j,0] + 0.5*(state.u_11[i,j,0]*state.uc_11[i,j,0] + state.v_11[i,j,0]*state.vc_11[i,j,0]))) + ulocal*(state.u_31[i-1,j,0]-state.u_11[i,j,0])) / 2
            tend.fv_1[i,j,0] = ((0.0 + 0.0) + ulocal*(state.v_31[i-1,j,0]-state.v_11[i,j,0])) / 2
        
            tend.qh_1[i,j,0] = (state.h_31[i-1,j,0] + state.h_11[i,j,0]) / 2
            tend.qu_1[i,j,0] = (state.u_31[i-1,j,0] + state.u_11[i,j,0]) / 2
            tend.qv_1[i,j,0] = (state.v_31[i-1,j,0] + state.v_11[i,j,0]) / 2

    #1st-order derivative

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_31[i-1,j,0])+ math.sqrt(state.jab5_31[i-1,j,0] *cv.gra * state.h_31[i-1,j,0] / state.jab_31[i-1,j,0])
            dql = state.h_11[i-1,j,0] - 4 * state.h_21[i-1,j,0] + 3 * state.h_31[i-1,j,0]
            dqr = -(3*state.h_11[i,j,0] - 4*state.h_21[i,j,0] + state.h_31[i,j,0])
            dfl = (state.uc_11[i-1,j,0] * state.h_11[i-1,j,0]) - 4*(state.uc_21[i-1,j,0] * state.h_21[i-1,j,0]) + 3*(state.uc_31[i-1,j,0] * state.h_31[i-1,j,0])
            dfr = -(3*(state.uc_11[i,j,0] * state.h_11[i,j,0]) - 4*(state.uc_21[i,j,0] * state.h_21[i,j,0]) + (state.uc_31[i,j,0] * state.h_31[i,j,0]))
    
            tend.fh_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_31[i-1,j,0])+ math.sqrt(state.jab5_31[i-1,j,0] *cv.gra * state.h_31[i-1,j,0] / state.jab_31[i-1,j,0])
            dql = state.u_11[i-1,j,0] - 4 * state.u_21[i-1,j,0] + 3 * state.u_31[i-1,j,0]
            dqr = -(3*state.u_11[i,j,0] - 4*state.u_21[i,j,0] + state.u_31[i,j,0])
            dfl = (cv.gra*state.h_11[i-1,j,0] / state.jab_11[i-1,j,0] + 0.5*(state.u_11[i-1,j,0]*state.uc_11[i-1,j,0] + state.v_11[i-1,j,0]*state.vc_11[i-1,j,0])) - 4*(cv.gra*state.h_21[i-1,j,0] / state.jab_21[i-1,j,0] + 0.5*(state.u_21[i-1,j,0]*state.uc_21[i-1,j,0] + state.v_21[i-1,j,0]*state.vc_21[i-1,j,0])) + 3*(cv.gra*state.h_31[i-1,j,0] / state.jab_31[i-1,j,0] + 0.5*(state.u_31[i-1,j,0]*state.uc_31[i-1,j,0] + state.v_31[i-1,j,0]*state.vc_31[i-1,j,0]))
            dfr = -(3*(cv.gra*state.h_11[i,j,0] / state.jab_11[i,j,0] + 0.5*(state.u_11[i,j,0]*state.uc_11[i,j,0] + state.v_11[i,j,0]*state.vc_11[i,j,0])) - 4*(cv.gra*state.h_21[i,j,0] / state.jab_21[i,j,0] + 0.5*(state.u_21[i,j,0]*state.uc_21[i,j,0] + state.v_21[i,j,0]*state.vc_21[i,j,0])) + (cv.gra*state.h_31[i,j,0] / state.jab_31[i,j,0] + 0.5*(state.u_31[i,j,0]*state.uc_31[i,j,0] + state.v_31[i,j,0]*state.vc_31[i,j,0])))
            
            tend.fu_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_31[i-1,j,0])+ math.sqrt(state.jab5_31[i-1,j,0] *cv.gra * state.h_31[i-1,j,0] / state.jab_31[i-1,j,0])
            dql = state.v_11[i-1,j,0] - 4 * state.v_21[i-1,j,0] + 3 * state.v_31[i-1,j,0]
            dqr = -(3*state.v_11[i,j,0] - 4*state.v_21[i,j,0] + state.v_31[i,j,0])
            dfl = 0.0
            dfr = 0.0
            
            tend.fv_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d x
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qv_2[i,j,0] / cv.dx
            sr = tend.qv_2[i+1,j,0] / cv.dx
            sc = 3.0 * (tend.qv_1[i+1,j,0] - tend.qv_1[i,j,0]) / cv.dx / 2.0 - tend.qv_2[i,j,0] / cv.dx / 4.0 - tend.qv_2[i+1,j,0] / cv.dx / 4.0

            tend.dxh_11[i,j,0] = -tend.fh_2[i,j,0] / cv.dx
            tend.dxh_21[i,j,0] = -3.0 * (tend.fh_1[i+1,j,0] - tend.fh_1[i,j,0]) /cv.dx / 2.0 + tend.fh_2[i,j,0] / cv.dx / 4.0 + tend.fh_2[i+1,j,0] / cv.dx / 4.0
            tend.dxh_31[i,j,0] = -tend.fh_2[i+1,j,0] / cv.dx

            tend.dxu_11[i,j,0] = -tend.fu_2[i,j,0] / cv.dx + state.vc_11[i,j,0] * sl
            tend.dxu_21[i,j,0] = 3.0*(tend.fu_1[i+1,j,0]-tend.fu_1[i,j,0])/cv.dx/2.0 + tend.fu_2[i,j,0]/cv.dx/4.0 + tend.fu_2[i+1,j,0]/cv.dx/4.0 + state.vc_21[i,j,0] * sc
            tend.dxu_31[i,j,0] = -tend.fu_2[i+1,j,0] / cv.dx + state.vc_31[i,j,0] * sr

            tend.dxv_11[i,j,0] = -tend.fv_2[i,j,0] / cv.dx - state.uc_11[i,j,0] * sl
            tend.dxv_21[i,j,0] = -3.0*(tend.fv_1[i+1,j,0] - tend.fv_1[i,j,0]) / cv.dx / 2.0 + tend.fv_2[i,j,0] / cv.dx / 4.0 + tend.fv_2[i+1,j,0] / cv.dx / 4.0 - state.uc_21[i,j,0] * sc
            tend.dxv_31[i,j,0] = -tend.fv_2[i+1,j,0] / cv.dx - state.uc_31[i,j,0] * sr

    #------------------------------------------------------------------------
    #yp=2
    #计算f/q _1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_32[i-1,j,0])+ math.sqrt(state.jab5_32[i-1,j,0] *cv.gra * state.h_32[i-1,j,0] / state.jab_32[i-1,j,0])
            tend.fh_1[i,j,0] = (((state.uc_32[i-1,j,0] * state.h_32[i-1,j,0]) + (state.uc_12[i,j,0] * state.h_12[i,j,0])) + ulocal*(state.h_32[i-1,j,0]-state.h_12[i,j,0])) / 2
            tend.fu_1[i,j,0] = (((cv.gra*state.h_32[i-1,j,0] / state.jab_32[i-1,j,0] + 0.5*(state.u_32[i-1,j,0]*state.uc_32[i-1,j,0] + state.v_32[i-1,j,0]*state.vc_32[i-1,j,0])) + (cv.gra*state.h_12[i,j,0] / state.jab_12[i,j,0] + 0.5*(state.u_12[i,j,0]*state.uc_12[i,j,0] + state.v_12[i,j,0]*state.vc_12[i,j,0]))) + ulocal*(state.u_32[i-1,j,0]-state.u_12[i,j,0])) / 2
            tend.fv_1[i,j,0] = ((0.0 + 0.0) + ulocal*(state.v_32[i-1,j,0]-state.v_12[i,j,0])) / 2
        
            tend.qh_1[i,j,0] = (state.h_32[i-1,j,0] + state.h_12[i,j,0]) / 2
            tend.qu_1[i,j,0] = (state.u_32[i-1,j,0] + state.u_12[i,j,0]) / 2
            tend.qv_1[i,j,0] = (state.v_32[i-1,j,0] + state.v_12[i,j,0]) / 2

    #1st-order derivative

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_32[i-1,j,0])+ math.sqrt(state.jab5_32[i-1,j,0] *cv.gra * state.h_32[i-1,j,0] / state.jab_32[i-1,j,0])
            dql = state.h_12[i-1,j,0] - 4 * state.h_22[i-1,j,0] + 3 * state.h_32[i-1,j,0]
            dqr = -(3*state.h_12[i,j,0] - 4*state.h_22[i,j,0] + state.h_32[i,j,0])
            dfl = (state.uc_12[i-1,j,0] * state.h_12[i-1,j,0]) - 4*(state.uc_22[i-1,j,0] * state.h_22[i-1,j,0]) + 3*(state.uc_32[i-1,j,0] * state.h_32[i-1,j,0])
            dfr = -(3*(state.uc_12[i,j,0] * state.h_12[i,j,0]) - 4*(state.uc_22[i,j,0] * state.h_22[i,j,0]) + (state.uc_32[i,j,0] * state.h_32[i,j,0]))
    
            tend.fh_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_32[i-1,j,0])+ math.sqrt(state.jab5_32[i-1,j,0] *cv.gra * state.h_32[i-1,j,0] / state.jab_32[i-1,j,0])
            dql = state.u_12[i-1,j,0] - 4 * state.u_22[i-1,j,0] + 3 * state.u_32[i-1,j,0]
            dqr = -(3*state.u_12[i,j,0] - 4*state.u_22[i,j,0] + state.u_32[i,j,0])
            dfl = (cv.gra*state.h_12[i-1,j,0] / state.jab_12[i-1,j,0] + 0.5*(state.u_12[i-1,j,0]*state.uc_12[i-1,j,0] + state.v_12[i-1,j,0]*state.vc_12[i-1,j,0])) - 4*(cv.gra*state.h_22[i-1,j,0] / state.jab_22[i-1,j,0] + 0.5*(state.u_22[i-1,j,0]*state.uc_22[i-1,j,0] + state.v_22[i-1,j,0]*state.vc_22[i-1,j,0])) + 3*(cv.gra*state.h_32[i-1,j,0] / state.jab_32[i-1,j,0] + 0.5*(state.u_32[i-1,j,0]*state.uc_32[i-1,j,0] + state.v_32[i-1,j,0]*state.vc_32[i-1,j,0]))
            dfr = -(3*(cv.gra*state.h_12[i,j,0] / state.jab_12[i,j,0] + 0.5*(state.u_12[i,j,0]*state.uc_12[i,j,0] + state.v_12[i,j,0]*state.vc_12[i,j,0])) - 4*(cv.gra*state.h_22[i,j,0] / state.jab_22[i,j,0] + 0.5*(state.u_22[i,j,0]*state.uc_22[i,j,0] + state.v_22[i,j,0]*state.vc_22[i,j,0])) + (cv.gra*state.h_32[i,j,0] / state.jab_32[i,j,0] + 0.5*(state.u_32[i,j,0]*state.uc_32[i,j,0] + state.v_32[i,j,0]*state.vc_32[i,j,0])))
            
            tend.fu_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_32[i-1,j,0])+ math.sqrt(state.jab5_32[i-1,j,0] *cv.gra * state.h_32[i-1,j,0] / state.jab_32[i-1,j,0])
            dql = state.v_12[i-1,j,0] - 4 * state.v_22[i-1,j,0] + 3 * state.v_32[i-1,j,0]
            dqr = -(3*state.v_12[i,j,0] - 4*state.v_22[i,j,0] + state.v_32[i,j,0])
            dfl = 0.0
            dfr = 0.0
            
            tend.fv_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d x
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qv_2[i,j,0] / cv.dx
            sr = tend.qv_2[i+1,j,0] / cv.dx
            sc = 3.0 * (tend.qv_1[i+1,j,0] - tend.qv_1[i,j,0]) / cv.dx / 2.0 - tend.qv_2[i,j,0] / cv.dx / 4.0 - tend.qv_2[i+1,j,0] / cv.dx / 4.0

            tend.dxh_12[i,j,0] = -tend.fh_2[i,j,0] / cv.dx
            tend.dxh_22[i,j,0] = -3.0 * (tend.fh_1[i+1,j,0] - tend.fh_1[i,j,0]) /cv.dx / 2.0 + tend.fh_2[i,j,0] / cv.dx / 4.0 + tend.fh_2[i+1,j,0] / cv.dx / 4.0
            tend.dxh_32[i,j,0] = -tend.fh_2[i+1,j,0] / cv.dx

            tend.dxu_12[i,j,0] = -tend.fu_2[i,j,0] / cv.dx + state.vc_12[i,j,0] * sl
            tend.dxu_22[i,j,0] = 3.0*(tend.fu_1[i+1,j,0]-tend.fu_1[i,j,0])/cv.dx/2.0 + tend.fu_2[i,j,0]/cv.dx/4.0 + tend.fu_2[i+1,j,0]/cv.dx/4.0 + state.vc_22[i,j,0] * sc
            tend.dxu_32[i,j,0] = -tend.fu_2[i+1,j,0] / cv.dx + state.vc_32[i,j,0] * sr

            tend.dxv_12[i,j,0] = -tend.fv_2[i,j,0] / cv.dx - state.uc_12[i,j,0] * sl
            tend.dxv_22[i,j,0] = -3.0*(tend.fv_1[i+1,j,0] - tend.fv_1[i,j,0]) / cv.dx / 2.0 + tend.fv_2[i,j,0] / cv.dx / 4.0 + tend.fv_2[i+1,j,0] / cv.dx / 4.0 - state.uc_22[i,j,0] * sc
            tend.dxv_32[i,j,0] = -tend.fv_2[i+1,j,0] / cv.dx - state.uc_32[i,j,0] * sr

    #---------------------------------------------------------------------------------------------------
    #yp=3
    #计算f/q _1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_33[i-1,j,0])+ math.sqrt(state.jab5_33[i-1,j,0] *cv.gra * state.h_33[i-1,j,0] / state.jab_33[i-1,j,0])
            tend.fh_1[i,j,0] = (((state.uc_33[i-1,j,0] * state.h_33[i-1,j,0]) + (state.uc_13[i,j,0] * state.h_13[i,j,0])) + ulocal*(state.h_33[i-1,j,0]-state.h_13[i,j,0])) / 2
            tend.fu_1[i,j,0] = (((cv.gra*state.h_33[i-1,j,0] / state.jab_33[i-1,j,0] + 0.5*(state.u_33[i-1,j,0]*state.uc_33[i-1,j,0] + state.v_33[i-1,j,0]*state.vc_33[i-1,j,0])) + (cv.gra*state.h_13[i,j,0] / state.jab_13[i,j,0] + 0.5*(state.u_13[i,j,0]*state.uc_13[i,j,0] + state.v_13[i,j,0]*state.vc_13[i,j,0]))) + ulocal*(state.u_33[i-1,j,0]-state.u_13[i,j,0])) / 2
            tend.fv_1[i,j,0] = ((0.0 + 0.0) + ulocal*(state.v_33[i-1,j,0]-state.v_13[i,j,0])) / 2
        
            tend.qh_1[i,j,0] = (state.h_33[i-1,j,0] + state.h_13[i,j,0]) / 2
            tend.qu_1[i,j,0] = (state.u_33[i-1,j,0] + state.u_13[i,j,0]) / 2
            tend.qv_1[i,j,0] = (state.v_33[i-1,j,0] + state.v_13[i,j,0]) / 2

    #1st-order derivative

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_33[i-1,j,0])+ math.sqrt(state.jab5_33[i-1,j,0] *cv.gra * state.h_33[i-1,j,0] / state.jab_33[i-1,j,0])
            dql = state.h_13[i-1,j,0] - 4 * state.h_23[i-1,j,0] + 3 * state.h_33[i-1,j,0]
            dqr = -(3*state.h_13[i,j,0] - 4*state.h_23[i,j,0] + state.h_33[i,j,0])
            dfl = (state.uc_13[i-1,j,0] * state.h_13[i-1,j,0]) - 4*(state.uc_23[i-1,j,0] * state.h_23[i-1,j,0]) + 3*(state.uc_33[i-1,j,0] * state.h_33[i-1,j,0])
            dfr = -(3*(state.uc_13[i,j,0] * state.h_13[i,j,0]) - 4*(state.uc_23[i,j,0] * state.h_23[i,j,0]) + (state.uc_33[i,j,0] * state.h_33[i,j,0]))
    
            tend.fh_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_33[i-1,j,0])+ math.sqrt(state.jab5_33[i-1,j,0] *cv.gra * state.h_33[i-1,j,0] / state.jab_33[i-1,j,0])
            dql = state.u_13[i-1,j,0] - 4 * state.u_23[i-1,j,0] + 3 * state.u_33[i-1,j,0]
            dqr = -(3*state.u_13[i,j,0] - 4*state.u_23[i,j,0] + state.u_33[i,j,0])
            dfl = (cv.gra*state.h_13[i-1,j,0] / state.jab_13[i-1,j,0] + 0.5*(state.u_13[i-1,j,0]*state.uc_13[i-1,j,0] + state.v_13[i-1,j,0]*state.vc_13[i-1,j,0])) - 4*(cv.gra*state.h_23[i-1,j,0] / state.jab_23[i-1,j,0] + 0.5*(state.u_23[i-1,j,0]*state.uc_23[i-1,j,0] + state.v_23[i-1,j,0]*state.vc_23[i-1,j,0])) + 3*(cv.gra*state.h_33[i-1,j,0] / state.jab_33[i-1,j,0] + 0.5*(state.u_33[i-1,j,0]*state.uc_33[i-1,j,0] + state.v_33[i-1,j,0]*state.vc_33[i-1,j,0]))
            dfr = -(3*(cv.gra*state.h_13[i,j,0] / state.jab_13[i,j,0] + 0.5*(state.u_13[i,j,0]*state.uc_13[i,j,0] + state.v_13[i,j,0]*state.vc_13[i,j,0])) - 4*(cv.gra*state.h_23[i,j,0] / state.jab_23[i,j,0] + 0.5*(state.u_23[i,j,0]*state.uc_23[i,j,0] + state.v_23[i,j,0]*state.vc_23[i,j,0])) + (cv.gra*state.h_33[i,j,0] / state.jab_33[i,j,0] + 0.5*(state.u_33[i,j,0]*state.uc_33[i,j,0] + state.v_33[i,j,0]*state.vc_33[i,j,0])))
            
            tend.fu_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ulocal = math.fabs(state.uc_33[i-1,j,0])+ math.sqrt(state.jab5_33[i-1,j,0] *cv.gra * state.h_33[i-1,j,0] / state.jab_33[i-1,j,0])
            dql = state.v_13[i-1,j,0] - 4 * state.v_23[i-1,j,0] + 3 * state.v_33[i-1,j,0]
            dqr = -(3*state.v_13[i,j,0] - 4*state.v_23[i,j,0] + state.v_33[i,j,0])
            dfl = 0.0
            dfr = 0.0
            
            tend.fv_2[i,j,0] = (dfl + dfr) + ulocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d x
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qv_2[i,j,0] / cv.dx
            sr = tend.qv_2[i+1,j,0] / cv.dx
            sc = 3.0 * (tend.qv_1[i+1,j,0] - tend.qv_1[i,j,0]) / cv.dx / 2.0 - tend.qv_2[i,j,0] / cv.dx / 4.0 - tend.qv_2[i+1,j,0] / cv.dx / 4.0

            tend.dxh_13[i,j,0] = -tend.fh_2[i,j,0] / cv.dx
            tend.dxh_23[i,j,0] = -3.0 * (tend.fh_1[i+1,j,0] - tend.fh_1[i,j,0]) /cv.dx / 2.0 + tend.fh_2[i,j,0] / cv.dx / 4.0 + tend.fh_2[i+1,j,0] / cv.dx / 4.0
            tend.dxh_33[i,j,0] = -tend.fh_2[i+1,j,0] / cv.dx

            tend.dxu_13[i,j,0] = -tend.fu_2[i,j,0] / cv.dx + state.vc_13[i,j,0] * sl
            tend.dxu_23[i,j,0] = 3.0*(tend.fu_1[i+1,j,0]-tend.fu_1[i,j,0])/cv.dx/2.0 + tend.fu_2[i,j,0]/cv.dx/4.0 + tend.fu_2[i+1,j,0]/cv.dx/4.0 + state.vc_23[i,j,0] * sc
            tend.dxu_33[i,j,0] = -tend.fu_2[i+1,j,0] / cv.dx + state.vc_33[i,j,0] * sr

            tend.dxv_13[i,j,0] = -tend.fv_2[i,j,0] / cv.dx - state.uc_13[i,j,0] * sl
            tend.dxv_23[i,j,0] = -3.0*(tend.fv_1[i+1,j,0] - tend.fv_1[i,j,0]) / cv.dx / 2.0 + tend.fv_2[i,j,0] / cv.dx / 4.0 + tend.fv_2[i+1,j,0] / cv.dx / 4.0 - state.uc_23[i,j,0] * sc
            tend.dxv_33[i,j,0] = -tend.fv_2[i+1,j,0] / cv.dx - state.uc_33[i,j,0] * sr

@space_op
def updateY(state: phv.HybridStateField, tend: phv.HybridTendField):
    #ToRemember y下标没改
    
    #xp=1
    #计算f/q _1  
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_13[i,j-1,0])+ math.sqrt(state.jab7_13[i,j-1,0] *cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0])
            tend.fh_1[i,j,0] = ((state.vc_13[i,j-1,0] * state.h_13[i,j-1,0]) + (state.vc_11[i,j,0] * state.h_11[i,j,0])) + vlocal * (state.h_13[i,j-1,0] - state.h_11[i,j,0] )/2.0
            tend.fu_1[i,j,0] = (0.0 + 0.0) + vlocal * (state.u_13[i,j-1,0] - state.u_11[i,j,0]) / 2.0
            tend.fv_1[i,j,0] = ((cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0] + 0.5*(state.u_13[i,j-1,0]*state.uc_13[i,j-1,0] + state.v_13[i,j-1,0]*state.vc_13[i,j-1,0])) + (cv.gra * state.h_11[i,j,0] / state.jab_11[i,j,0] + 0.5*(state.u_11[i,j,0]*state.uc_11[i,j,0] + state.v_11[i,j,0]*state.vc_11[i,j,0]))) + vlocal * (state.v_13[i,j-1,0] - state.v_11[i,j,0]) / 2.0

            tend.qh_1[i,j,0] = (state.h_13[i,j-1,0] + state.h_11[i,j,0]) / 2.0
            tend.qu_1[i,j,0] = (state.u_13[i,j-1,0] + state.u_11[i,j,0]) / 2.0
            tend.qv_1[i,j,0] = (state.v_13[i,j-1,0] + state.v_11[i,j,0]) / 2.0

    #1st-order derivative           

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_13[i,j-1,0])+ math.sqrt(state.jab7_13[i,j-1,0] *cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0])
            dql = state.h_11[i,j-1,0] - 4 * state.h_12[i,j-1,0] + 3 * state.h_13[i,j-1,0]
            dqr = - (3*state.h_11[i,j,0] - 4 * state.h_12[i,j,0] + state.h_13[i,j,0])
            dfl = (state.vc_11[i,j-1,0] * state.h_11[i,j-1,0]) - 4*(state.vc_12[i,j-1,0] * state.h_12[i,j-1,0]) + 3*(state.vc_13[i,j-1,0] * state.h_13[i,j-1,0])
            dfr = -(3*(state.vc_11[i,j,0] * state.h_11[i,j,0]) - 4*(state.vc_12[i,j,0] * state.h_12[i,j,0]) + (state.vc_13[i,j,0] * state.h_13[i,j,0]) )

            tend.fh_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_13[i,j-1,0])+ math.sqrt(state.jab7_13[i,j-1,0] *cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0])
            dql = state.u_11[i,j-1,0] - 4 * state.u_12[i,j-1,0] + 3 * state.u_13[i,j-1,0]
            dqr = - (3*state.u_11[i,j,0] - 4 * state.u_12[i,j,0] + state.u_13[i,j,0])
            dfl = 0.0
            dfr = 0.0

            tend.fu_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_13[i,j-1,0])+ math.sqrt(state.jab7_13[i,j-1,0] *cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0])
            dql = state.v_11[i,j-1,0] - 4 * state.v_12[i,j-1,0] + 3 * state.v_13[i,j-1,0]
            dqr = - (3*state.v_11[i,j,0] - 4 * state.v_12[i,j,0] + state.v_13[i,j,0])
            dfl = (cv.gra * state.h_11[i,j-1,0] / state.jab_11[i,j-1,0] + 0.5*(state.u_11[i,j-1,0]*state.uc_11[i,j-1,0] + state.v_11[i,j-1,0]*state.vc_11[i,j-1,0])) - 4*(cv.gra * state.h_12[i,j-1,0] / state.jab_12[i,j-1,0] + 0.5*(state.u_12[i,j-1,0]*state.uc_12[i,j-1,0] + state.v_12[i,j-1,0]*state.vc_12[i,j-1,0])) + 3*(cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0] + 0.5*(state.u_13[i,j-1,0]*state.uc_13[i,j-1,0] + state.v_13[i,j-1,0]*state.vc_13[i,j-1,0]))
            dfr = -(3*(cv.gra * state.h_11[i,j,0] / state.jab_11[i,j,0] + 0.5*(state.u_11[i,j,0]*state.uc_11[i,j,0] + state.v_11[i,j,0]*state.vc_11[i,j,0])) - 4*(cv.gra * state.h_12[i,j,0] / state.jab_12[i,j,0] + 0.5*(state.u_12[i,j,0]*state.uc_12[i,j,0] + state.v_12[i,j,0]*state.vc_12[i,j,0])) + (cv.gra * state.h_13[i,j,0] / state.jab_13[i,j,0] + 0.5*(state.u_13[i,j,0]*state.uc_13[i,j,0] + state.v_13[i,j,0]*state.vc_13[i,j,0])) )

            tend.fv_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d y
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qu_2[i,j,0] / cv.dy
            sr = tend.qu_2[i,j+1,0] / cv.dy
            sc = 3.0*(tend.qu_1[i,j+1,0] - tend.qu_1[i,j,0]) /cv.dy / 2.0 - tend.qu_2[i,j,0] / cv.dy / 4.0 - tend.qu_2[i,j+1,0] / cv.dy / 4.0

            tend.dyh_11[i,j,0] = -tend.fh_2[i,j,0] / cv.dy
            tend.dyh_12[i,j,0] = -3*(tend.fh_1[i,j+1,0] - tend.fh_1[i,j,0]) / cv.dy / 2.0 + tend.fh_2[i,j,0] / cv.dy / 4.0 + tend.fh_2[i,j+1,0] / cv.dy / 4.0
            tend.dyh_13[i,j,0] = -tend.fh_2[i,j+1,0] / cv.dy

            tend.dyu_11[i,j,0] = -tend.fu_2[i,j,0] / cv.dy - state.vc_11[i,j,0] * sl
            tend.dyu_12[i,j,0] = -3.0*(tend.fu_1[i,j+1,0] - tend.fu_1[i,j,0]) / cv.dy / 2.0 + tend.fu_2[i,j,0] / cv.dy / 4.0 + tend.fu_2[i,j+1,0] / cv.dy / 4.0 - state.vc_12[i,j,0] * sc
            tend.dyu_13[i,j,0] = -tend.fu_2[i,j+1,0] / cv.dy - state.vc_13[i,j,0] * sr

            tend.dyv_11[i,j,0] = -tend.fv_2[i,j,0] / cv.dy + state.uc_11[i,j,0] * sl
            tend.dyv_12[i,j,0] = -3.0*(tend.fv_1[i,j+1,0] - tend.fv_1[i,j,0]) / cv.dy / 2.0 + tend.fv_2[i,j,0] / cv.dy / 4.0 + tend.fv_2[i,j+1,0] / cv.dy / 4.0 - state.uc_12[i,j,0] * sc
            tend.dyv_13[i,j,0] = -tend.fv_2[i,j+1,0] / cv.dy - state.uc_13[i,j,0] * sr

    #----------------------------------------------------------------------
    #xp=2
    #计算f/q _1  
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_23[i,j-1,0])+ math.sqrt(state.jab7_23[i,j-1,0] *cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0])
            tend.fh_1[i,j,0] = ((state.vc_23[i,j-1,0] * state.h_23[i,j-1,0]) + (state.vc_21[i,j,0] * state.h_21[i,j,0])) + vlocal * (state.h_23[i,j-1,0] - state.h_21[i,j,0] )/2.0
            tend.fu_1[i,j,0] = (0.0 + 0.0) + vlocal * (state.u_23[i,j-1,0] - state.u_21[i,j,0]) / 2.0
            tend.fv_1[i,j,0] = ((cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0] + 0.5*(state.u_23[i,j-1,0]*state.uc_23[i,j-1,0] + state.v_23[i,j-1,0]*state.vc_23[i,j-1,0])) + (cv.gra * state.h_21[i,j,0] / state.jab_21[i,j,0] + 0.5*(state.u_21[i,j,0]*state.uc_21[i,j,0] + state.v_21[i,j,0]*state.vc_21[i,j,0]))) + vlocal * (state.v_23[i,j-1,0] - state.v_21[i,j,0]) / 2.0

            tend.qh_1[i,j,0] = (state.h_23[i,j-1,0] + state.h_21[i,j,0]) / 2.0
            tend.qu_1[i,j,0] = (state.u_23[i,j-1,0] + state.u_21[i,j,0]) / 2.0
            tend.qv_1[i,j,0] = (state.v_23[i,j-1,0] + state.v_21[i,j,0]) / 2.0

    #1st-order derivative           

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_23[i,j-1,0])+ math.sqrt(state.jab7_23[i,j-1,0] *cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0])
            dql = state.h_21[i,j-1,0] - 4 * state.h_22[i,j-1,0] + 3 * state.h_23[i,j-1,0]
            dqr = - (3*state.h_21[i,j,0] - 4 * state.h_22[i,j,0] + state.h_23[i,j,0])
            dfl = (state.vc_21[i,j-1,0] * state.h_21[i,j-1,0]) - 4*(state.vc_22[i,j-1,0] * state.h_22[i,j-1,0]) + 3*(state.vc_23[i,j-1,0] * state.h_23[i,j-1,0])
            dfr = -(3*(state.vc_21[i,j,0] * state.h_21[i,j,0]) - 4*(state.vc_22[i,j,0] * state.h_22[i,j,0]) + (state.vc_23[i,j,0] * state.h_23[i,j,0]) )

            tend.fh_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_23[i,j-1,0])+ math.sqrt(state.jab7_23[i,j-1,0] *cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0])
            dql = state.u_21[i,j-1,0] - 4 * state.u_22[i,j-1,0] + 3 * state.u_23[i,j-1,0]
            dqr = - (3*state.u_21[i,j,0] - 4 * state.u_22[i,j,0] + state.u_23[i,j,0])
            dfl = 0.0
            dfr = 0.0

            tend.fu_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_23[i,j-1,0])+ math.sqrt(state.jab7_23[i,j-1,0] *cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0])
            dql = state.v_21[i,j-1,0] - 4 * state.v_22[i,j-1,0] + 3 * state.v_23[i,j-1,0]
            dqr = - (3*state.v_21[i,j,0] - 4 * state.v_22[i,j,0] + state.v_23[i,j,0])
            dfl = (cv.gra * state.h_21[i,j-1,0] / state.jab_21[i,j-1,0] + 0.5*(state.u_21[i,j-1,0]*state.uc_21[i,j-1,0] + state.v_21[i,j-1,0]*state.vc_21[i,j-1,0])) - 4*(cv.gra * state.h_22[i,j-1,0] / state.jab_22[i,j-1,0] + 0.5*(state.u_22[i,j-1,0]*state.uc_22[i,j-1,0] + state.v_22[i,j-1,0]*state.vc_22[i,j-1,0])) + 3*(cv.gra * state.h_23[i,j-1,0] / state.jab_23[i,j-1,0] + 0.5*(state.u_23[i,j-1,0]*state.uc_23[i,j-1,0] + state.v_23[i,j-1,0]*state.vc_23[i,j-1,0]))
            dfr = -(3*(cv.gra * state.h_21[i,j,0] / state.jab_21[i,j,0] + 0.5*(state.u_21[i,j,0]*state.uc_21[i,j,0] + state.v_21[i,j,0]*state.vc_21[i,j,0])) - 4*(cv.gra * state.h_22[i,j,0] / state.jab_22[i,j,0] + 0.5*(state.u_22[i,j,0]*state.uc_22[i,j,0] + state.v_22[i,j,0]*state.vc_22[i,j,0])) + (cv.gra * state.h_23[i,j,0] / state.jab_23[i,j,0] + 0.5*(state.u_23[i,j,0]*state.uc_23[i,j,0] + state.v_23[i,j,0]*state.vc_23[i,j,0])) )

            tend.fv_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d y
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qu_2[i,j,0] / cv.dy
            sr = tend.qu_2[i,j+1,0] / cv.dy
            sc = 3.0*(tend.qu_1[i,j+1,0] - tend.qu_1[i,j,0]) /cv.dy / 2.0 - tend.qu_2[i,j,0] / cv.dy / 4.0 - tend.qu_2[i,j+1,0] / cv.dy / 4.0

            tend.dyh_11[i,j,0] = -tend.fh_2[i,j,0] / cv.dy
            tend.dyh_12[i,j,0] = -3*(tend.fh_1[i,j+1,0] - tend.fh_1[i,j,0]) / cv.dy / 2.0 + tend.fh_2[i,j,0] / cv.dy / 4.0 + tend.fh_2[i,j+1,0] / cv.dy / 4.0
            tend.dyh_13[i,j,0] = -tend.fh_2[i,j+1,0] / cv.dy

            tend.dyu_11[i,j,0] = -tend.fu_2[i,j,0] / cv.dy - state.vc_21[i,j,0] * sl
            tend.dyu_12[i,j,0] = -3.0*(tend.fu_1[i,j+1,0] - tend.fu_1[i,j,0]) / cv.dy / 2.0 + tend.fu_2[i,j,0] / cv.dy / 4.0 + tend.fu_2[i,j+1,0] / cv.dy / 4.0 - state.vc_22[i,j,0] * sc
            tend.dyu_13[i,j,0] = -tend.fu_2[i,j+1,0] / cv.dy - state.vc_23[i,j,0] * sr

            tend.dyv_11[i,j,0] = -tend.fv_2[i,j,0] / cv.dy + state.uc_21[i,j,0] * sl
            tend.dyv_12[i,j,0] = -3.0*(tend.fv_1[i,j+1,0] - tend.fv_1[i,j,0]) / cv.dy / 2.0 + tend.fv_2[i,j,0] / cv.dy / 4.0 + tend.fv_2[i,j+1,0] / cv.dy / 4.0 - state.uc_22[i,j,0] * sc
            tend.dyv_13[i,j,0] = -tend.fv_2[i,j+1,0] / cv.dy - state.uc_23[i,j,0] * sr

    #----------------------------------------------------------------------
    #xp=3
    #计算f/q _1  
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_33[i,j-1,0])+ math.sqrt(state.jab7_33[i,j-1,0] *cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0])
            tend.fh_1[i,j,0] = ((state.vc_33[i,j-1,0] * state.h_33[i,j-1,0]) + (state.vc_31[i,j,0] * state.h_31[i,j,0])) + vlocal * (state.h_33[i,j-1,0] - state.h_31[i,j,0] )/2.0
            tend.fu_1[i,j,0] = (0.0 + 0.0) + vlocal * (state.u_33[i,j-1,0] - state.u_31[i,j,0]) / 2.0
            tend.fv_1[i,j,0] = ((cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0] + 0.5*(state.u_33[i,j-1,0]*state.uc_33[i,j-1,0] + state.v_33[i,j-1,0]*state.vc_33[i,j-1,0])) + (cv.gra * state.h_31[i,j,0] / state.jab_31[i,j,0] + 0.5*(state.u_31[i,j,0]*state.uc_31[i,j,0] + state.v_31[i,j,0]*state.vc_31[i,j,0]))) + vlocal * (state.v_33[i,j-1,0] - state.v_31[i,j,0]) / 2.0

            tend.qh_1[i,j,0] = (state.h_33[i,j-1,0] + state.h_31[i,j,0]) / 2.0
            tend.qu_1[i,j,0] = (state.u_33[i,j-1,0] + state.u_31[i,j,0]) / 2.0
            tend.qv_1[i,j,0] = (state.v_33[i,j-1,0] + state.v_31[i,j,0]) / 2.0

    #1st-order derivative           

    #计算f/q _2 h
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_33[i,j-1,0])+ math.sqrt(state.jab7_33[i,j-1,0] *cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0])
            dql = state.h_31[i,j-1,0] - 4 * state.h_32[i,j-1,0] + 3 * state.h_33[i,j-1,0]
            dqr = - (3*state.h_31[i,j,0] - 4 * state.h_32[i,j,0] + state.h_33[i,j,0])
            dfl = (state.vc_31[i,j-1,0] * state.h_31[i,j-1,0]) - 4*(state.vc_32[i,j-1,0] * state.h_32[i,j-1,0]) + 3*(state.vc_33[i,j-1,0] * state.h_33[i,j-1,0])
            dfr = -(3*(state.vc_31[i,j,0] * state.h_31[i,j,0]) - 4*(state.vc_32[i,j,0] * state.h_32[i,j,0]) + (state.vc_33[i,j,0] * state.h_33[i,j,0]) )

            tend.fh_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qh_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 u
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_33[i,j-1,0])+ math.sqrt(state.jab7_33[i,j-1,0] *cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0])
            dql = state.u_31[i,j-1,0] - 4 * state.u_32[i,j-1,0] + 3 * state.u_33[i,j-1,0]
            dqr = - (3*state.u_31[i,j,0] - 4 * state.u_32[i,j,0] + state.u_33[i,j,0])
            dfl = 0.0
            dfr = 0.0

            tend.fu_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qu_2[i,j,0] = (dql + dqr) / 2

    #计算f/q _2 v
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            vlocal = math.fabs(state.vc_33[i,j-1,0])+ math.sqrt(state.jab7_33[i,j-1,0] *cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0])
            dql = state.v_31[i,j-1,0] - 4 * state.v_32[i,j-1,0] + 3 * state.v_33[i,j-1,0]
            dqr = - (3*state.v_31[i,j,0] - 4 * state.v_32[i,j,0] + state.v_33[i,j,0])
            dfl = (cv.gra * state.h_31[i,j-1,0] / state.jab_31[i,j-1,0] + 0.5*(state.u_31[i,j-1,0]*state.uc_31[i,j-1,0] + state.v_31[i,j-1,0]*state.vc_31[i,j-1,0])) - 4*(cv.gra * state.h_32[i,j-1,0] / state.jab_32[i,j-1,0] + 0.5*(state.u_32[i,j-1,0]*state.uc_32[i,j-1,0] + state.v_32[i,j-1,0]*state.vc_32[i,j-1,0])) + 3*(cv.gra * state.h_33[i,j-1,0] / state.jab_33[i,j-1,0] + 0.5*(state.u_33[i,j-1,0]*state.uc_33[i,j-1,0] + state.v_33[i,j-1,0]*state.vc_33[i,j-1,0]))
            dfr = -(3*(cv.gra * state.h_31[i,j,0] / state.jab_31[i,j,0] + 0.5*(state.u_31[i,j,0]*state.uc_31[i,j,0] + state.v_31[i,j,0]*state.vc_31[i,j,0])) - 4*(cv.gra * state.h_32[i,j,0] / state.jab_32[i,j,0] + 0.5*(state.u_32[i,j,0]*state.uc_32[i,j,0] + state.v_32[i,j,0]*state.vc_32[i,j,0])) + (cv.gra * state.h_33[i,j,0] / state.jab_33[i,j,0] + 0.5*(state.u_33[i,j,0]*state.uc_33[i,j,0] + state.v_33[i,j,0]*state.vc_33[i,j,0])) )

            tend.fv_2[i,j,0] = (dfl + dfr) + vlocal*(dql - dqr) / 2
            tend.qv_2[i,j,0] = (dql + dqr) / 2

    #计算d y
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            sl = tend.qu_2[i,j,0] / cv.dy
            sr = tend.qu_2[i,j+1,0] / cv.dy
            sc = 3.0*(tend.qu_1[i,j+1,0] - tend.qu_1[i,j,0]) /cv.dy / 2.0 - tend.qu_2[i,j,0] / cv.dy / 4.0 - tend.qu_2[i,j+1,0] / cv.dy / 4.0

            tend.dyh_11[i,j,0] = -tend.fh_2[i,j,0] / cv.dy
            tend.dyh_12[i,j,0] = -3*(tend.fh_1[i,j+1,0] - tend.fh_1[i,j,0]) / cv.dy / 2.0 + tend.fh_2[i,j,0] / cv.dy / 4.0 + tend.fh_2[i,j+1,0] / cv.dy / 4.0
            tend.dyh_13[i,j,0] = -tend.fh_2[i,j+1,0] / cv.dy

            tend.dyu_11[i,j,0] = -tend.fu_2[i,j,0] / cv.dy - state.vc_31[i,j,0] * sl
            tend.dyu_12[i,j,0] = -3.0*(tend.fu_1[i,j+1,0] - tend.fu_1[i,j,0]) / cv.dy / 2.0 + tend.fu_2[i,j,0] / cv.dy / 4.0 + tend.fu_2[i,j+1,0] / cv.dy / 4.0 - state.vc_32[i,j,0] * sc
            tend.dyu_13[i,j,0] = -tend.fu_2[i,j+1,0] / cv.dy - state.vc_33[i,j,0] * sr

            tend.dyv_11[i,j,0] = -tend.fv_2[i,j,0] / cv.dy + state.uc_31[i,j,0] * sl
            tend.dyv_12[i,j,0] = -3.0*(tend.fv_1[i,j+1,0] - tend.fv_1[i,j,0]) / cv.dy / 2.0 + tend.fv_2[i,j,0] / cv.dy / 4.0 + tend.fv_2[i,j+1,0] / cv.dy / 4.0 - state.uc_32[i,j,0] * sc
            tend.dyv_13[i,j,0] = -tend.fv_2[i,j+1,0] / cv.dy - state.uc_33[i,j,0] * sr

@space_op
def updateUV(state: phv.HybridStateField, tend: phv.HybridTendField, mesh: gm.HybridMeshField):
    ret = Vector2(0.0,0.0)
    pypanel = Int32(0)
    #i=1 j=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_1[i,0,0], mesh.y_1[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_11[i,j,0] = tend.dxv_11[i,j,0] - state.uc_11[i,j,0] * state.jab_11[i,j,0] * f
            tend.dyu_11[i,j,0] = tend.dyu_11[i,j,0] - state.vc_11[i,j,0] * state.jab_11[i,j,0] * f
    
    #i=1 j=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_1[i,0,0], mesh.y_2[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_12[i,j,0] = tend.dxv_12[i,j,0] - state.uc_12[i,j,0] * state.jab_12[i,j,0] * f
            tend.dyu_12[i,j,0] = tend.dyu_12[i,j,0] - state.vc_12[i,j,0] * state.jab_12[i,j,0] * f

    #i=1 j=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_1[i,0,0], mesh.y_3[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_13[i,j,0] = tend.dxv_13[i,j,0] - state.uc_13[i,j,0] * state.jab_13[i,j,0] * f
            tend.dyu_13[i,j,0] = tend.dyu_13[i,j,0] - state.vc_13[i,j,0] * state.jab_13[i,j,0] * f

    #i=2 j=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_2[i,0,0], mesh.y_1[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_21[i,j,0] = tend.dxv_21[i,j,0] - state.uc_21[i,j,0] * state.jab_21[i,j,0] * f
            tend.dyu_21[i,j,0] = tend.dyu_21[i,j,0] - state.vc_21[i,j,0] * state.jab_21[i,j,0] * f
    
    #i=2 j=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_2[i,0,0], mesh.y_2[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_22[i,j,0] = tend.dxv_22[i,j,0] - state.uc_22[i,j,0] * state.jab_22[i,j,0] * f
            tend.dyu_22[i,j,0] = tend.dyu_22[i,j,0] - state.vc_22[i,j,0] * state.jab_22[i,j,0] * f

    #i=2 j=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_2[i,0,0], mesh.y_3[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_23[i,j,0] = tend.dxv_23[i,j,0] - state.uc_23[i,j,0] * state.jab_23[i,j,0] * f
            tend.dyu_23[i,j,0] = tend.dyu_23[i,j,0] - state.vc_23[i,j,0] * state.jab_23[i,j,0] * f

    #i=3 j=1
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_3[i,0,0], mesh.y_1[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_31[i,j,0] = tend.dxv_31[i,j,0] - state.uc_31[i,j,0] * state.jab_31[i,j,0] * f
            tend.dyu_31[i,j,0] = tend.dyu_31[i,j,0] - state.vc_31[i,j,0] * state.jab_31[i,j,0] * f
    
    #i=3 j=2
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_3[i,0,0], mesh.y_2[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_32[i,j,0] = tend.dxv_32[i,j,0] - state.uc_32[i,j,0] * state.jab_32[i,j,0] * f
            tend.dyu_32[i,j,0] = tend.dyu_32[i,j,0] - state.vc_32[i,j,0] * state.jab_32[i,j,0] * f

    #i=3 j=3
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            ret = mf.pprop2sp(pypanel, mesh.x_3[i,0,0], mesh.y_3[0,j,0])
            lamb = ret.x
            theta = ret.y
            f = 2.0 * cv.omega * (- math.cos(lamb) * math.cos(theta) * math.sin(0.0)+ math.sin(theta) * math.cos(0.0))
            tend.dxv_33[i,j,0] = tend.dxv_33[i,j,0] - state.uc_33[i,j,0] * state.jab_33[i,j,0] * f
            tend.dyu_33[i,j,0] = tend.dyu_33[i,j,0] - state.vc_33[i,j,0] * state.jab_33[i,j,0] * f

@space_op
def update_rk1(state_old: phv.HybridStateField, state_new: phv.HybridStateField, tend1: phv.HybridTendField, dt: Float64):
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_11[i,j,0] = state_old.h_11[i,j,0] + (tend1.dxh_11[i,j,0] + tend1.dyh_11[i,j,0]) * dt
            state_new.u_11[i,j,0] = state_old.u_11[i,j,0] + (tend1.dxu_11[i,j,0] + tend1.dyu_11[i,j,0]) * dt
            state_new.v_11[i,j,0] = state_old.v_11[i,j,0] + (tend1.dxv_11[i,j,0] + tend1.dyv_11[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_12[i,j,0] = state_old.h_12[i,j,0] + (tend1.dxh_12[i,j,0] + tend1.dyh_12[i,j,0]) * dt
            state_new.u_12[i,j,0] = state_old.u_12[i,j,0] + (tend1.dxu_12[i,j,0] + tend1.dyu_12[i,j,0]) * dt
            state_new.v_12[i,j,0] = state_old.v_12[i,j,0] + (tend1.dxv_12[i,j,0] + tend1.dyv_12[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_13[i,j,0] = state_old.h_13[i,j,0] + (tend1.dxh_13[i,j,0] + tend1.dyh_13[i,j,0]) * dt
            state_new.u_13[i,j,0] = state_old.u_13[i,j,0] + (tend1.dxu_13[i,j,0] + tend1.dyu_13[i,j,0]) * dt
            state_new.v_13[i,j,0] = state_old.v_13[i,j,0] + (tend1.dxv_13[i,j,0] + tend1.dyv_13[i,j,0]) * dt

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_21[i,j,0] = state_old.h_21[i,j,0] + (tend1.dxh_21[i,j,0] + tend1.dyh_21[i,j,0]) * dt
            state_new.u_21[i,j,0] = state_old.u_21[i,j,0] + (tend1.dxu_21[i,j,0] + tend1.dyu_21[i,j,0]) * dt
            state_new.v_21[i,j,0] = state_old.v_21[i,j,0] + (tend1.dxv_21[i,j,0] + tend1.dyv_21[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_22[i,j,0] = state_old.h_22[i,j,0] + (tend1.dxh_22[i,j,0] + tend1.dyh_22[i,j,0]) * dt
            state_new.u_22[i,j,0] = state_old.u_22[i,j,0] + (tend1.dxu_22[i,j,0] + tend1.dyu_22[i,j,0]) * dt
            state_new.v_22[i,j,0] = state_old.v_22[i,j,0] + (tend1.dxv_22[i,j,0] + tend1.dyv_22[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_23[i,j,0] = state_old.h_23[i,j,0] + (tend1.dxh_23[i,j,0] + tend1.dyh_23[i,j,0]) * dt
            state_new.u_23[i,j,0] = state_old.u_23[i,j,0] + (tend1.dxu_23[i,j,0] + tend1.dyu_23[i,j,0]) * dt
            state_new.v_23[i,j,0] = state_old.v_23[i,j,0] + (tend1.dxv_23[i,j,0] + tend1.dyv_23[i,j,0]) * dt

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_31[i,j,0] = state_old.h_31[i,j,0] + (tend1.dxh_31[i,j,0] + tend1.dyh_31[i,j,0]) * dt
            state_new.u_31[i,j,0] = state_old.u_31[i,j,0] + (tend1.dxu_31[i,j,0] + tend1.dyu_31[i,j,0]) * dt
            state_new.v_31[i,j,0] = state_old.v_31[i,j,0] + (tend1.dxv_31[i,j,0] + tend1.dyv_31[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_32[i,j,0] = state_old.h_32[i,j,0] + (tend1.dxh_32[i,j,0] + tend1.dyh_32[i,j,0]) * dt
            state_new.u_32[i,j,0] = state_old.u_32[i,j,0] + (tend1.dxu_32[i,j,0] + tend1.dyu_32[i,j,0]) * dt
            state_new.v_32[i,j,0] = state_old.v_32[i,j,0] + (tend1.dxv_32[i,j,0] + tend1.dyv_32[i,j,0]) * dt
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_33[i,j,0] = state_old.h_33[i,j,0] + (tend1.dxh_33[i,j,0] + tend1.dyh_33[i,j,0]) * dt
            state_new.u_33[i,j,0] = state_old.u_33[i,j,0] + (tend1.dxu_33[i,j,0] + tend1.dyu_33[i,j,0]) * dt
            state_new.v_33[i,j,0] = state_old.v_33[i,j,0] + (tend1.dxv_33[i,j,0] + tend1.dyv_33[i,j,0]) * dt

@space_op
def update_rk2(state_old: phv.HybridStateField, state_new: phv.HybridStateField, tend1: phv.HybridTendField, tend2: phv.HybridTendField, dt: Float64):
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_11[i,j,0] = state_old.h_11[i,j,0] + (tend1.dxh_11[i,j,0] + tend1.dyh_11[i,j,0] + tend2.dxh_11[i,j,0] + tend2.dyh_11[i,j,0]) * dt / 4.0
            state_new.u_11[i,j,0] = state_old.u_11[i,j,0] + (tend1.dxu_11[i,j,0] + tend1.dyu_11[i,j,0] + tend2.dxu_11[i,j,0] + tend2.dyu_11[i,j,0]) * dt / 4.0
            state_new.v_11[i,j,0] = state_old.v_11[i,j,0] + (tend1.dxv_11[i,j,0] + tend1.dyv_11[i,j,0] + tend2.dxv_11[i,j,0] + tend2.dyv_11[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_12[i,j,0] = state_old.h_12[i,j,0] + (tend1.dxh_12[i,j,0] + tend1.dyh_12[i,j,0] + tend2.dxh_12[i,j,0] + tend2.dyh_12[i,j,0]) * dt / 4.0
            state_new.u_12[i,j,0] = state_old.u_12[i,j,0] + (tend1.dxu_12[i,j,0] + tend1.dyu_12[i,j,0] + tend2.dxu_12[i,j,0] + tend2.dyu_12[i,j,0]) * dt / 4.0
            state_new.v_12[i,j,0] = state_old.v_12[i,j,0] + (tend1.dxv_12[i,j,0] + tend1.dyv_12[i,j,0] + tend2.dxv_12[i,j,0] + tend2.dyv_12[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_13[i,j,0] = state_old.h_13[i,j,0] + (tend1.dxh_13[i,j,0] + tend1.dyh_13[i,j,0] + tend2.dxh_13[i,j,0] + tend2.dyh_13[i,j,0]) * dt / 4.0
            state_new.u_13[i,j,0] = state_old.u_13[i,j,0] + (tend1.dxu_13[i,j,0] + tend1.dyu_13[i,j,0] + tend2.dxu_13[i,j,0] + tend2.dyu_13[i,j,0]) * dt / 4.0
            state_new.v_13[i,j,0] = state_old.v_13[i,j,0] + (tend1.dxv_13[i,j,0] + tend1.dyv_13[i,j,0] + tend2.dxv_13[i,j,0] + tend2.dyv_13[i,j,0]) * dt / 4.0

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_21[i,j,0] = state_old.h_21[i,j,0] + (tend1.dxh_21[i,j,0] + tend1.dyh_21[i,j,0] + tend2.dxh_21[i,j,0] + tend2.dyh_21[i,j,0]) * dt / 4.0
            state_new.u_21[i,j,0] = state_old.u_21[i,j,0] + (tend1.dxu_21[i,j,0] + tend1.dyu_21[i,j,0] + tend2.dxu_21[i,j,0] + tend2.dyu_21[i,j,0]) * dt / 4.0
            state_new.v_21[i,j,0] = state_old.v_21[i,j,0] + (tend1.dxv_21[i,j,0] + tend1.dyv_21[i,j,0] + tend2.dxv_21[i,j,0] + tend2.dyv_21[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_22[i,j,0] = state_old.h_22[i,j,0] + (tend1.dxh_22[i,j,0] + tend1.dyh_22[i,j,0] + tend2.dxh_22[i,j,0] + tend2.dyh_22[i,j,0]) * dt / 4.0
            state_new.u_22[i,j,0] = state_old.u_22[i,j,0] + (tend1.dxu_22[i,j,0] + tend1.dyu_22[i,j,0] + tend2.dxu_22[i,j,0] + tend2.dyu_22[i,j,0]) * dt / 4.0
            state_new.v_22[i,j,0] = state_old.v_22[i,j,0] + (tend1.dxv_22[i,j,0] + tend1.dyv_22[i,j,0] + tend2.dxv_22[i,j,0] + tend2.dyv_22[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_23[i,j,0] = state_old.h_23[i,j,0] + (tend1.dxh_23[i,j,0] + tend1.dyh_23[i,j,0] + tend2.dxh_23[i,j,0] + tend2.dyh_23[i,j,0]) * dt / 4.0
            state_new.u_23[i,j,0] = state_old.u_23[i,j,0] + (tend1.dxu_23[i,j,0] + tend1.dyu_23[i,j,0] + tend2.dxu_23[i,j,0] + tend2.dyu_23[i,j,0]) * dt / 4.0
            state_new.v_23[i,j,0] = state_old.v_23[i,j,0] + (tend1.dxv_23[i,j,0] + tend1.dyv_23[i,j,0] + tend2.dxv_23[i,j,0] + tend2.dyv_23[i,j,0]) * dt / 4.0

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_31[i,j,0] = state_old.h_31[i,j,0] + (tend1.dxh_31[i,j,0] + tend1.dyh_31[i,j,0] + tend2.dxh_31[i,j,0] + tend2.dyh_31[i,j,0]) * dt / 4.0
            state_new.u_31[i,j,0] = state_old.u_31[i,j,0] + (tend1.dxu_31[i,j,0] + tend1.dyu_31[i,j,0] + tend2.dxu_31[i,j,0] + tend2.dyu_31[i,j,0]) * dt / 4.0
            state_new.v_31[i,j,0] = state_old.v_31[i,j,0] + (tend1.dxv_31[i,j,0] + tend1.dyv_31[i,j,0] + tend2.dxv_31[i,j,0] + tend2.dyv_31[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_32[i,j,0] = state_old.h_32[i,j,0] + (tend1.dxh_32[i,j,0] + tend1.dyh_32[i,j,0] + tend2.dxh_32[i,j,0] + tend2.dyh_32[i,j,0]) * dt / 4.0
            state_new.u_32[i,j,0] = state_old.u_32[i,j,0] + (tend1.dxu_32[i,j,0] + tend1.dyu_32[i,j,0] + tend2.dxu_32[i,j,0] + tend2.dyu_32[i,j,0]) * dt / 4.0
            state_new.v_32[i,j,0] = state_old.v_32[i,j,0] + (tend1.dxv_32[i,j,0] + tend1.dyv_32[i,j,0] + tend2.dxv_32[i,j,0] + tend2.dyv_32[i,j,0]) * dt / 4.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_33[i,j,0] = state_old.h_33[i,j,0] + (tend1.dxh_33[i,j,0] + tend1.dyh_33[i,j,0] + tend2.dxh_33[i,j,0] + tend2.dyh_33[i,j,0]) * dt / 4.0
            state_new.u_33[i,j,0] = state_old.u_33[i,j,0] + (tend1.dxu_33[i,j,0] + tend1.dyu_33[i,j,0] + tend2.dxu_33[i,j,0] + tend2.dyu_33[i,j,0]) * dt / 4.0
            state_new.v_33[i,j,0] = state_old.v_33[i,j,0] + (tend1.dxv_33[i,j,0] + tend1.dyv_33[i,j,0] + tend2.dxv_33[i,j,0] + tend2.dyv_33[i,j,0]) * dt / 4.0

@space_op
def update_rk3(state_old: phv.HybridStateField, state_new: phv.HybridStateField, tend1: phv.HybridTendField, tend2: phv.HybridTendField, tend3: phv.HybridTendField, dt: Float64):
       
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_11[i,j,0] = state_old.h_11[i,j,0] + (tend1.dxh_11[i,j,0] + tend1.dyh_11[i,j,0] + tend2.dxh_11[i,j,0] + tend2.dyh_11[i,j,0] + 4.0 * tend3.dxh_11[i,j,0] + 4.0 * tend3.dyh_11[i,j,0]) * dt / 6.0
            state_new.u_11[i,j,0] = state_old.u_11[i,j,0] + (tend1.dxu_11[i,j,0] + tend1.dyu_11[i,j,0] + tend2.dxu_11[i,j,0] + tend2.dyu_11[i,j,0] + 4.0 * tend3.dxu_11[i,j,0] + 4.0 * tend3.dyu_11[i,j,0]) * dt / 6.0
            state_new.v_11[i,j,0] = state_old.v_11[i,j,0] + (tend1.dxv_11[i,j,0] + tend1.dyv_11[i,j,0] + tend2.dxv_11[i,j,0] + tend2.dyv_11[i,j,0] + 4.0 * tend3.dxv_11[i,j,0] + 4.0 * tend3.dyv_11[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_12[i,j,0] = state_old.h_12[i,j,0] + (tend1.dxh_12[i,j,0] + tend1.dyh_12[i,j,0] + tend2.dxh_12[i,j,0] + tend2.dyh_12[i,j,0] + 4.0 * tend3.dxh_12[i,j,0] + 4.0 * tend3.dyh_12[i,j,0]) * dt / 6.0
            state_new.u_12[i,j,0] = state_old.u_12[i,j,0] + (tend1.dxu_12[i,j,0] + tend1.dyu_12[i,j,0] + tend2.dxu_12[i,j,0] + tend2.dyu_12[i,j,0] + 4.0 * tend3.dxu_12[i,j,0] + 4.0 * tend3.dyu_12[i,j,0]) * dt / 6.0
            state_new.v_12[i,j,0] = state_old.v_12[i,j,0] + (tend1.dxv_12[i,j,0] + tend1.dyv_12[i,j,0] + tend2.dxv_12[i,j,0] + tend2.dyv_12[i,j,0] + 4.0 * tend3.dxv_12[i,j,0] + 4.0 * tend3.dyv_12[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_13[i,j,0] = state_old.h_13[i,j,0] + (tend1.dxh_13[i,j,0] + tend1.dyh_13[i,j,0] + tend2.dxh_13[i,j,0] + tend2.dyh_13[i,j,0] + 4.0 * tend3.dxh_13[i,j,0] + 4.0 * tend3.dyh_13[i,j,0]) * dt / 6.0
            state_new.u_13[i,j,0] = state_old.u_13[i,j,0] + (tend1.dxu_13[i,j,0] + tend1.dyu_13[i,j,0] + tend2.dxu_13[i,j,0] + tend2.dyu_13[i,j,0] + 4.0 * tend3.dxu_13[i,j,0] + 4.0 * tend3.dyu_13[i,j,0]) * dt / 6.0
            state_new.v_13[i,j,0] = state_old.v_13[i,j,0] + (tend1.dxv_13[i,j,0] + tend1.dyv_13[i,j,0] + tend2.dxv_13[i,j,0] + tend2.dyv_13[i,j,0] + 4.0 * tend3.dxv_13[i,j,0] + 4.0 * tend3.dyv_13[i,j,0]) * dt / 6.0

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_21[i,j,0] = state_old.h_21[i,j,0] + (tend1.dxh_21[i,j,0] + tend1.dyh_21[i,j,0] + tend2.dxh_21[i,j,0] + tend2.dyh_21[i,j,0] + 4.0 * tend3.dxh_21[i,j,0] + 4.0 * tend3.dyh_21[i,j,0]) * dt / 6.0
            state_new.u_21[i,j,0] = state_old.u_21[i,j,0] + (tend1.dxu_21[i,j,0] + tend1.dyu_21[i,j,0] + tend2.dxu_21[i,j,0] + tend2.dyu_21[i,j,0] + 4.0 * tend3.dxu_21[i,j,0] + 4.0 * tend3.dyu_21[i,j,0]) * dt / 6.0
            state_new.v_21[i,j,0] = state_old.v_21[i,j,0] + (tend1.dxv_21[i,j,0] + tend1.dyv_21[i,j,0] + tend2.dxv_21[i,j,0] + tend2.dyv_21[i,j,0] + 4.0 * tend3.dxv_21[i,j,0] + 4.0 * tend3.dyv_21[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_22[i,j,0] = state_old.h_22[i,j,0] + (tend1.dxh_22[i,j,0] + tend1.dyh_22[i,j,0] + tend2.dxh_22[i,j,0] + tend2.dyh_22[i,j,0] + 4.0 * tend3.dxh_22[i,j,0] + 4.0 * tend3.dyh_22[i,j,0]) * dt / 6.0
            state_new.u_22[i,j,0] = state_old.u_22[i,j,0] + (tend1.dxu_22[i,j,0] + tend1.dyu_22[i,j,0] + tend2.dxu_22[i,j,0] + tend2.dyu_22[i,j,0] + 4.0 * tend3.dxu_22[i,j,0] + 4.0 * tend3.dyu_22[i,j,0]) * dt / 6.0
            state_new.v_22[i,j,0] = state_old.v_22[i,j,0] + (tend1.dxv_22[i,j,0] + tend1.dyv_22[i,j,0] + tend2.dxv_22[i,j,0] + tend2.dyv_22[i,j,0] + 4.0 * tend3.dxv_22[i,j,0] + 4.0 * tend3.dyv_22[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_23[i,j,0] = state_old.h_23[i,j,0] + (tend1.dxh_23[i,j,0] + tend1.dyh_23[i,j,0] + tend2.dxh_23[i,j,0] + tend2.dyh_23[i,j,0] + 4.0 * tend3.dxh_23[i,j,0] + 4.0 * tend3.dyh_23[i,j,0]) * dt / 6.0
            state_new.u_23[i,j,0] = state_old.u_23[i,j,0] + (tend1.dxu_23[i,j,0] + tend1.dyu_23[i,j,0] + tend2.dxu_23[i,j,0] + tend2.dyu_23[i,j,0] + 4.0 * tend3.dxu_23[i,j,0] + 4.0 * tend3.dyu_23[i,j,0]) * dt / 6.0
            state_new.v_23[i,j,0] = state_old.v_23[i,j,0] + (tend1.dxv_23[i,j,0] + tend1.dyv_23[i,j,0] + tend2.dxv_23[i,j,0] + tend2.dyv_23[i,j,0] + 4.0 * tend3.dxv_23[i,j,0] + 4.0 * tend3.dyv_23[i,j,0]) * dt / 6.0

    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_31[i,j,0] = state_old.h_31[i,j,0] + (tend1.dxh_31[i,j,0] + tend1.dyh_31[i,j,0] + tend2.dxh_31[i,j,0] + tend2.dyh_31[i,j,0] + 4.0 * tend3.dxh_31[i,j,0] + 4.0 * tend3.dyh_31[i,j,0]) * dt / 6.0
            state_new.u_31[i,j,0] = state_old.u_31[i,j,0] + (tend1.dxu_31[i,j,0] + tend1.dyu_31[i,j,0] + tend2.dxu_31[i,j,0] + tend2.dyu_31[i,j,0] + 4.0 * tend3.dxu_31[i,j,0] + 4.0 * tend3.dyu_31[i,j,0]) * dt / 6.0
            state_new.v_31[i,j,0] = state_old.v_31[i,j,0] + (tend1.dxv_31[i,j,0] + tend1.dyv_31[i,j,0] + tend2.dxv_31[i,j,0] + tend2.dyv_31[i,j,0] + 4.0 * tend3.dxv_31[i,j,0] + 4.0 * tend3.dyv_31[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_32[i,j,0] = state_old.h_32[i,j,0] + (tend1.dxh_32[i,j,0] + tend1.dyh_32[i,j,0] + tend2.dxh_32[i,j,0] + tend2.dyh_32[i,j,0] + 4.0 * tend3.dxh_32[i,j,0] + 4.0 * tend3.dyh_32[i,j,0]) * dt / 6.0
            state_new.u_32[i,j,0] = state_old.u_32[i,j,0] + (tend1.dxu_32[i,j,0] + tend1.dyu_32[i,j,0] + tend2.dxu_32[i,j,0] + tend2.dyu_32[i,j,0] + 4.0 * tend3.dxu_32[i,j,0] + 4.0 * tend3.dyu_32[i,j,0]) * dt / 6.0
            state_new.v_32[i,j,0] = state_old.v_32[i,j,0] + (tend1.dxv_32[i,j,0] + tend1.dyv_32[i,j,0] + tend2.dxv_32[i,j,0] + tend2.dyv_32[i,j,0] + 4.0 * tend3.dxv_32[i,j,0] + 4.0 * tend3.dyv_32[i,j,0]) * dt / 6.0
    for j in range(0, gc.ny):
        for i in range(0,gc.nx):
            state_new.h_33[i,j,0] = state_old.h_33[i,j,0] + (tend1.dxh_33[i,j,0] + tend1.dyh_33[i,j,0] + tend2.dxh_33[i,j,0] + tend2.dyh_33[i,j,0] + 4.0 * tend3.dxh_33[i,j,0] + 4.0 * tend3.dyh_33[i,j,0]) * dt / 6.0
            state_new.u_33[i,j,0] = state_old.u_33[i,j,0] + (tend1.dxu_33[i,j,0] + tend1.dyu_33[i,j,0] + tend2.dxu_33[i,j,0] + tend2.dyu_33[i,j,0] + 4.0 * tend3.dxu_33[i,j,0] + 4.0 * tend3.dyu_33[i,j,0]) * dt / 6.0
            state_new.v_33[i,j,0] = state_old.v_33[i,j,0] + (tend1.dxv_33[i,j,0] + tend1.dyv_33[i,j,0] + tend2.dxv_33[i,j,0] + tend2.dyv_33[i,j,0] + 4.0 * tend3.dxv_33[i,j,0] + 4.0 * tend3.dyv_33[i,j,0]) * dt / 6.0

    # fl11 = state.vc_11[i,j-1,0] * state.h_11[i,j-1,0]
    # fl12 = 0.0
    # fl13 = cv.gra * state.h_11[i,j-1,0] / state.jab_11[i,j-1,0] + 0.5*(state.u_11[i,j-1,0]*state.uc_11[i,j-1,0] + state.v_11[i,j-1,0]*state.vc_11[i,j-1,0])
    # fr11 = state.vc_11[i,j,0] * state.h_11[i,j,0]
    # fr12 = 0.0
    # fr13 = cv.gra * state.h_11[i,j,0] / state.jab_11[i,j,0] + 0.5*(state.u_11[i,j,0]*state.uc_11[i,j,0] + state.v_11[i,j,0]*state.vc_11[i,j,0])

    # fl21 = state.vc_12[i,j-1,0] * state.h_12[i,j-1,0]
    # fl22 = 0.0
    # fl23 = cv.gra * state.h_12[i,j-1,0] / state.jab_12[i,j-1,0] + 0.5*(state.u_12[i,j-1,0]*state.uc_12[i,j-1,0] + state.v_12[i,j-1,0]*state.vc_12[i,j-1,0])
    # fr21 = state.vc_12[i,j,0] * state.h_12[i,j,0]
    # fr22 = 0.0
    # fr23 = cv.gra * state.h_12[i,j,0] / state.jab_12[i,j,0] + 0.5*(state.u_12[i,j,0]*state.uc_12[i,j,0] + state.v_12[i,j,0]*state.vc_12[i,j,0])

    # fl31 = state.vc_13[i,j-1,0] * state.h_13[i,j-1,0]
    # fl32 = 0.0
    # fl33 = cv.gra * state.h_13[i,j-1,0] / state.jab_13[i,j-1,0] + 0.5*(state.u_13[i,j-1,0]*state.uc_13[i,j-1,0] + state.v_13[i,j-1,0]*state.vc_13[i,j-1,0])
    # fr31 = state.vc_13[i,j,0] * state.h_13[i,j,0]
    # fr32 = 0.0
    # fr33 = cv.gra * state.h_13[i,j,0] / state.jab_13[i,j,0] + 0.5*(state.u_13[i,j,0]*state.uc_13[i,j,0] + state.v_13[i,j,0]*state.vc_13[i,j,0])