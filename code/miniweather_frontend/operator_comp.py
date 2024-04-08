import math

from frontend.pytoy.lang import space_op, Float64, Int32, extern_func
from frontend.pytoy.lang.dtype_ext import LonLatField

import miniweather_scuderia.physical_variable as phv
import miniweather_scuderia.const_value as cv
import miniweather_scuderia.global_config as gc
import miniweather_scuderia.mesh as gm
import miniweather_scuderia.math_formual as mf
from miniweather_scuderia.utils import InitVector, Vector4

@space_op
def StateInitial(state:phv.HybridStateField, state_tmp:phv.HybridStateField, staticv:phv.HybridStaticField, mesh:gm.HybridMeshField, DirPara:phv.HybridDirPara):
    
    #Initialize the cell-averaged fluid state via Gauss-Legendre quadrature
    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx):
                #Skip Initialize the state to zero,it has done after allocate
                for kk in range(0, cv.npoints):
                    for ii in range(0, cv.npoints):
                        #No Local Array. Crying~
                        if (ii == 0): 
                            qpi = cv.qpoints_1
                            qwi = cv.qweights_1
                        elif (ii == 1):
                            qpi = cv.qpoints_2
                            qwi = cv.qweights_2
                        elif (ii == 2):
                            qpi = cv.qpoints_3
                            qwi = cv.qweights_3
                        if (kk == 0): 
                            qpk = cv.qpoints_1
                            qwk = cv.qweights_1
                        elif (kk == 1):
                            qpk = cv.qpoints_2
                            qwk = cv.qweights_2
                        elif (kk == 2):
                            qpk = cv.qpoints_3
                            qwk = cv.qweights_3

                        x = mesh.xp[i,0,0]  + (qpi-0.5) * gc.dx
                        z = mesh.zpf[0,0,k] + (qpk-0.5) * gc.dz

                        InitV = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        if (gc.initial == 0):
                            InitV = mf.collision(x,z)
                        if (gc.initial == 1):
                            InitV = mf.thermal(x,z)
                        if (gc.initial == 2):
                            InitV = mf.gravity_waves(x,z)

                        state.rho[i,j,k] = state.rho[i,j,k] + InitV.r * qwi * qwk
                        state.u[i,j,k]   = state.u[i,j,k]   + (InitV.r+InitV.hr) * InitV.u * qwi * qwk
                        state.w[i,j,k]   = state.w[i,j,k]   + (InitV.r+InitV.hr) * InitV.w * qwi * qwk
                        state.pt[i,j,k]  = state.pt[i,j,k]  + ((InitV.r+InitV.hr) * (InitV.t+InitV.ht) - InitV.hr*InitV.ht) * qwi * qwk              

    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state_tmp.rho[i,j,k] = state.rho[i,j,k]
                state_tmp.u[i,j,k] = state.u[i,j,k]
                state_tmp.w[i,j,k] = state.w[i,j,k]
                state_tmp.pt[i,j,k] = state.pt[i,j,k]
    
    #Compute the hydrostatic background state over vertical cell averages
    for k in range(-2,gc.nz+2):
        for kk in range(0, cv.npoints):
            if (kk == 0): 
                qpk = cv.qpoints_1
                qwk = cv.qweights_1
            elif (kk == 1):
                qpk = cv.qpoints_2
                qwk = cv.qweights_2
            elif (kk == 2):
                qpk = cv.qpoints_3
                qwk = cv.qweights_3
            
            z = mesh.zpf[0,0,k]
            InitV = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

            if (gc.initial == 0):
                InitV = mf.collision(0.0,z)
            if (gc.initial == 1):
                InitV = mf.thermal(0.0,z)
            if (gc.initial == 2):
                InitV = mf.gravity_waves(0.0,z)
            staticv.hy_dens_cell[0,0,k]       = staticv.hy_dens_cell[0,0,k]       + InitV.hr * qwk
            staticv.hy_dens_theta_cell[0,0,k] = staticv.hy_dens_theta_cell[0,0,k] + InitV.hr * InitV.ht * qwk
    
    #Compute the hydrostatic background state at vertical cell interfaces
    for k in range(0,gc.nz+1):
        z = mesh.zph[0,0,k]
        InitV = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if (gc.initial == 0):
            InitV = mf.collision(0.0,z)
        if (gc.initial == 1):
            InitV = mf.thermal(0.0,z)
        if (gc.initial == 2):
            InitV = mf.gravity_waves(0.0,z)
        staticv.hy_dens_int[0,0,k] = InitV.hr
        staticv.hy_dens_theta_int[0,0,k] = InitV.hr * InitV.ht
        staticv.hy_pressure_int[0,0,k] = cv.c0 * (InitV.hr * InitV.ht)**cv.gamm

    DirPara.xd = True
    DirPara.zd = False
    DirPara.stepnum = 0

@space_op
def SwitchDirection(DirPara:phv.HybridDirPara):
    if (DirPara.xd):
        DirPara.xd = False
        DirPara.zd = True
    else:
        DirPara.xd = True
        DirPara.zd = False 
    DirPara.stepnum = 0

@space_op
def ComputeTendenciesX(state:phv.HybridStateField, staticv:phv.HybridStaticField, flux:phv.HybridFluxField, tend:phv.HybridTendField, dt:Float64):
    #Compute the hyperviscosity coefficient
    hv_coef = -cv.hv_beta * gc.dx / (16*dt)    

    vals = Vector4(0.0,0.0,0.0,0.0)
    d3_vals = Vector4(0.0,0.0,0.0,0.0)

    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx+1):
                #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
                vals.r = -state.rho[i-2,j,k]/12 + 7*state.rho[i-1,j,k]/12 + 7*state.rho[i,j,k]/12 - state.rho[i+1,j,k]/12
                d3_vals.r = -state.rho[i-2,j,k] + 3*state.rho[i-1,j,k] - 3*state.rho[i,j,k] + state.rho[i+1,j,k]
                vals.u = -state.u[i-2,j,k]/12 + 7*state.u[i-1,j,k]/12 + 7*state.u[i,j,k]/12 - state.u[i+1,j,k]/12
                d3_vals.u = -state.u[i-2,j,k] + 3*state.u[i-1,j,k] - 3*state.u[i,j,k] + state.u[i+1,j,k]
                vals.w = -state.w[i-2,j,k]/12 + 7*state.w[i-1,j,k]/12 + 7*state.w[i,j,k]/12 - state.w[i+1,j,k]/12
                d3_vals.w = -state.w[i-2,j,k] + 3*state.w[i-1,j,k] - 3*state.w[i,j,k] + state.w[i+1,j,k]
                vals.t = -state.pt[i-2,j,k]/12 + 7*state.pt[i-1,j,k]/12 + 7*state.pt[i,j,k]/12 - state.pt[i+1,j,k]/12
                d3_vals.t = -state.pt[i-2,j,k] + 3*state.pt[i-1,j,k] - 3*state.pt[i,j,k] + state.pt[i+1,j,k]

                #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
                r = vals.r + staticv.hy_dens_cell[0,0,k]
                u = vals.u / r
                w = vals.w / r
                t = (vals.t + staticv.hy_dens_theta_cell[0,0,k]) / r
                p = cv.c0 * ((r*t) ** cv.gamm)

                #Compute the flux vector
                flux.frho[i,j,k] = r * u - hv_coef * d3_vals.r
                flux.fu[i,j,k] = r * u * u + p - hv_coef * d3_vals.u
                flux.fw[i,j,k] = r * u * w - hv_coef * d3_vals.w
                flux.fpt[i,j,k] = r * u * t - hv_coef * d3_vals.t

    #Use the fluxes to compute tendencies for each cell
    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx):
                tend.drho[i,j,k] = -(flux.frho[i+1,j,k] - flux.frho[i,j,k]) / gc.dx
                tend.du[i,j,k] = -(flux.fu[i+1,j,k] - flux.fu[i,j,k]) / gc.dx
                tend.dw[i,j,k] = -(flux.fw[i+1,j,k] - flux.fw[i,j,k]) / gc.dx
                tend.dpt[i,j,k] = -(flux.fpt[i+1,j,k] - flux.fpt[i,j,k]) / gc.dx

@space_op
def SetBoundaryZ(state:phv.HybridStateField, staticv:phv.HybridStaticField):
    for k in range(-2,-1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.rho[i,j,k] = state.rho[i,j,k+2]
    for k in range(-1,0):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.rho[i,j,k] = state.rho[i,j,k+1]
    for k in range(gc.nz,gc.nz+1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.rho[i,j,k] = state.rho[i,j,k-1]
    for k in range(gc.nz+1,gc.nz+2):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.rho[i,j,k] = state.rho[i,j,k-2]

    for k in range(-2,-1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.u[i,j,k] = state.u[i,j,k+2] / staticv.hy_dens_cell[0,0,k+2] * staticv.hy_dens_cell[0,0,k]
    for k in range(-1,0):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.u[i,j,k] = state.u[i,j,k+1] / staticv.hy_dens_cell[0,0,k+1] * staticv.hy_dens_cell[0,0,k]
    for k in range(gc.nz,gc.nz+1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.u[i,j,k] = state.u[i,j,k-1] / staticv.hy_dens_cell[0,0,k-1] * staticv.hy_dens_cell[0,0,k]
    for k in range(gc.nz+1,gc.nz+2):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.u[i,j,k] = state.u[i,j,k-2] / staticv.hy_dens_cell[0,0,k-2] * staticv.hy_dens_cell[0,0,k]

    for k in range(-2,0):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.w[i,j,k] = 0.0
    for k in range(gc.nz,gc.nz+2):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.w[i,j,k] = 0.0

    for k in range(-2,-1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.pt[i,j,k] = state.pt[i,j,k+2]
    for k in range(-1,0):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.pt[i,j,k] = state.pt[i,j,k+1]
    for k in range(gc.nz,gc.nz+1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.pt[i,j,k] = state.pt[i,j,k-1]
    for k in range(gc.nz+1,gc.nz+2):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state.pt[i,j,k] = state.pt[i,j,k-2]               

@space_op
def ComputeTendenciesZ(state:phv.HybridStateField, staticv:phv.HybridStaticField, flux:phv.HybridFluxField, tend:phv.HybridTendField, dt:Float64):
    #Compute the hyperviscosity coefficient
    hv_coef = -cv.hv_beta * gc.dz / (16*dt)   

    vals = Vector4(0.0,0.0,0.0,0.0)
    d3_vals = Vector4(0.0,0.0,0.0,0.0)

    #k == 0, w = 0
    for k in range(0,1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
                vals.r = -state.rho[i,j,k-2]/12 + 7*state.rho[i,j,k-1]/12 + 7*state.rho[i,j,k]/12 - state.rho[i,j,k+1]/12
                d3_vals.r = -state.rho[i,j,k-2] + 3*state.rho[i,j,k-1] - 3*state.rho[i,j,k] + state.rho[i,j,k+1]
                vals.u = -state.u[i,j,k-2]/12 + 7*state.u[i,j,k-1]/12 + 7*state.u[i,j,k]/12 - state.u[i,j,k+1]/12
                d3_vals.u = -state.u[i,j,k-2] + 3*state.u[i,j,k-1] - 3*state.u[i,j,k] + state.u[i,j,k+1]
                vals.w = -state.w[i,j,k-2]/12 + 7*state.w[i,j,k-1]/12 + 7*state.w[i,j,k]/12 - state.w[i,j,k+1]/12
                d3_vals.w = -state.w[i,j,k-2] + 3*state.w[i,j,k-1] - 3*state.w[i,j,k] + state.w[i,j,k+1]
                vals.t = -state.pt[i,j,k-2]/12 + 7*state.pt[i,j,k-1]/12 + 7*state.pt[i,j,k]/12 - state.pt[i,j,k+1]/12
                d3_vals.t = -state.pt[i,j,k-2] + 3*state.pt[i,j,k-1] - 3*state.pt[i,j,k] + state.pt[i,j,k+1]

                #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
                r = vals.r + staticv.hy_dens_int[0,0,k]
                u = vals.u / r
                w = 0.0
                d3_vals.r = 0.0
                t = (vals.t + staticv.hy_dens_theta_int[0,0,k]) / r
                p = cv.c0 * ((r*t) ** cv.gamm) - staticv.hy_pressure_int[0,0,k]

                #Compute the flux vector
                flux.frho[i,j,k] = r * w - hv_coef * d3_vals.r
                flux.fu[i,j,k] = r * w * u - hv_coef * d3_vals.u
                flux.fw[i,j,k] = r * w * w + p - hv_coef * d3_vals.w
                flux.fpt[i,j,k] = r * w * t - hv_coef * d3_vals.t

    for k in range(1,gc.nz+1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
                vals.r = -state.rho[i,j,k-2]/12 + 7*state.rho[i,j,k-1]/12 + 7*state.rho[i,j,k]/12 - state.rho[i,j,k+1]/12
                d3_vals.r = -state.rho[i,j,k-2] + 3*state.rho[i,j,k-1] - 3*state.rho[i,j,k] + state.rho[i,j,k+1]
                vals.u = -state.u[i,j,k-2]/12 + 7*state.u[i,j,k-1]/12 + 7*state.u[i,j,k]/12 - state.u[i,j,k+1]/12
                d3_vals.u = -state.u[i,j,k-2] + 3*state.u[i,j,k-1] - 3*state.u[i,j,k] + state.u[i,j,k+1]
                vals.w = -state.w[i,j,k-2]/12 + 7*state.w[i,j,k-1]/12 + 7*state.w[i,j,k]/12 - state.w[i,j,k+1]/12
                d3_vals.w = -state.w[i,j,k-2] + 3*state.w[i,j,k-1] - 3*state.w[i,j,k] + state.w[i,j,k+1]
                vals.t = -state.pt[i,j,k-2]/12 + 7*state.pt[i,j,k-1]/12 + 7*state.pt[i,j,k]/12 - state.pt[i,j,k+1]/12
                d3_vals.t = -state.pt[i,j,k-2] + 3*state.pt[i,j,k-1] - 3*state.pt[i,j,k] + state.pt[i,j,k+1]

                #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
                r = vals.r + staticv.hy_dens_int[0,0,k]
                u = vals.u / r
                w = vals.w / r
                t = (vals.t + staticv.hy_dens_theta_int[0,0,k]) / r
                p = cv.c0 * ((r*t) ** cv.gamm) - staticv.hy_pressure_int[0,0,k]

                #Compute the flux vector
                flux.frho[i,j,k] = r * w - hv_coef * d3_vals.r
                flux.fu[i,j,k] = r * w * u - hv_coef * d3_vals.u
                flux.fw[i,j,k] = r * w * w + p - hv_coef * d3_vals.w
                flux.fpt[i,j,k] = r * w * t - hv_coef * d3_vals.t

    #k == nz, w = 0
    for k in range(gc.nz,gc.nz+1):
        for j in range(0,1):
            for i in range(0,gc.nx):
                #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
                vals.r = -state.rho[i,j,k-2]/12 + 7*state.rho[i,j,k-1]/12 + 7*state.rho[i,j,k]/12 - state.rho[i,j,k+1]/12
                d3_vals.r = -state.rho[i,j,k-2] + 3*state.rho[i,j,k-1] - 3*state.rho[i,j,k] + state.rho[i,j,k+1]
                vals.u = -state.u[i,j,k-2]/12 + 7*state.u[i,j,k-1]/12 + 7*state.u[i,j,k]/12 - state.u[i,j,k+1]/12
                d3_vals.u = -state.u[i,j,k-2] + 3*state.u[i,j,k-1] - 3*state.u[i,j,k] + state.u[i,j,k+1]
                vals.w = -state.w[i,j,k-2]/12 + 7*state.w[i,j,k-1]/12 + 7*state.w[i,j,k]/12 - state.w[i,j,k+1]/12
                d3_vals.w = -state.w[i,j,k-2] + 3*state.w[i,j,k-1] - 3*state.w[i,j,k] + state.w[i,j,k+1]
                vals.t = -state.pt[i,j,k-2]/12 + 7*state.pt[i,j,k-1]/12 + 7*state.pt[i,j,k]/12 - state.pt[i,j,k+1]/12
                d3_vals.t = -state.pt[i,j,k-2] + 3*state.pt[i,j,k-1] - 3*state.pt[i,j,k] + state.pt[i,j,k+1]

                #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
                r = vals.r + staticv.hy_dens_int[0,0,k]
                u = vals.u / r
                w = 0.0
                d3_vals.r = 0.0
                t = (vals.t + staticv.hy_dens_theta_int[0,0,k]) / r
                p = cv.c0 * ((r*t) ** cv.gamm) - staticv.hy_pressure_int[0,0,k]

                #Compute the flux vector
                flux.frho[i,j,k] = r * w - hv_coef * d3_vals.r
                flux.fu[i,j,k] = r * w * u - hv_coef * d3_vals.u
                flux.fw[i,j,k] = r * w * w + p - hv_coef * d3_vals.w
                flux.fpt[i,j,k] = r * w * t - hv_coef * d3_vals.t

    #Use the fluxes to compute tendencies for each cell
    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx):
                tend.drho[i,j,k] = -(flux.frho[i,j,k+1] - flux.frho[i,j,k]) / gc.dz
                tend.du[i,j,k] = -(flux.fu[i,j,k+1] - flux.fu[i,j,k]) / gc.dz
                tend.dw[i,j,k] = -(flux.fw[i,j,k+1] - flux.fw[i,j,k]) / gc.dz
                tend.dw[i,j,k] = tend.dw[i,j,k] - state.rho[i,j,k] * cv.grav
                tend.dpt[i,j,k] = -(flux.fpt[i,j,k+1] - flux.fpt[i,j,k]) / gc.dz

@space_op
def UpdateState(state_init:phv.HybridStateField, state_out:phv.HybridStateField, staticv:phv.HybridStaticField, tend:phv.HybridTendField, dt:Float64):
    
    #ToRemember Add DATA_SPEC_GRAVITY_WAVES
    
    for k in range(0,gc.nz):
        for j in range(0,1):
            for i in range(0,gc.nx):
                state_out.rho[i,j,k] = state_init.rho[i,j,k] + dt * tend.drho[i,j,k]
                state_out.u[i,j,k] = state_init.u[i,j,k] + dt * tend.du[i,j,k]
                state_out.w[i,j,k] = state_init.w[i,j,k] + dt * tend.dw[i,j,k]
                state_out.pt[i,j,k] = state_init.pt[i,j,k] + dt * tend.dpt[i,j,k]


                                












   
                        

                        