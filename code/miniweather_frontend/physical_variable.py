import miniweather_scuderia.global_config as gc
from frontend.pytoy.lang import Float64, Bool, Int32, hybrid   
from frontend.pytoy.lang.dtype_ext import LonLatField 

@hybrid
class HybridStateField:
    u    : LonLatField[Float64] #momentum in the x-direction ("rho * u")
    w    : LonLatField[Float64] #momentum in the z-direction ("rho * w")
    rho  : LonLatField[Float64] #density ("rho")
    pt   : LonLatField[Float64] #density * potential temperature ("rho * theta")

u = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
w = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
rho = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
pt = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))

state = HybridStateField(u,w,rho,pt)

@hybrid
class HybridStaticField:
    hy_dens_cell        : LonLatField[Float64]
    hy_dens_theta_cell  : LonLatField[Float64]
    hy_dens_int         : LonLatField[Float64]
    hy_dens_theta_int   : LonLatField[Float64]
    hy_pressure_int     : LonLatField[Float64]

#ToCheck Lev Full or Half
hy_dens_cell = LonLatField(Float64,(1,1,gc.nz),(True,True,True),const = True)
hy_dens_theta_cell = LonLatField(Float64,(1,1,gc.nz),(True,True,True),const = True)
hy_dens_int = LonLatField(Float64,(1,1,gc.nz+1),(True,True,False),const = True)
hy_dens_theta_int = LonLatField(Float64,(1,1,gc.nz+1),(True,True,False),const = True)
hy_pressure_int = LonLatField(Float64,(1,1,gc.nz+1),(True,True,False),const = True)

staticv = HybridStaticField(hy_dens_cell, hy_dens_theta_cell, hy_dens_int, hy_dens_theta_int, hy_pressure_int)

@hybrid
class HybridFluxField:
    fu    : LonLatField[Float64]
    fw    : LonLatField[Float64]
    frho  : LonLatField[Float64]
    fpt   : LonLatField[Float64]

fu = LonLatField(Float64,(gc.nx,1,gc.nz),(False,True,False))
fw = LonLatField(Float64,(gc.nx,1,gc.nz),(False,True,False))
frho = LonLatField(Float64,(gc.nx,1,gc.nz),(False,True,False))
fpt = LonLatField(Float64,(gc.nx,1,gc.nz),(False,True,False))

flux = HybridFluxField(fu,fw,frho,fpt)



@hybrid
class HybridTendField:
    du    : LonLatField[Float64]
    dw    : LonLatField[Float64]
    drho  : LonLatField[Float64]
    dpt   : LonLatField[Float64]

du = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
dw = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
drho = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
dpt = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))

tend = HybridTendField(du,dw,drho,dpt)

@hybrid
class HybridOutField:
    dens    : LonLatField[Float64]
    uwnd    : LonLatField[Float64]
    wwnd  : LonLatField[Float64]
    theta   : LonLatField[Float64] 

dens = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
uwnd = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
wwnd = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))
theta = LonLatField(Float64,(gc.nx,1,gc.nz),(True,True,True))

stateout = HybridOutField(dens,uwnd,wwnd,theta)

@hybrid
class HybridOutPara:
    step : Int32
    interval : Int32

OutPara = HybridOutPara(0,2000)

@hybrid
class HybridDirPara:
    xd : Bool
    zd : Bool
    stepnum : Int32

DirPara = HybridDirPara(False,False,0)
