from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init  
from pytoy.lang.dtype_ext import LonLatField 
import math

import miniweather_scuderia.const_value as cv
import miniweather_scuderia.global_config as gc
from miniweather_scuderia.utils import MeshVector3, MeshVector6, disect



@hybrid
class HybridMeshField:
    xp : LonLatField[Float64]
    zpf: LonLatField[Float64] #full
    zph: LonLatField[Float64] #half

xp = LonLatField(Float64, (gc.nx,1,1),  (True,True,True), const = True )
zpf = LonLatField(Float64, (1,1,gc.nz),  (True,True,True), const = True )
zph = LonLatField(Float64, (1,1,gc.nz+1),  (True,True,False), const = True )

mesh = HybridMeshField(xp,zpf,zph)

@space_op 
@const_init
def MeshInit(mesh:HybridMeshField):

    #ToRemember 转换全域下标,将计算的初始网格下标对齐为0,MeshVector6前三是计算全域大小,后三是每一维的需要移的大小,[z,y,x]
    t = 0
    loopinfo = MeshVector6(1, 1, gc.nx,0,0,-2)
    for i in range(-1,1):
        v_inx = disect(t,loopinfo)
        t = t+1
        mesh.xp[i,0,0] = (v_inx.z + 0.5) * gc.dx
    for i in range(1,gc.nx):
        mesh.xp[i,0,0] = mesh.xp[i-1,0,0] + gc.dx

    t = 0
    loopinfo = MeshVector6(gc.nz, 1, 1,-2,0,0)
    for k in range(-2,1):
        v_inx = disect(t,loopinfo)
        t = t+1
        mesh.zpf[0,0,k] = (v_inx.x + 0.5) * gc.dz
    for k in range(1,gc.nz+2):
        mesh.zpf[0,0,k] = mesh.zpf[0,0,k-1] + gc.dz

    t = 0
    loopinfo = MeshVector6(gc.nz, 1, 1,-2,0,0)
    for k in range(-2,1):
        v_inx = disect(t,loopinfo)
        t = t+1
        mesh.zph[0,0,k] = (v_inx.x) * gc.dz
    for k in range(1,gc.nz+1):
        mesh.zph[0,0,k] = mesh.zph[0,0,k-1] + gc.dz
    
