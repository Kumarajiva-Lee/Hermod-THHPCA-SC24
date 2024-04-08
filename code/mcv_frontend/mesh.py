import math

from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init  
from pytoy.lang.dtype_ext import LonLatField, CubedSphereField 

import frontend.mcv_stroud.global_config as gc
import frontend.mcv_stroud.const_value as cv
import frontend.mcv_stroud.math_formual as mf
from frontend.mcv_stroud.utils import MeshVector3, MeshVector6, disect

@hybrid
class HybridMeshField:
    #网格距
    x_1:CubedSphereField[Float64]
    x_2:CubedSphereField[Float64]
    x_3:CubedSphereField[Float64]

    y_1:CubedSphereField[Float64]
    y_2:CubedSphereField[Float64]
    y_3:CubedSphereField[Float64]


    # xh : CubedSphereField[Float64]
    # xf : CubedSphereField[Float64]
    # yh : CubedSphereField[Float64]
    # yf : CubedSphereField[Float64]
    #Jacobi变换
    # lt_ff : CubedSphereField[Float64]
    # rt_ff : CubedSphereField[Float64]
    # rb_ff : CubedSphereField[Float64]



x_1 = CubedSphereField(Float64, (gc.nx,1,1),  (True,True,True), const = True )
x_2 = CubedSphereField(Float64, (gc.nx,1,1),  (True,True,True), const = True )
x_3 = CubedSphereField(Float64, (gc.nx,1,1),  (True,True,True), const = True )

y_1 = CubedSphereField(Float64, (1,gc.ny,1),  (True,True,True), const = True )
y_2 = CubedSphereField(Float64, (1,gc.ny,1),  (True,True,True), const = True )
y_3 = CubedSphereField(Float64, (1,gc.ny,1),  (True,True,True), const = True )



# lt_ff = CubedSphereField(Float64, (gc.nx+1,gc.ny+1,1),  (False,False,True), const = True, usepanel = True)
# rt_ff = CubedSphereField(Float64, (gc.nx+1,gc.ny+1,1),  (False,False,True), const = True, usepanel = True )
# rb_ff = CubedSphereField(Float64, (gc.nx+1,gc.ny+1,1),  (False,False,True), const = True, usepanel = True )

mesh = HybridMeshField(x_1,x_2,x_3,y_1,y_2,y_3)

@space_op 
@const_init
def MeshInit(mesh:HybridMeshField):

    for i in range(-1,gc.nx):
        mesh.x_1[i,0,0] = -cv.r * cv.pi / 4.0 + (i-1) * cv.dx
        mesh.x_2[i,0,0] = mesh.x_1[i,0,0] + cv.dxi
        mesh.x_3[i,0,0] = mesh.x_2[i,0,0] + cv.dxi

    for j in range(-1,gc.ny):
        mesh.y_1[0,j,0] = -cv.r * cv.pi / 4.0 + (j-1) * cv.dy
        mesh.y_2[0,j,0] = mesh.y_1[0,j,0] + cv.dyi
        mesh.y_3[0,j,0] = mesh.y_2[0,j,0] + cv.dyi

    # for i in range(-1,gc.nx-1):
    #     mesh.xh[i,0,0] = -cv.r * cv.pi / 4.0 + (i-0.5) * cv.dx       
    # for j in range(-1,gc.ny):
    #     mesh.yf[0,j,0] = -cv.r * cv.pi / 4.0 + (j-1) * cv.dy
    # for j in range(-1,gc.ny-1):
    #     mesh.yh[0,j,0] = -cv.r * cv.pi / 4.0 + (j-0.5) * cv.dy
    

    # if pypanel == 1 or pypanel == 2: 
    #     for i in range(-1,gc.nx):
    #         mesh.xf[i,0,0] = -cv.r * cv.pi / 4.0 + (i-1) * cv.dx
    #     for i in range(-1,gc.nx-1):
    #         mesh.xh[i,0,0] = -cv.r * cv.pi / 4.0 + (i-0.5) * cv.dx       
    #     for j in range(-1,gc.ny):
    #         mesh.yf[0,j,0] = -cv.r * cv.pi / 4.0 + (j-1) * cv.dy
    #     for j in range(-1,gc.ny-1):
    #         mesh.yh[0,j,0] = -cv.r * cv.pi / 4.0 + (j-0.5) * cv.dy
    # elif pypanel == 4:
    #     #(x,y) -> (y, n-x)
    #     for i in range(-1,gc.nx):
    #         mesh.xf[i,0,0] = -cv.r * cv.pi / 4.0 + (i-1) * cv.dy
    #     for i in range(-1,gc.nx-1):
    #         mesh.xh[i,0,0] = -cv.r * cv.pi / 4.0 + (i-0.5) * cv.dy      
    #     for j in range(-1,gc.ny):
    #         mesh.yf[0,j,0] = -cv.r * cv.pi / 4.0 + (gc.ny + 2 - j - 1) * cv.dx
    #     for j in range(-1,gc.ny-1):
    #         mesh.yh[0,j,0] = -cv.r * cv.pi / 4.0 + (gc.ny-1 + 2 - j-0.5) * cv.dx




