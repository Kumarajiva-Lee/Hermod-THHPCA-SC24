from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init
from frontend.pytoy.lang.dtype_ext import LonLatField
import gmcore_smallball.mesh as gm


@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 1, parallel = True)
def Filter_Init(mesh:gm.HybridMeshField, dt:Float64):
    pass

@space_op
def filterinit(mesh:gm.HybridMeshField, dt:Float64):
    Filter_Init(mesh,dt)

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 1, parallel = True, name="filter")
def Filter_On_Cell(id:Int32, is3d:Int32, xfield: LonLatField[Float64]):
    pass

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 1, parallel = True, name="filter")
def Filter_on_lon_edge(id:Int32, xfield: LonLatField[Float64]):
    pass

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 1, parallel = True, name="filter")
def Filter_on_lat_edge(id:Int32, xfield: LonLatField[Float64]):
    pass