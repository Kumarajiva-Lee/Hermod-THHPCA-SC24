from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init
import gmcore_smallball.physical_variable as phv
import gmcore_smallball.global_config as gc

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 0, parallel = True)
def ncOpenFile(file_path:String, is_write:Int32)->Int32:
    pass

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 0, parallel = True)
def ncCloseFile(file_id:Int32):
    pass

@extern_func(avg_malloc = 0, avg_mem = 0, avg_flops = 0, parallel = True)
def ncDefDim(file_id:Int32, num_lev:Int32, num_lat:Int32, num_lon:Int32):
    pass

@space_op
def ncIO_Read(state:phv.HybridStateField,staticv:phv.HybridStaticField):
    nid = ncOpenFile(gc.ncIn, 0)
    state.u_lon.ncRead(nid)
    state.v_lat.ncRead(nid)
    state.pt.ncRead(nid)
    state.phs.ncRead(nid)
    staticv.gzs.ncRead(nid)
    ncCloseFile(nid)

@space_op
def ncIO_Write(state:phv.HybridStateField,staticv:phv.HybridStaticField):
    nid = ncOpenFile("gmcore_result.nc", 1)
    ncDefDim(nid, gc.nlev, gc.nlat, gc.nlon)
    state.u_lon.ncWrite(nid)
    state.v_lat.ncWrite(nid)
    state.pt.ncWrite(nid)
    state.phs.ncWrite(nid)
    staticv.gzs.ncWrite(nid)
    ncCloseFile(nid)


