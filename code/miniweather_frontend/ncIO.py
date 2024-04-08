from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init
import miniweather_scuderia.const_value as cv
import miniweather_scuderia.global_config as gc
import miniweather_scuderia.physical_variable as phv

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
def OutputPrepare(state:phv.HybridStateField, stateout:phv.HybridOutField, staticv:phv.HybridStaticField, OutPara:phv.HybridOutPara):
    
    OutPara.step = OutPara.step + 1

    if OutPara.step == OutPara.interval:

        for k in range(0,gc.nz):
            for j in range(0,1):
                for i in range(0,gc.nx):
                    stateout.dens[i,j,k] = state.rho[i,j,k]
                    stateout.uwnd[i,j,k] = state.u[i,j,k] / (staticv.hy_dens_cell[0,0,k] + state.rho[i,j,k])
                    stateout.wwnd[i,j,k] = state.w[i,j,k] / (staticv.hy_dens_cell[0,0,k] + state.rho[i,j,k])
                    stateout.theta[i,j,k] = (state.pt[i,j,k] + staticv.hy_dens_theta_cell[0,0,k]) / (staticv.hy_dens_cell[0,0,k] + state.rho[i,j,k]) - staticv.hy_dens_theta_cell[0,0,k] / staticv.hy_dens_cell[0,0,k]

        ncIO_Write(stateout)

        OutPara.step = 0



@space_op
def ncIO_Write(state:phv.HybridOutField):
    nid = ncOpenFile("output.nc", 1)
    ncDefDim(nid, gc.nz, 1, gc.nx)
    state.dens.ncWrite(nid)
    state.uwnd.ncWrite(nid)
    state.wwnd.ncWrite(nid)
    state.theta.ncWrite(nid)
    ncCloseFile(nid)
