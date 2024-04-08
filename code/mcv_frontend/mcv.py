from pytoy.lang import time_op, space_op, Float64, Int32, Field, hybrid, const_init, func, extern_func, const_init
from pytoy.codegen.process import ToBackend

import frontend.mcv_stroud.global_config as gc
import frontend.mcv_stroud.const_value as cv
from frontend.mcv_stroud.mesh import mesh
from frontend.mcv_stroud.physical_variable import state, tend
import frontend.mcv_stroud.time_operator as to


to.MeshInit(mesh)
to.OperatorInit(state, mesh)
to.Rk3(state, tend, mesh, gc.dt)

ToBackend()