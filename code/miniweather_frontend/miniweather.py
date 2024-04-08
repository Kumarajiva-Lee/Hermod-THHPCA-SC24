import copy
import math
from typing import Dict, List, cast
from pytoy.lang import time_op, space_op, Float64, Int32, Field, hybrid, const_init, func, extern_func, const_init
from pytoy.lang.dtype_ext import LonLatField
from pytoy.lang.annot import Float64
from pytoy.optim.proc.globalanalyse import Proc
from pytoy.configuration import GlobalConfiguration as gc
from pytoy.codegen.process import ToBackend


import miniweather_scuderia.time_operator as to
import miniweather_scuderia.global_config as gc
from miniweather_scuderia.mesh import mesh
from miniweather_scuderia.physical_variable import state, staticv, stateout, tend, flux, DirPara, OutPara


to.MeshInit(mesh)
to.timeOperatorInit(state, staticv, mesh, DirPara)
to.performTimestep(state, staticv, flux, tend, mesh, DirPara, stateout, OutPara, gc.dt, gc.dt/2, gc.dt/3)


ToBackend()