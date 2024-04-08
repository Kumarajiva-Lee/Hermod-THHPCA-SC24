import  gmcore_smallball.global_config as gc
import  gmcore_smallball.time_operator as to
from gmcore_smallball.physical_variable import state,staticv,tend,tendPara, adv, advptPara, advmPara
from  gmcore_smallball.mesh import mesh

from pytoy.codegen.process import ToBackend



to.MeshInit(mesh)
to.timeOperatorInit(state,staticv,adv,advptPara,mesh,gc.dt_dyn)

to.wrk3d(state,staticv,tend,tendPara,adv,advptPara,advmPara,mesh,gc.dt_dyn,gc.dt_dyn/2.0,gc.dt_dyn/3.0)


ToBackend()