from pytoy.lang import time_op, const_init, Float64, Int32

import gmcore_smallball.physical_variable as phv
import gmcore_smallball.space_operator as so
import gmcore_smallball.mesh as gm
import gmcore_smallball.global_config as gc
import gmcore_smallball.operator_comp as oc
# import vert_coord as vc
import gmcore_smallball.ncIO as nc
import gmcore_smallball.filter as filter

@time_op(1,1)
@const_init
#def timeOperatorInit(state:phv.HybridStateField,tend:phv.HybridTendField,static:phv.HybridStaticField,mesh:gm.HybridMeshField,filter_field:filter.FilterField):
def MeshInit(mesh:gm.HybridMeshField):
    gm.LatLonMeshInit(mesh)


@time_op(1,1)
def timeOperatorInit(state:phv.HybridStateField, staticv:phv.HybridStaticField, adv: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    nc.ncIO_Read(state,staticv)
    #nc.ncIO_Write(state,staticv)
    so.spaceOperatorInit(staticv[-1],mesh)
    so.gzlevInit(state[0],staticv[-1])
    so.gzlevInit(state[1],staticv[-1])
    so.gzlevInit(state[2],staticv[-1])
    so.uvInit(state[0])
    so.operatorPrepareNull(state[0],staticv[-1],adv[-2],advPara,mesh,dt)
    so.advPrepare(state[0],adv[-1],adv[-2])
    filter.filterinit(mesh,dt)

#在timeop中无法运算,故传入多个dt参数,dtd2 = dt / 2
@time_op(30,1)
def wrk3d(state:phv.HybridStateField, staticv:phv.HybridStaticField, tend:phv.HybridAdvField, tendPara:phv.HybridTendPara, adv: phv.HybridAdvField, advptPara: phv.HybridAdvPara, advmPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64 , dtd2: Float64, dtd3: Float64):
    so.stepForwardBackward(state[0],state[0],state[1],staticv[-1],tend[0],tend[1],tendPara,adv[-2],advptPara,mesh,dtd3)
    so.stepForwardBackward(state[0],state[1],state[2],staticv[-1],tend[0],tend[1],tendPara,adv[-2],advptPara,mesh,dtd2)
    so.stepForwardBackward(state[0],state[2],state[1],staticv[-1],tend[0],tend[1],tendPara,adv[-2],advptPara,mesh,dt)
    so.damp_run(state[1],tend[1],mesh,dt)
    oc.c2a(state[1])
    oc.timeadvance()
    oc.adv_run(state[0],state[1],adv[-1],advmPara,mesh,dt)
    oc.Diagnose_State(state[0])

