from pytoy.lang import time_op, const_init, Float64, Int32

import frontend.mcv_stroud.mesh as gm
import frontend.mcv_stroud.physical_variable as phv
import frontend.mcv_stroud.space_operator as sp
import frontend.mcv_stroud.operator_comp as oc

@time_op(1,1)
@const_init
def MeshInit(mesh:gm.HybridMeshField):
    gm.MeshInit(mesh)

@time_op(1,1)
def OperatorInit(state: phv.HybridStateField, mesh: gm.HybridMeshField):
    sp.spaceOperatorInit(state[0],mesh)

@time_op(30,1)
def Rk3(state: phv.HybridStateField, tend: phv.HybridTendField, mesh: gm.HybridMeshField, dt:Float64):
    oc.updateState(state[0],state[2])
    sp.mcvupdateXY(state[0],tend[0],mesh)
    oc.update_rk1(state[0], state[2], tend[0], dt)
    sp.mcvupdateXY(state[2],tend[1],mesh)
    oc.update_rk2(state[0], state[2], tend[0], tend[1], dt)
    sp.mcvupdateXY(state[2],tend[2],mesh)
    oc.update_rk3(state[0], state[1], tend[0], tend[1], tend[2], dt)
    # sp.mcvupdateXY(state[0],tend[0],mesh)
    # oc.update_rk1(state[0], state[2], tend[0], dt)
    # sp.mcvupdateXY(state[0],tend[1],mesh)
    # oc.update_rk2(state[0], state[2], tend[0], tend[1], dt)
    # sp.mcvupdateXY(state[0],tend[2],mesh)
    # oc.update_rk3(state[0], state[1], tend[0], tend[1], tend[2], dt)


