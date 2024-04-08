from frontend.pytoy.lang import space_op, Float64, Int32

import frontend.mcv_stroud.const_value as cv
import frontend.mcv_stroud.global_config as gc
import frontend.mcv_stroud.physical_variable as phv
import frontend.mcv_stroud.mesh as gm
import frontend.mcv_stroud.operator_comp as oc



@space_op
def spaceOperatorInit(state:phv.HybridStateField,  mesh:gm.HybridMeshField):
    oc.StateInitial(state, mesh)

@space_op
def mcvupdateXY(state: phv.HybridStateField, tend: phv.HybridTendField, mesh: gm.HybridMeshField):
    oc.updateX(state, tend)
    oc.updateY(state, tend)
    # oc.updateUV(state, tend, mesh)
