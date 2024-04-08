from frontend.pytoy.lang import space_op, Float64, Int32

import miniweather_scuderia.physical_variable as phv
import miniweather_scuderia.operator_comp as oc
import miniweather_scuderia.mesh as gm
import miniweather_scuderia.global_config as gc

@space_op
def spaceOperatorInit(state:phv.HybridStateField, state_tmp:phv.HybridStateField, staticv:phv.HybridStaticField, mesh:gm.HybridMeshField, DirPara:phv.HybridDirPara):
    oc.StateInitial(state, state_tmp, staticv, mesh, DirPara)

@space_op
def SemiDiscreteStep(state_init:phv.HybridStateField, state_forcing:phv.HybridStateField, state_out:phv.HybridStateField, staticv:phv.HybridStaticField, flux:phv.HybridFluxField, tend:phv.HybridTendField,  DirPara:phv.HybridDirPara, dt: Float64):
    if DirPara.xd:
        if DirPara.stepnum < 3:
            oc.ComputeTendenciesX(state_forcing,staticv,flux,tend,dt)
        else:
            oc.SetBoundaryZ(state_forcing,staticv)
            oc.ComputeTendenciesZ(state_forcing,staticv,flux,tend,dt)
    else:
        if DirPara.stepnum < 3:
            oc.SetBoundaryZ(state_forcing,staticv)
            oc.ComputeTendenciesZ(state_forcing,staticv,flux,tend,dt)
        else:
            oc.ComputeTendenciesX(state_forcing,staticv,flux,tend,dt)

    oc.UpdateState(state_init,state_out,staticv,tend,dt)

    DirPara.stepnum = DirPara.stepnum + 1

