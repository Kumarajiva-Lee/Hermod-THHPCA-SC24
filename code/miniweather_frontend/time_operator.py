from pytoy.lang import time_op, const_init, Float64, Int32

import miniweather_scuderia.physical_variable as phv
import miniweather_scuderia.space_operatpr as so
import miniweather_scuderia.operator_comp as oc
import miniweather_scuderia.mesh as gm
import miniweather_scuderia.ncIO as nc

@time_op(1,1)
@const_init
def MeshInit(mesh:gm.HybridMeshField):
    gm.MeshInit(mesh)

@time_op(1,1)
def timeOperatorInit(state:phv.HybridStateField,staticv:phv.HybridStaticField, mesh:gm.HybridMeshField, DirPara:phv.HybridDirPara):
    so.spaceOperatorInit(state[0],state[1],staticv[-1],mesh,DirPara)

@time_op(20,1)
def performTimestep(state:phv.HybridStateField, staticv:phv.HybridStaticField, flux:phv.HybridFluxField, tend:phv.HybridTendField, mesh:gm.HybridMeshField, DirPara:phv.HybridDirPara, stateoutput:phv.HybridOutField, OutPara:phv.HybridOutPara, dt:Float64 , dtd2: Float64, dtd3: Float64):
    so.SemiDiscreteStep(state[0],state[0],state[1],staticv[-1],flux[-1],tend[-1],DirPara,dtd3)
    so.SemiDiscreteStep(state[0],state[1],state[1],staticv[-1],flux[-1],tend[-1],DirPara,dtd2)
    so.SemiDiscreteStep(state[0],state[1],state[0],staticv[-1],flux[-1],tend[-1],DirPara,dt)

    so.SemiDiscreteStep(state[0],state[0],state[1],staticv[-1],flux[-1],tend[-1],DirPara,dtd3)
    so.SemiDiscreteStep(state[0],state[1],state[1],staticv[-1],flux[-1],tend[-1],DirPara,dtd2)
    so.SemiDiscreteStep(state[0],state[1],state[0],staticv[-1],flux[-1],tend[-1],DirPara,dt)

    oc.SwitchDirection(DirPara)

    nc.OutputPrepare(state[0],stateoutput[0],staticv,OutPara)