from frontend.pytoy.lang import space_op, Float64, Int32

import gmcore_smallball.global_config as gc
import gmcore_smallball.mesh as gm
import gmcore_smallball.physical_variable as phv
import gmcore_smallball.operator_comp as oc

@space_op
def spaceOperatorInit(staticv:phv.HybridStaticField,mesh:gm.HybridMeshField):
    oc.prepareStatic(staticv,mesh)

#需要同时预处理所有时间维度,故单独处理
@space_op
def gzlevInit(state:phv.HybridStateField,staticv:phv.HybridStaticField):
    oc.preparegzlev(state,staticv)

@space_op
def uvInit(state:phv.HybridStateField):
    oc.c2a(state)

@space_op
def advPrepare(state:phv.HybridStateField, advm: phv.HybridAdvField, advpt: phv.HybridAdvField):
    oc.copy_old_m(state,advm)
    oc.copy_old_m(state,advpt)

@space_op
def stepForwardBackward(old_state:phv.HybridStateField, star_state:phv.HybridStateField, new_state:phv.HybridStateField, staticv:phv.HybridStaticField, tend1:phv.HybridTendField, tend2:phv.HybridTendField, tendPara: phv.HybridTendPara, advpt: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    spaceoperatorForward(old_state,star_state,new_state,tend1,tend2,tendPara,advpt,advPara,mesh,dt)
    oc.update_state(old_state,new_state,tend1,tendPara,mesh,dt)
    spaceoperatorBackward(old_state,star_state,new_state,staticv,tend2,tend1,tendPara,advpt,advPara,mesh,dt)
    oc.update_state(old_state,new_state,tend2,tendPara,mesh,dt)

@space_op
def spaceoperatorForward(old_state:phv.HybridStateField, star_state:phv.HybridStateField, new_state:phv.HybridStateField, tend1:phv.HybridTendField, tend2:phv.HybridTendField, tendPara:phv.HybridTendPara, advpt: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    oc.reset_flag(tendPara)
    operatorPrepareForward(star_state,advpt,advPara,mesh,dt)
    oc.calc_grad_mf(star_state,tend1,mesh)
    oc.calc_dphsdt(tend1)
    oc.calc_we_lev(star_state,tend1,advpt,advPara,mesh,dt)
    oc.calc_wedudlev_wedvdlev(star_state,tend1)
    oc.calc_grad_ptf(star_state,tend1,advpt,mesh,dt)
    oc.calc_coriolis(star_state,tend1,mesh)
    oc.calc_grad_ke(star_state,tend1,mesh)
    oc.calc_tend_forward(tend1,tendPara)

@space_op
def spaceoperatorBackward(old_state:phv.HybridStateField, star_state:phv.HybridStateField, new_state:phv.HybridStateField, staticv:phv.HybridStaticField, tend1:phv.HybridTendField, tend2:phv.HybridTendField, tendPara:phv.HybridTendPara, advpt: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    oc.reset_flag(tendPara)
    operatorPrepareBackward(new_state,staticv)
    oc.pgf_lin97(new_state,tend1,mesh)
    oc.calc_tend_backward(tend1,tend2,tendPara)

@space_op
def operatorPrepareNull(state:phv.HybridStateField,staticv:phv.HybridStaticField, advpt: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    oc.calc_ph(state,mesh)
    oc.calc_m(state,mesh)
    oc.calc_t(state)
    oc.calc_mf(state,advpt,advPara,mesh,dt)
    oc.calc_ke(state,mesh)
    oc.calc_pv(state,mesh)
    oc.interp_pv(state)
    oc.calc_div(state,mesh)
    oc.calc_gz_lev(state,staticv)

@space_op
def operatorPrepareForward(state:phv.HybridStateField, adv: phv.HybridAdvField, advPara: phv.HybridAdvPara, mesh:gm.HybridMeshField, dt:Float64):
    oc.calc_mf(state,adv,advPara,mesh,dt)
    oc.calc_ke(state,mesh)
    oc.calc_div(state,mesh)
    oc.calc_pv(state,mesh)
    oc.interp_pv(state)

@space_op
def operatorPrepareBackward(state:phv.HybridStateField, staticv:phv.HybridStaticField):
    oc.calc_t(state)
    oc.calc_gz_lev(state,staticv)

@space_op
def damp_run(state:phv.HybridStateField, tend:phv.HybridTendField, mesh:gm.HybridMeshField, dt:Float64):
    oc.div_damp_run(state,mesh)
    oc.smag_damp_run(state,tend,mesh,dt)
    oc.pole_damp_run(state,mesh)