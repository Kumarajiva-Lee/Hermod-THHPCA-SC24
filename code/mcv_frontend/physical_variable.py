from frontend.pytoy.lang import Float64, Bool, Int32, hybrid   
from frontend.pytoy.lang.dtype_ext import LonLatField, CubedSphereField

import frontend.mcv_stroud.global_config as gc

@hybrid
class HybridStateField:

    #前i后j
    h_11: CubedSphereField[Float64]
    u_11: CubedSphereField[Float64]
    v_11: CubedSphereField[Float64]
    jab_11:CubedSphereField[Float64]
    jab5_11:CubedSphereField[Float64]
    jab6_11:CubedSphereField[Float64]
    jab7_11:CubedSphereField[Float64]
    uc_11:CubedSphereField[Float64]
    vc_11:CubedSphereField[Float64]

    h_12: CubedSphereField[Float64]
    u_12: CubedSphereField[Float64]
    v_12: CubedSphereField[Float64]
    jab_12:CubedSphereField[Float64]
    jab5_12:CubedSphereField[Float64]
    jab6_12:CubedSphereField[Float64]
    jab7_12:CubedSphereField[Float64]
    uc_12:CubedSphereField[Float64]
    vc_12:CubedSphereField[Float64]

    h_13: CubedSphereField[Float64]
    u_13: CubedSphereField[Float64]
    v_13: CubedSphereField[Float64]
    jab_13:CubedSphereField[Float64]
    jab5_13:CubedSphereField[Float64]
    jab6_13:CubedSphereField[Float64]
    jab7_13:CubedSphereField[Float64]
    uc_13:CubedSphereField[Float64]
    vc_13:CubedSphereField[Float64]

    h_21: CubedSphereField[Float64]
    u_21: CubedSphereField[Float64]
    v_21: CubedSphereField[Float64]
    jab_21:CubedSphereField[Float64]
    jab5_21:CubedSphereField[Float64]
    jab6_21:CubedSphereField[Float64]
    jab7_21:CubedSphereField[Float64]
    uc_21:CubedSphereField[Float64]
    vc_21:CubedSphereField[Float64]

    h_22: CubedSphereField[Float64]
    u_22: CubedSphereField[Float64]
    v_22: CubedSphereField[Float64]
    jab_22:CubedSphereField[Float64]
    jab5_22:CubedSphereField[Float64]
    jab6_22:CubedSphereField[Float64]
    jab7_22:CubedSphereField[Float64]
    uc_22:CubedSphereField[Float64]
    vc_22:CubedSphereField[Float64]

    h_23: CubedSphereField[Float64]
    u_23: CubedSphereField[Float64]
    v_23: CubedSphereField[Float64]
    jab_23:CubedSphereField[Float64]
    jab5_23:CubedSphereField[Float64]
    jab6_23:CubedSphereField[Float64]
    jab7_23:CubedSphereField[Float64]
    uc_23:CubedSphereField[Float64]
    vc_23:CubedSphereField[Float64]

    h_31: CubedSphereField[Float64]
    u_31: CubedSphereField[Float64]
    v_31: CubedSphereField[Float64]
    jab_31:CubedSphereField[Float64]
    jab5_31:CubedSphereField[Float64]
    jab6_31:CubedSphereField[Float64]
    jab7_31:CubedSphereField[Float64]
    uc_31:CubedSphereField[Float64]
    vc_31:CubedSphereField[Float64]

    h_32: CubedSphereField[Float64]
    u_32: CubedSphereField[Float64]
    v_32: CubedSphereField[Float64]
    jab_32:CubedSphereField[Float64]
    jab5_32:CubedSphereField[Float64]
    jab6_32:CubedSphereField[Float64]
    jab7_32:CubedSphereField[Float64]
    uc_32:CubedSphereField[Float64]
    vc_32:CubedSphereField[Float64]

    h_33: CubedSphereField[Float64]
    u_33: CubedSphereField[Float64]
    v_33: CubedSphereField[Float64]
    jab_33:CubedSphereField[Float64]
    jab5_33:CubedSphereField[Float64]
    jab6_33:CubedSphereField[Float64]
    jab7_33:CubedSphereField[Float64]
    uc_33:CubedSphereField[Float64]
    vc_33:CubedSphereField[Float64]


h_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

h_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
u_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
v_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab5_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab6_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
jab7_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
uc_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
vc_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

state = HybridStateField(h_11, u_11, v_11, jab_11, jab5_11, jab6_11, jab7_11, uc_11, vc_11,
                         h_12, u_12, v_12, jab_12, jab5_12, jab6_12, jab7_12, uc_12, vc_12,
                         h_13, u_13, v_13, jab_13, jab5_13, jab6_13, jab7_13, uc_13, vc_13,
                         h_21, u_21, v_21, jab_21, jab5_21, jab6_21, jab7_21, uc_21, vc_21,
                         h_22, u_22, v_22, jab_22, jab5_22, jab6_22, jab7_22, uc_22, vc_22,
                         h_23, u_23, v_23, jab_23, jab5_23, jab6_23, jab7_23, uc_23, vc_23,
                         h_31, u_31, v_31, jab_31, jab5_31, jab6_31, jab7_31, uc_31, vc_31,
                         h_32, u_32, v_32, jab_32, jab5_32, jab6_32, jab7_32, uc_32, vc_32,
                         h_33, u_33, v_33, jab_33, jab5_33, jab6_33, jab7_33, uc_33, vc_33)


# state = HybridStateField(h_ff, u_ff, v_ff, jab_ff, jab5_ff, jab6_ff, jab7_ff, uc_ff, vc_ff,
#                          h_hf, u_hf, v_hf, jab_hf, jab5_hf, jab6_hf, jab7_hf, uc_hf, vc_hf,
#                          h_fh, u_fh, v_fh, jab_fh, jab5_fh, jab6_fh, jab7_fh, uc_fh, vc_fh,
#                          h_hh, u_hh, v_hh, jab_hh, jab5_hh, jab6_hh, jab7_hh, uc_hh, vc_hh)

@hybrid
class HybridTendField:
    #前i后j
    #x
    dxh_11: CubedSphereField[Float64]
    dxu_11: CubedSphereField[Float64]
    dxv_11: CubedSphereField[Float64]
    dxh_21: CubedSphereField[Float64]
    dxu_21: CubedSphereField[Float64]
    dxv_21: CubedSphereField[Float64]
    dxh_31: CubedSphereField[Float64]
    dxu_31: CubedSphereField[Float64]
    dxv_31: CubedSphereField[Float64]

    dxh_12: CubedSphereField[Float64]
    dxu_12: CubedSphereField[Float64]
    dxv_12: CubedSphereField[Float64]
    dxh_22: CubedSphereField[Float64]
    dxu_22: CubedSphereField[Float64]
    dxv_22: CubedSphereField[Float64]
    dxh_32: CubedSphereField[Float64]
    dxu_32: CubedSphereField[Float64]
    dxv_32: CubedSphereField[Float64]

    dxh_13: CubedSphereField[Float64]
    dxu_13: CubedSphereField[Float64]
    dxv_13: CubedSphereField[Float64]
    dxh_23: CubedSphereField[Float64]
    dxu_23: CubedSphereField[Float64]
    dxv_23: CubedSphereField[Float64]
    dxh_33: CubedSphereField[Float64]
    dxu_33: CubedSphereField[Float64]
    dxv_33: CubedSphereField[Float64]

    #y
    dyh_11: CubedSphereField[Float64]
    dyu_11: CubedSphereField[Float64]
    dyv_11: CubedSphereField[Float64]
    dyh_21: CubedSphereField[Float64]
    dyu_21: CubedSphereField[Float64]
    dyv_21: CubedSphereField[Float64]
    dyh_31: CubedSphereField[Float64]
    dyu_31: CubedSphereField[Float64]
    dyv_31: CubedSphereField[Float64]

    dyh_12: CubedSphereField[Float64]
    dyu_12: CubedSphereField[Float64]
    dyv_12: CubedSphereField[Float64]
    dyh_22: CubedSphereField[Float64]
    dyu_22: CubedSphereField[Float64]
    dyv_22: CubedSphereField[Float64]
    dyh_32: CubedSphereField[Float64]
    dyu_32: CubedSphereField[Float64]
    dyv_32: CubedSphereField[Float64]

    dyh_13: CubedSphereField[Float64]
    dyu_13: CubedSphereField[Float64]
    dyv_13: CubedSphereField[Float64]
    dyh_23: CubedSphereField[Float64]
    dyu_23: CubedSphereField[Float64]
    dyv_23: CubedSphereField[Float64]
    dyh_33: CubedSphereField[Float64]
    dyu_33: CubedSphereField[Float64]
    dyv_33: CubedSphereField[Float64]

    #f和q只有一个方向,第二个数对应底mcv代码中f/q的最后一维
    fh_1:CubedSphereField[Float64]
    fh_2:CubedSphereField[Float64]
    qh_1:CubedSphereField[Float64]
    qh_2:CubedSphereField[Float64]
    fu_1:CubedSphereField[Float64]
    fu_2:CubedSphereField[Float64]
    qu_1:CubedSphereField[Float64]
    qu_2:CubedSphereField[Float64]
    fv_1:CubedSphereField[Float64]
    fv_2:CubedSphereField[Float64]
    qv_1:CubedSphereField[Float64]
    qv_2:CubedSphereField[Float64]

#x
dxh_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

dxh_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

dxh_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxh_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxu_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dxv_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

#y
dyh_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_11 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_21 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_31 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

dyh_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_12 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_22 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_32 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

dyh_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_13 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_23 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyh_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyu_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
dyv_33 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

fh_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
fh_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qh_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qh_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

fu_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
fu_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qu_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qu_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

fv_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
fv_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qv_1 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))
qv_2 = CubedSphereField(Float64,(gc.nx,gc.ny,1),(True,True,True))

tend = HybridTendField(dxh_11,dxu_11,dxv_11,dxh_21,dxu_21,dxv_21,dxh_31,dxu_31,dxv_31,
                       dxh_12,dxu_12,dxv_12,dxh_22,dxu_22,dxv_22,dxh_32,dxu_32,dxv_32,
                       dxh_13,dxu_13,dxv_13,dxh_23,dxu_23,dxv_23,dxh_33,dxu_33,dxv_33,
                       dyh_11,dyu_11,dyv_11,dyh_21,dyu_21,dyv_21,dyh_31,dyu_31,dyv_31,
                       dyh_12,dyu_12,dyv_12,dyh_22,dyu_22,dyv_22,dyh_32,dyu_32,dyv_32,
                       dyh_13,dyu_13,dyv_13,dyh_23,dyu_23,dyv_23,dyh_33,dyu_33,dyv_33,
                       fh_1,fh_2,qh_1,qh_2,fu_1,fu_2,qu_1,qu_2,fv_1,fv_2,qv_1,qv_2)

