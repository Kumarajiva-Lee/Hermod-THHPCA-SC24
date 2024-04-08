#ifndef PHYSICAL_VARIABLE_H_INCLUDED
#define PHYSICAL_VARIABLE_H_INCLUDED 1

#include<stdbool.h>
struct HybridMeshField{ 
  double*** x_1;
  double*** x_2;
  double*** x_3;
  double*** y_1;
  double*** y_2;
  double*** y_3;
};
struct HybridStateField{ 
  double*** h_11;
  double*** u_11;
  double*** v_11;
  double*** jab_11;
  double*** jab5_11;
  double*** jab6_11;
  double*** jab7_11;
  double*** uc_11;
  double*** vc_11;
  double*** h_12;
  double*** u_12;
  double*** v_12;
  double*** jab_12;
  double*** jab5_12;
  double*** jab6_12;
  double*** jab7_12;
  double*** uc_12;
  double*** vc_12;
  double*** h_13;
  double*** u_13;
  double*** v_13;
  double*** jab_13;
  double*** jab5_13;
  double*** jab6_13;
  double*** jab7_13;
  double*** uc_13;
  double*** vc_13;
  double*** h_21;
  double*** u_21;
  double*** v_21;
  double*** jab_21;
  double*** jab5_21;
  double*** jab6_21;
  double*** jab7_21;
  double*** uc_21;
  double*** vc_21;
  double*** h_22;
  double*** u_22;
  double*** v_22;
  double*** jab_22;
  double*** jab5_22;
  double*** jab6_22;
  double*** jab7_22;
  double*** uc_22;
  double*** vc_22;
  double*** h_23;
  double*** u_23;
  double*** v_23;
  double*** jab_23;
  double*** jab5_23;
  double*** jab6_23;
  double*** jab7_23;
  double*** uc_23;
  double*** vc_23;
  double*** h_31;
  double*** u_31;
  double*** v_31;
  double*** jab_31;
  double*** jab5_31;
  double*** jab6_31;
  double*** jab7_31;
  double*** uc_31;
  double*** vc_31;
  double*** h_32;
  double*** u_32;
  double*** v_32;
  double*** jab_32;
  double*** jab5_32;
  double*** jab6_32;
  double*** jab7_32;
  double*** uc_32;
  double*** vc_32;
  double*** h_33;
  double*** u_33;
  double*** v_33;
  double*** jab_33;
  double*** jab5_33;
  double*** jab6_33;
  double*** jab7_33;
  double*** uc_33;
  double*** vc_33;
};
struct Vector2{
      double x;
      double y;
};
struct HybridTendField{ 
  double*** dxh_11;
  double*** dxu_11;
  double*** dxv_11;
  double*** dxh_21;
  double*** dxu_21;
  double*** dxv_21;
  double*** dxh_31;
  double*** dxu_31;
  double*** dxv_31;
  double*** dxh_12;
  double*** dxu_12;
  double*** dxv_12;
  double*** dxh_22;
  double*** dxu_22;
  double*** dxv_22;
  double*** dxh_32;
  double*** dxu_32;
  double*** dxv_32;
  double*** dxh_13;
  double*** dxu_13;
  double*** dxv_13;
  double*** dxh_23;
  double*** dxu_23;
  double*** dxv_23;
  double*** dxh_33;
  double*** dxu_33;
  double*** dxv_33;
  double*** dyh_11;
  double*** dyu_11;
  double*** dyv_11;
  double*** dyh_21;
  double*** dyu_21;
  double*** dyv_21;
  double*** dyh_31;
  double*** dyu_31;
  double*** dyv_31;
  double*** dyh_12;
  double*** dyu_12;
  double*** dyv_12;
  double*** dyh_22;
  double*** dyu_22;
  double*** dyv_22;
  double*** dyh_32;
  double*** dyu_32;
  double*** dyv_32;
  double*** dyh_13;
  double*** dyu_13;
  double*** dyv_13;
  double*** dyh_23;
  double*** dyu_23;
  double*** dyv_23;
  double*** dyh_33;
  double*** dyu_33;
  double*** dyv_33;
  double*** fh_1;
  double*** fh_2;
  double*** qh_1;
  double*** qh_2;
  double*** fu_1;
  double*** fu_2;
  double*** qu_1;
  double*** qu_2;
  double*** fv_1;
  double*** fv_2;
  double*** qv_1;
  double*** qv_2;
};
struct MeshVector3{
  int x;
  int y;
  int z;
};
struct MeshVector6{
  int lx;
  int ly;
  int lz;
  int sx;
  int sy;
  int sz;
};
struct Vector4{
  double a11;
  double a12;
  double a21;
  double a22;
};

struct HybridMeshField mesh[1];
struct HybridMeshField global_mesh[1];
struct HybridStateField state[3];
struct HybridTendField tend[3];

#define async_h_11 1
#define async_u_11 2
#define async_v_11 3
#define async_jab_11 4
#define async_jab5_11 5
#define async_jab6_11 6
#define async_jab7_11 7
#define async_uc_11 8
#define async_vc_11 9
#define async_h_12 10
#define async_u_12 11
#define async_v_12 12
#define async_jab_12 13
#define async_jab5_12 14
#define async_jab6_12 15
#define async_jab7_12 16
#define async_uc_12 17
#define async_vc_12 18
#define async_h_13 19
#define async_u_13 20
#define async_v_13 21
#define async_jab_13 22
#define async_jab5_13 23
#define async_jab6_13 24
#define async_jab7_13 25
#define async_uc_13 26
#define async_vc_13 27
#define async_h_21 28
#define async_u_21 29
#define async_v_21 30
#define async_jab_21 31
#define async_jab5_21 32
#define async_jab6_21 33
#define async_jab7_21 34
#define async_uc_21 35
#define async_vc_21 36
#define async_h_22 37
#define async_u_22 38
#define async_v_22 39
#define async_jab_22 40
#define async_jab5_22 41
#define async_jab6_22 42
#define async_jab7_22 43
#define async_uc_22 44
#define async_vc_22 45
#define async_h_23 46
#define async_u_23 47
#define async_v_23 48
#define async_jab_23 49
#define async_jab5_23 50
#define async_jab6_23 51
#define async_jab7_23 52
#define async_uc_23 53
#define async_vc_23 54
#define async_h_31 55
#define async_u_31 56
#define async_v_31 57
#define async_jab_31 58
#define async_jab5_31 59
#define async_jab6_31 60
#define async_jab7_31 61
#define async_uc_31 62
#define async_vc_31 63
#define async_h_32 64
#define async_u_32 65
#define async_v_32 66
#define async_jab_32 67
#define async_jab5_32 68
#define async_jab6_32 69
#define async_jab7_32 70
#define async_uc_32 71
#define async_vc_32 72
#define async_h_33 73
#define async_u_33 74
#define async_v_33 75
#define async_jab_33 76
#define async_jab5_33 77
#define async_jab6_33 78
#define async_jab7_33 79
#define async_uc_33 80
#define async_vc_33 81
#define async_dxh_11 82
#define async_dxu_11 83
#define async_dxv_11 84
#define async_dxh_21 85
#define async_dxu_21 86
#define async_dxv_21 87
#define async_dxh_31 88
#define async_dxu_31 89
#define async_dxv_31 90
#define async_dxh_12 91
#define async_dxu_12 92
#define async_dxv_12 93
#define async_dxh_22 94
#define async_dxu_22 95
#define async_dxv_22 96
#define async_dxh_32 97
#define async_dxu_32 98
#define async_dxv_32 99
#define async_dxh_13 100
#define async_dxu_13 101
#define async_dxv_13 102
#define async_dxh_23 103
#define async_dxu_23 104
#define async_dxv_23 105
#define async_dxh_33 106
#define async_dxu_33 107
#define async_dxv_33 108
#define async_dyh_11 109
#define async_dyu_11 110
#define async_dyv_11 111
#define async_dyh_21 112
#define async_dyu_21 113
#define async_dyv_21 114
#define async_dyh_31 115
#define async_dyu_31 116
#define async_dyv_31 117
#define async_dyh_12 118
#define async_dyu_12 119
#define async_dyv_12 120
#define async_dyh_22 121
#define async_dyu_22 122
#define async_dyv_22 123
#define async_dyh_32 124
#define async_dyu_32 125
#define async_dyv_32 126
#define async_dyh_13 127
#define async_dyu_13 128
#define async_dyv_13 129
#define async_dyh_23 130
#define async_dyu_23 131
#define async_dyv_23 132
#define async_dyh_33 133
#define async_dyu_33 134
#define async_dyv_33 135
#define async_fh_1 136
#define async_fh_2 137
#define async_qh_1 138
#define async_qh_2 139
#define async_fu_1 140
#define async_fu_2 141
#define async_qu_1 142
#define async_qu_2 143
#define async_fv_1 144
#define async_fv_2 145
#define async_qv_1 146
#define async_qv_2 147

extern double TimeBeg;
extern double TimeEnd;
extern double CommTime;
extern double CompTime;
extern double MeshInitTime;
extern double StateInitialh_11Time;
extern double StateInitialh_21Time;
extern double StateInitialh_31Time;
extern double StateInitialh_12Time;
extern double StateInitialh_22Time;
extern double StateInitialh_32Time;
extern double StateInitialh_13Time;
extern double StateInitialh_23Time;
extern double StateInitialh_33Time;
extern double updateStateh_11Time;
extern double updateStateh_12Time;
extern double updateStateh_13Time;
extern double updateStateh_21Time;
extern double updateStateh_22Time;
extern double updateStateh_23Time;
extern double updateStateh_31Time;
extern double updateStateh_32Time;
extern double updateStateh_33Time;
extern double updateXfh_1Time;
extern double updateXfh_2Time;
extern double updateXfu_2Time;
extern double updateXfv_2Time;
extern double updateXdxh_11Time;
extern double updateXdxh_12Time;
extern double updateXdxh_13Time;
extern double updateYfh_1Time;
extern double updateYfh_2Time;
extern double updateYfu_2Time;
extern double updateYfv_2Time;
extern double updateYdyh_11Time;
extern double update_rk1h_11Time;
extern double update_rk1h_12Time;
extern double update_rk1h_13Time;
extern double update_rk1h_21Time;
extern double update_rk1h_22Time;
extern double update_rk1h_23Time;
extern double update_rk1h_31Time;
extern double update_rk1h_32Time;
extern double update_rk1h_33Time;
extern double update_rk2h_11Time;
extern double update_rk2h_12Time;
extern double update_rk2h_13Time;
extern double update_rk2h_21Time;
extern double update_rk2h_22Time;
extern double update_rk2h_23Time;
extern double update_rk2h_31Time;
extern double update_rk2h_32Time;
extern double update_rk2h_33Time;
extern double update_rk3h_11Time;
extern double update_rk3h_12Time;
extern double update_rk3h_13Time;
extern double update_rk3h_21Time;
extern double update_rk3h_22Time;
extern double update_rk3h_23Time;
extern double update_rk3h_31Time;
extern double update_rk3h_32Time;
extern double update_rk3h_33Time;

void ProfilingOutPut(int id);

void PhysicalVariableInit();
void PhysicalVariableFinish();

void MeshInit_cp(struct HybridMeshField* global_mesh, struct HybridMeshField* mesh);
#endif
