typedef struct{
  int fnx,fny,fnz;
  int hnx,hny,hnz;
  int ghx,ghy,ghz;
}global_info;





typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_0_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *stateph;
  double *stateqv;
  double *statet;
} calc_t_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *advptuu;
} accum_uv_cell_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *advptvv;
} accum_uv_cell_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
  int nstep;
} accum_uv_cell_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
  int nstep;
} accum_uv_cell_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
} accum_uv_cell_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
} accum_uv_cell_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statetmpsum;
} accum_uv_cell_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *advptdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_0_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lon;
  double *stateu_lon;
  double *statemfx_lon;
} calc_mf_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lat;
  double *statev_lat;
  double *statemfy_lat;
} calc_mf_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *advptmfx;
} accum_mf_cell_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *advptmfy;
} accum_mf_cell_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
  int nstep;
} accum_mf_cell_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
  int nstep;
} accum_mf_cell_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
} accum_mf_cell_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
} accum_mf_cell_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *statemfx_lat;
  double *statem_lat;
  double *stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *statemfy_lon;
  double *statem_lon;
  double *statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *statetmpsum;
} calc_ke_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *stateke;
} calc_ke_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lat;
  double *statetmpsum;
  double *le_lat;
} calc_vor_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *statevor;
  double *area_cell;
} calc_vor_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statevor;
  double *statem_vtx;
  double *statepv;
  double *half_f;
} calc_pv_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lon;
  double *stateu_lon;
  double *statepv;
  double *statepv_lon;
} interp_pv_upwind_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lat;
  double *statev_lat;
  double *statepv;
  double *statepv_lat;
} interp_pv_upwind_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *statetmpsum;
} calc_div_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *statediv;
  double *le_lat;
  double *area_cell;
} calc_div_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statet;
  double *stateph_lev;
  double *staticvgzs;
  double *stategz_lev;
} calc_gz_lev_0_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advmold_m;
} copy_old_m_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advptold_m;
} copy_old_m_0_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_0_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_0_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *advptuu;
} accum_uv_cell_0_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *advptvv;
} accum_uv_cell_0_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
  int nstep;
} accum_uv_cell_0_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
  int nstep;
} accum_uv_cell_0_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_0_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_0_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
} accum_uv_cell_0_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
} accum_uv_cell_0_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_0_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_0_56_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_0_57_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statetmpsum;
} accum_uv_cell_0_58_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *advptdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_0_59_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lon;
  double *star_stateu_lon;
  double *star_statemfx_lon;
} calc_mf_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lat;
  double *star_statev_lat;
  double *star_statemfy_lat;
} calc_mf_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *advptmfx;
} accum_mf_cell_0_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *advptmfy;
} accum_mf_cell_0_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
  int nstep;
} accum_mf_cell_0_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
  int nstep;
} accum_mf_cell_0_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
} accum_mf_cell_0_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
} accum_mf_cell_0_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statemfx_lat;
  double *star_statem_lat;
  double *star_stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_0_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statemfy_lon;
  double *star_statem_lon;
  double *star_statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_0_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_stateu_lon;
  double *star_stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_statetmpsum;
} calc_ke_0_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_stateke;
} calc_ke_0_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_statetmpsum;
} calc_div_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_statediv;
  double *le_lat;
  double *area_cell;
} calc_div_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lat;
  double *star_statetmpsum;
  double *le_lat;
} calc_vor_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_statevor;
  double *area_cell;
} calc_vor_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statevor;
  double *star_statem_vtx;
  double *star_statepv;
  double *half_f;
} calc_pv_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lon;
  double *star_stateu_lon;
  double *star_statepv;
  double *star_statepv_lon;
} interp_pv_upwind_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lat;
  double *star_statev_lat;
  double *star_statepv;
  double *star_statepv_lat;
} interp_pv_upwind_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *tend1dmfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_mf_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *tend1dmfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_mf_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statetmpsum;
} calc_grad_mf_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *tend1dmfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_mf_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe0;
  double *advptm0;
  double *advptwe;
  double *advptmm;
} accum_we_lev_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statem_lev;
  double *advptwe;
  double *advptmm;
} accum_we_lev_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
  double *advptwe0;
  double *advptm0;
  int nstep;
} accum_we_lev_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
} accum_we_lev_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *advptmm;
  double *advptcflz;
  double dt;
} accum_we_lev_0_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lon;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lon;
} interp_lev_edge_to_lev_lon_edge_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lat;
  double *area_lat_north;
  double *area_lat_south;
  double *area_lat;
} interp_lev_edge_to_lev_lat_edge_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_inner_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *star_statept;
  double *advptuu;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_inner_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptvv;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_inner_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *star_stateptf_lon;
  double *advptdivx;
  double *star_stateptf_lat;
  double *advptdivy;
  double *advptqx;
  double *advptqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *star_statetmpsum;
} ffsl_calc_tracer_hflx_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *star_statetmpsum;
  double *advptdivy;
  double *advptqx;
  double *advptqy;
  double *le_lat;
  double *area_cell;
  double dt;
} ffsl_calc_tracer_hflx_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptqy;
  double *advptqx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_outer_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *advptqy;
  double *advptmfx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_outer_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptmfy;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_outer_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lon;
  double *tend1dptfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_ptf_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *tend1dptfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_ptf_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *star_statetmpsum;
} calc_grad_ptf_0_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *tend1dptfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_ptf_0_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
} calc_grad_ptf_0_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lev;
  double *tend1dptfdlev;
} calc_grad_ptf_0_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statepv_lat;
  double *star_statepv_lon;
  double *tend1qhu;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_coriolis_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statepv_lon;
  double *star_statepv_lat;
  double *tend1qhv;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_coriolis_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlon;
  double *de_lon;
} calc_grad_ke_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlat;
  double *de_lat;
} calc_grad_ke_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhv;
  double *tend1dkedlon;
  double *tend1wedudlev;
  double *tend1du;
} calc_tend_forward_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhu;
  double *tend1dkedlat;
  double *tend1wedvdlev;
  double *tend1dv;
} calc_tend_forward_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dptfdlon;
  double *tend1dptfdlat;
  double *tend1dptfdlev;
  double *tend1dpt;
} calc_tend_forward_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_0_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_0_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_0_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_0_36_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_37_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_0_38_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_39_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend1dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_0_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend1dgz;
  double *new_stategz;
  double dt;
} update_state_0_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_0_40_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_41_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_0_42_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_43_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend1du;
  double *new_stateu_lon;
  double dt;
} update_state_0_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend1dv;
  double *new_statev_lat;
  double dt;
} update_state_0_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statept;
  double *new_stateph;
  double *new_stateqv;
  double *new_statet;
} calc_t_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statet;
  double *new_stateph_lev;
  double *staticvgzs;
  double *new_stategz_lev;
} calc_gz_lev_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateqm;
  double *new_stateph_exn_lev;
  double *new_stategz_lev;
  double *tend2pgf_lon;
  double *de_lon;
} pgf_lin97_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1du;
  double *tend2pgf_lon;
  double *tend2du;
} calc_tend_backward_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dv;
  double *tend2pgf_lat;
  double *tend2dv;
} calc_tend_backward_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_0_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_0_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_0_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_0_44_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_0_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend2dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_0_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend2dgz;
  double *new_stategz;
  double dt;
} update_state_0_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_0_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_0_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_0_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend2du;
  double *new_stateu_lon;
  double dt;
} update_state_0_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend2dv;
  double *new_statev_lat;
  double dt;
} update_state_0_23_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statesmag_t;
  double *de_lon;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *statesmag_s;
  double *le_lat;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lon;
  double *de_lon;
  double *le_lon;
} smag_damp_run_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lat;
  double *le_lat;
  double *de_lat;
} smag_damp_run_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lon;
  double *stateu_lon;
  double *tendsmag_dudt;
  double *de_lon;
  double *de_lat;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *tendsmag_dudt;
  double dt;
} smag_damp_run_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lat;
  double *statev_lat;
  double *tendsmag_dvdt;
  double *le_lat;
} smag_damp_run_0_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lat;
  double *statev_lat;
  double *tendsmag_dvdt;
  double *le_lat;
  double *le_lon;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_0_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *tendsmag_dvdt;
  double dt;
} smag_damp_run_0_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} trickypt_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_0_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_0_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_0_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_0_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_0_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_0_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_0_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_0_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} pole_damp_run_0_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_0_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advu0;
  double *advuu;
} accum_uv_cell_0_60_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advv0;
  double *advvv;
} accum_uv_cell_0_61_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldu_lon;
  double *advuu;
} accum_uv_cell_0_62_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldv_lat;
  double *advvv;
} accum_uv_cell_0_63_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
  int nstep;
} accum_uv_cell_0_64_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
  int nstep;
} accum_uv_cell_0_65_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advu0;
} accum_uv_cell_0_66_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advv0;
} accum_uv_cell_0_67_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
} accum_uv_cell_0_68_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
} accum_uv_cell_0_69_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_0_70_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_0_71_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advvv;
  double *advdivx;
  double *advdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_0_72_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldtmpsum;
} accum_uv_cell_0_73_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldtmpsum;
  double *advdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_0_74_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfx_lon;
  double *advmfx;
} accum_mf_cell_0_24_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfy_lat;
  double *advmfy;
} accum_mf_cell_0_25_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
  int nstep;
} accum_mf_cell_0_26_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
  int nstep;
} accum_mf_cell_0_27_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
} accum_mf_cell_0_28_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
} accum_mf_cell_0_29_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe0;
  double *advm0;
  double *advwe;
  double *advmm;
} accum_we_lev_0_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldwe_lev;
  double *state_oldm_lev;
  double *advwe;
  double *advmm;
} accum_we_lev_0_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
  double *advwe0;
  double *advm0;
  int nstep;
} accum_we_lev_0_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
} accum_we_lev_0_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *advmm;
  double *advcflz;
  double dt;
} accum_we_lev_0_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_inner_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *state_oldqv;
  double *advuu;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_inner_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advvv;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_inner_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqmf_lon;
  double *advdivx;
  double *advqmf_lat;
  double *advdivy;
  double *advqx;
  double *advqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqmf_lat;
  double *state_oldtmpsum;
} ffsl_calc_tracer_hflx_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *state_oldtmpsum;
  double *advdivy;
  double *advqx;
  double *advqy;
  double *le_lat;
  double *area_cell;
  double dt;
} ffsl_calc_tracer_hflx_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqy;
  double *advqx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_outer_0_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *advqy;
  double *advmfx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_outer_0_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advmfy;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_outer_0_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advold_m;
  double *state_oldqv;
  double *advqmf_lon;
  double *advqmf_lat;
  double *state_newqv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} adv_run_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqmf_lat;
  double *state_newtmpsum;
} adv_run_0_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advold_m;
  double *state_oldqv;
  double *state_newtmpsum;
  double *state_newqv;
  double *le_lat;
  double *area_cell;
} adv_run_0_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_0_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
  double *advqmf_lev;
} adv_run_0_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_0_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldm;
  double *advold_m;
} copy_old_m_0_2_info;






typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_1_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *stateph;
  double *stateqv;
  double *statet;
} calc_t_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *advptuu;
} accum_uv_cell_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *advptvv;
} accum_uv_cell_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
  int nstep;
} accum_uv_cell_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
  int nstep;
} accum_uv_cell_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
} accum_uv_cell_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
} accum_uv_cell_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lon;
  double *stateu_lon;
  double *statemfx_lon;
} calc_mf_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lat;
  double *statev_lat;
  double *statemfy_lat;
} calc_mf_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *advptmfx;
} accum_mf_cell_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *advptmfy;
} accum_mf_cell_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
  int nstep;
} accum_mf_cell_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
  int nstep;
} accum_mf_cell_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
} accum_mf_cell_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
} accum_mf_cell_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *statemfx_lat;
  double *statem_lat;
  double *stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *statemfy_lon;
  double *statem_lon;
  double *statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statevor;
  double *statem_vtx;
  double *statepv;
  double *half_f;
} calc_pv_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lon;
  double *stateu_lon;
  double *statepv;
  double *statepv_lon;
} interp_pv_upwind_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lat;
  double *statev_lat;
  double *statepv;
  double *statepv_lat;
} interp_pv_upwind_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statet;
  double *stateph_lev;
  double *staticvgzs;
  double *stategz_lev;
} calc_gz_lev_1_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advmold_m;
} copy_old_m_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advptold_m;
} copy_old_m_1_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_1_39_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_1_40_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *advptuu;
} accum_uv_cell_1_41_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *advptvv;
} accum_uv_cell_1_42_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
  int nstep;
} accum_uv_cell_1_43_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
  int nstep;
} accum_uv_cell_1_44_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_1_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_1_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
} accum_uv_cell_1_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
} accum_uv_cell_1_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_1_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_1_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_1_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lon;
  double *star_stateu_lon;
  double *star_statemfx_lon;
} calc_mf_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lat;
  double *star_statev_lat;
  double *star_statemfy_lat;
} calc_mf_1_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *advptmfx;
} accum_mf_cell_1_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *advptmfy;
} accum_mf_cell_1_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
  int nstep;
} accum_mf_cell_1_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
  int nstep;
} accum_mf_cell_1_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
} accum_mf_cell_1_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
} accum_mf_cell_1_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statemfx_lat;
  double *star_statem_lat;
  double *star_stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_1_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statemfy_lon;
  double *star_statem_lon;
  double *star_statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_1_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_stateu_lon;
  double *star_stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statevor;
  double *star_statem_vtx;
  double *star_statepv;
  double *half_f;
} calc_pv_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lon;
  double *star_stateu_lon;
  double *star_statepv;
  double *star_statepv_lon;
} interp_pv_upwind_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lat;
  double *star_statev_lat;
  double *star_statepv;
  double *star_statepv_lat;
} interp_pv_upwind_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *tend1dmfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_mf_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *tend1dmfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_mf_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe0;
  double *advptm0;
  double *advptwe;
  double *advptmm;
} accum_we_lev_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statem_lev;
  double *advptwe;
  double *advptmm;
} accum_we_lev_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
  double *advptwe0;
  double *advptm0;
  int nstep;
} accum_we_lev_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
} accum_we_lev_1_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *advptmm;
  double *advptcflz;
  double dt;
} accum_we_lev_1_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lon;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lon;
} interp_lev_edge_to_lev_lon_edge_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lat;
  double *area_lat_north;
  double *area_lat_south;
  double *area_lat;
} interp_lev_edge_to_lev_lat_edge_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_inner_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *star_statept;
  double *advptuu;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_inner_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptvv;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_inner_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *star_stateptf_lon;
  double *advptdivx;
  double *star_stateptf_lat;
  double *advptdivy;
  double *advptqx;
  double *advptqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptqy;
  double *advptqx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_outer_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *advptqy;
  double *advptmfx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_outer_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptmfy;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_outer_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lon;
  double *tend1dptfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_ptf_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *tend1dptfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_ptf_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
} calc_grad_ptf_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lev;
  double *tend1dptfdlev;
} calc_grad_ptf_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statepv_lat;
  double *star_statepv_lon;
  double *tend1qhu;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_coriolis_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statepv_lon;
  double *star_statepv_lat;
  double *tend1qhv;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_coriolis_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlon;
  double *de_lon;
} calc_grad_ke_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlat;
  double *de_lat;
} calc_grad_ke_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhv;
  double *tend1dkedlon;
  double *tend1wedudlev;
  double *tend1du;
} calc_tend_forward_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhu;
  double *tend1dkedlat;
  double *tend1wedvdlev;
  double *tend1dv;
} calc_tend_forward_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dptfdlon;
  double *tend1dptfdlat;
  double *tend1dptfdlev;
  double *tend1dpt;
} calc_tend_forward_1_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_1_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_1_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_1_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_1_36_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_37_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_1_38_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_39_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend1dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_1_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend1dgz;
  double *new_stategz;
  double dt;
} update_state_1_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_1_40_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_41_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_1_42_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_43_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend1du;
  double *new_stateu_lon;
  double dt;
} update_state_1_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend1dv;
  double *new_statev_lat;
  double dt;
} update_state_1_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statept;
  double *new_stateph;
  double *new_stateqv;
  double *new_statet;
} calc_t_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statet;
  double *new_stateph_lev;
  double *staticvgzs;
  double *new_stategz_lev;
} calc_gz_lev_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateqm;
  double *new_stateph_exn_lev;
  double *new_stategz_lev;
  double *tend2pgf_lon;
  double *de_lon;
} pgf_lin97_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1du;
  double *tend2pgf_lon;
  double *tend2du;
} calc_tend_backward_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dv;
  double *tend2pgf_lat;
  double *tend2dv;
} calc_tend_backward_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_1_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_1_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_1_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_1_44_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_1_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend2dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_1_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend2dgz;
  double *new_stategz;
  double dt;
} update_state_1_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_1_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_1_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_1_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend2du;
  double *new_stateu_lon;
  double dt;
} update_state_1_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend2dv;
  double *new_statev_lat;
  double dt;
} update_state_1_23_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statesmag_t;
  double *de_lon;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *statesmag_s;
  double *le_lat;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lon;
  double *de_lon;
  double *le_lon;
} smag_damp_run_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lat;
  double *le_lat;
  double *de_lat;
} smag_damp_run_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lon;
  double *stateu_lon;
  double *tendsmag_dudt;
  double *de_lon;
  double *de_lat;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_1_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *tendsmag_dudt;
  double dt;
} smag_damp_run_1_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lat;
  double *statev_lat;
  double *tendsmag_dvdt;
  double *le_lat;
  double *le_lon;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_1_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *tendsmag_dvdt;
  double dt;
} smag_damp_run_1_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} trickypt_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_1_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_1_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_1_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_1_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_1_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_1_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_1_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_1_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_1_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_1_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} pole_damp_run_1_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_1_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advu0;
  double *advuu;
} accum_uv_cell_1_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advv0;
  double *advvv;
} accum_uv_cell_1_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldu_lon;
  double *advuu;
} accum_uv_cell_1_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldv_lat;
  double *advvv;
} accum_uv_cell_1_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
  int nstep;
} accum_uv_cell_1_56_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
  int nstep;
} accum_uv_cell_1_57_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advu0;
} accum_uv_cell_1_58_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advv0;
} accum_uv_cell_1_59_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
} accum_uv_cell_1_60_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
} accum_uv_cell_1_61_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_1_62_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_1_63_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advvv;
  double *advdivx;
  double *advdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_1_64_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfx_lon;
  double *advmfx;
} accum_mf_cell_1_24_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfy_lat;
  double *advmfy;
} accum_mf_cell_1_25_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
  int nstep;
} accum_mf_cell_1_26_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
  int nstep;
} accum_mf_cell_1_27_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
} accum_mf_cell_1_28_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
} accum_mf_cell_1_29_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe0;
  double *advm0;
  double *advwe;
  double *advmm;
} accum_we_lev_1_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldwe_lev;
  double *state_oldm_lev;
  double *advwe;
  double *advmm;
} accum_we_lev_1_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
  double *advwe0;
  double *advm0;
  int nstep;
} accum_we_lev_1_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
} accum_we_lev_1_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *advmm;
  double *advcflz;
  double dt;
} accum_we_lev_1_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_inner_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *state_oldqv;
  double *advuu;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_inner_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advvv;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_inner_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqmf_lon;
  double *advdivx;
  double *advqmf_lat;
  double *advdivy;
  double *advqx;
  double *advqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqy;
  double *advqx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_outer_1_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *advqy;
  double *advmfx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_outer_1_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advmfy;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_outer_1_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advold_m;
  double *state_oldqv;
  double *advqmf_lon;
  double *advqmf_lat;
  double *state_newqv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} adv_run_1_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_1_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
  double *advqmf_lev;
} adv_run_1_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_1_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldm;
  double *advold_m;
} copy_old_m_1_2_info;






typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_2_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *stateph;
  double *stateqv;
  double *statet;
} calc_t_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *advptuu;
} accum_uv_cell_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *advptvv;
} accum_uv_cell_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
  int nstep;
} accum_uv_cell_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
  int nstep;
} accum_uv_cell_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *stateu_lon;
} accum_uv_cell_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statev_lat;
} accum_uv_cell_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *statetmpsum;
} accum_uv_cell_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *advptdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_2_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lon;
  double *stateu_lon;
  double *statemfx_lon;
} calc_mf_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem_lat;
  double *statev_lat;
  double *statemfy_lat;
} calc_mf_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *advptmfx;
} accum_mf_cell_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *advptmfy;
} accum_mf_cell_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
  int nstep;
} accum_mf_cell_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
  int nstep;
} accum_mf_cell_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *statemfx_lon;
} accum_mf_cell_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *statemfy_lat;
} accum_mf_cell_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfx_lon;
  double *statemfx_lat;
  double *statem_lat;
  double *stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statemfy_lat;
  double *statemfy_lon;
  double *statem_lon;
  double *statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *statetmpsum;
} calc_ke_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *stateke;
} calc_ke_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lat;
  double *statetmpsum;
  double *le_lat;
} calc_vor_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *statevor;
  double *area_cell;
} calc_vor_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statevor;
  double *statem_vtx;
  double *statepv;
  double *half_f;
} calc_pv_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lon;
  double *stateu_lon;
  double *statepv;
  double *statepv_lon;
} interp_pv_upwind_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lat;
  double *statev_lat;
  double *statepv;
  double *statepv_lat;
} interp_pv_upwind_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *statetmpsum;
} calc_div_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statetmpsum;
  double *statediv;
  double *le_lat;
  double *area_cell;
} calc_div_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statet;
  double *stateph_lev;
  double *staticvgzs;
  double *stategz_lev;
} calc_gz_lev_2_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advmold_m;
} copy_old_m_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *advptold_m;
} copy_old_m_2_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptu0;
  double *advptuu;
} accum_uv_cell_2_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptv0;
  double *advptvv;
} accum_uv_cell_2_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *advptuu;
} accum_uv_cell_2_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *advptvv;
} accum_uv_cell_2_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
  int nstep;
} accum_uv_cell_2_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
  int nstep;
} accum_uv_cell_2_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptu0;
} accum_uv_cell_2_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptv0;
} accum_uv_cell_2_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *star_stateu_lon;
} accum_uv_cell_2_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statev_lat;
} accum_uv_cell_2_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_2_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *advptcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_2_56_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptuu;
  double *advptvv;
  double *advptdivx;
  double *advptdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_2_57_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptvv;
  double *star_statetmpsum;
} accum_uv_cell_2_58_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *advptdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_2_59_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lon;
  double *star_stateu_lon;
  double *star_statemfx_lon;
} calc_mf_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statem_lat;
  double *star_statev_lat;
  double *star_statemfy_lat;
} calc_mf_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *advptmfx;
} accum_mf_cell_2_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *advptmfy;
} accum_mf_cell_2_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
  int nstep;
} accum_mf_cell_2_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
  int nstep;
} accum_mf_cell_2_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfx;
  double *star_statemfx_lon;
} accum_mf_cell_2_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptmfy;
  double *star_statemfy_lat;
} accum_mf_cell_2_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statemfx_lat;
  double *star_statem_lat;
  double *star_stateu_lat;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_mf_2_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statemfy_lon;
  double *star_statem_lon;
  double *star_statev_lon;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_mf_2_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_stateke;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lat_north;
  double *area_lat_south;
  double *area_cell;
} calc_ke_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_stateu_lon;
  double *star_stateke;
  double *area_lat_east;
  double *area_lat_west;
  double *area_lon_north;
  double *area_lon_south;
  double *area_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_cell;
} calc_ke_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_statetmpsum;
} calc_ke_2_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_stateke;
} calc_ke_2_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statediv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} calc_div_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lat;
  double *star_statetmpsum;
} calc_div_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_statediv;
  double *le_lat;
  double *area_cell;
} calc_div_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lon;
  double *star_statev_lat;
  double *star_statevor;
  double *de_lon;
  double *de_lat;
  double *area_vtx;
} calc_vor_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lat;
  double *star_statetmpsum;
  double *le_lat;
} calc_vor_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *star_statevor;
  double *area_cell;
} calc_vor_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statevor;
  double *star_statem_vtx;
  double *star_statepv;
  double *half_f;
} calc_pv_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statev_lon;
  double *star_stateu_lon;
  double *star_statepv;
  double *star_statepv_lon;
} interp_pv_upwind_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateu_lat;
  double *star_statev_lat;
  double *star_statepv;
  double *star_statepv_lat;
} interp_pv_upwind_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *tend1dmfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_mf_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *tend1dmfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_mf_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statetmpsum;
} calc_grad_mf_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *tend1dmfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_mf_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe0;
  double *advptm0;
  double *advptwe;
  double *advptmm;
} accum_we_lev_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statem_lev;
  double *advptwe;
  double *advptmm;
} accum_we_lev_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
  double *advptwe0;
  double *advptm0;
  int nstep;
} accum_we_lev_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *star_statewe_lev;
  double *advptmm;
  double *star_statem_lev;
} accum_we_lev_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptwe;
  double *advptmm;
  double *advptcflz;
  double dt;
} accum_we_lev_2_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lon;
  double *area_lon_west;
  double *area_lon_east;
  double *area_lon;
} interp_lev_edge_to_lev_lon_edge_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev;
  double *star_statewe_lev_lat;
  double *area_lat_north;
  double *area_lat_south;
  double *area_lat;
} interp_lev_edge_to_lev_lat_edge_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lon;
  double *star_stateu_lon;
  double *star_statem_lon;
  double *tend1wedudlev;
} calc_wedudlev_wedvdlev_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statewe_lev_lat;
  double *star_statev_lat;
  double *star_statem_lat;
  double *tend1wedvdlev;
} calc_wedudlev_wedvdlev_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_inner_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *star_statept;
  double *advptuu;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_inner_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptvv;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_inner_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *star_stateptf_lon;
  double *advptdivx;
  double *star_stateptf_lat;
  double *advptdivy;
  double *advptqx;
  double *advptqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *star_statetmpsum;
} ffsl_calc_tracer_hflx_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
  double *star_statetmpsum;
  double *advptdivy;
  double *advptqx;
  double *advptqy;
  double *le_lat;
  double *area_cell;
  double dt;
} ffsl_calc_tracer_hflx_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptqy;
  double *advptqx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
} hflx_ppm_outer_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcflx;
  double *advptqy;
  double *advptmfx;
  double *advptqlx;
  double *advptdqx;
  double *advptq6x;
  double *star_stateptf_lon;
} hflx_ppm_outer_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advptcfly;
  double *advptmfy;
  double *advptqly;
  double *advptdqy;
  double *advptq6y;
  double *star_stateptf_lat;
} hflx_ppm_outer_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lon;
  double *tend1dptfdlon;
  double *le_lon;
  double *area_cell;
} calc_grad_ptf_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *tend1dptfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_ptf_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lat;
  double *star_statetmpsum;
} calc_grad_ptf_2_14_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statetmpsum;
  double *tend1dptfdlat;
  double *le_lat;
  double *area_cell;
} calc_grad_ptf_2_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statept;
} calc_grad_ptf_2_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateptf_lev;
  double *tend1dptfdlev;
} calc_grad_ptf_2_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfx_lon;
  double *star_statepv_lat;
  double *star_statepv_lon;
  double *tend1qhu;
  double *half_tangent_wgt_0;
  double *half_tangent_wgt_1;
} calc_coriolis_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_statemfy_lat;
  double *star_statepv_lon;
  double *star_statepv_lat;
  double *tend1qhv;
  double *full_tangent_wgt_0;
  double *full_tangent_wgt_1;
} calc_coriolis_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlon;
  double *de_lon;
} calc_grad_ke_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *star_stateke;
  double *tend1dkedlat;
  double *de_lat;
} calc_grad_ke_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhv;
  double *tend1dkedlon;
  double *tend1wedudlev;
  double *tend1du;
} calc_tend_forward_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1qhu;
  double *tend1dkedlat;
  double *tend1wedvdlev;
  double *tend1dv;
} calc_tend_forward_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dptfdlon;
  double *tend1dptfdlat;
  double *tend1dptfdlev;
  double *tend1dpt;
} calc_tend_forward_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_2_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_2_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_2_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_2_36_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_37_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_2_38_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_39_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend1dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_2_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend1dgz;
  double *new_stategz;
  double dt;
} update_state_2_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_2_40_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_41_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_2_42_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_43_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend1du;
  double *new_stateu_lon;
  double dt;
} update_state_2_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend1dv;
  double *new_statev_lat;
  double dt;
} update_state_2_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statept;
  double *new_stateph;
  double *new_stateqv;
  double *new_statet;
} calc_t_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statet;
  double *new_stateph_lev;
  double *staticvgzs;
  double *new_stategz_lev;
} calc_gz_lev_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateqm;
  double *new_stateph_exn_lev;
  double *new_stategz_lev;
  double *tend2pgf_lon;
  double *de_lon;
} pgf_lin97_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1du;
  double *tend2pgf_lon;
  double *tend2du;
} calc_tend_backward_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *tend1dv;
  double *tend2pgf_lat;
  double *tend2dv;
} calc_tend_backward_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statephs;
  double *new_stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_2_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph_exn_lev;
} calc_ph_2_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
} calc_ph_2_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_2_44_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_45_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_2_46_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_47_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statept;
  double *old_statem;
  double *tend2dpt;
  double *new_statem;
  double *new_statept;
  double dt;
} update_state_2_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stategz;
  double *tend2dgz;
  double *new_stategz;
  double dt;
} update_state_2_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_statem;
} calc_m_2_48_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_49_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph;
  double *new_stateph_lev;
  double *new_statem_lev;
} calc_m_2_50_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_stateph_lev;
  double *new_stateph;
  double *new_statem_lev;
} calc_m_2_51_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lon;
} averageCellToLonEdge_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_lat;
} averageCellToLatEdge_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *new_statem;
  double *new_statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_12_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_stateu_lon;
  double *tend2du;
  double *new_stateu_lon;
  double dt;
} update_state_2_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *old_statev_lat;
  double *tend2dv;
  double *new_statev_lat;
  double dt;
} update_state_2_23_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *statesmag_t;
  double *de_lon;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *stateu_lon;
  double *statesmag_s;
  double *le_lat;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lon;
  double *de_lon;
  double *le_lon;
} smag_damp_run_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statesmag_t;
  double *statesmag_s;
  double *statekmh_lat;
  double *le_lat;
  double *de_lat;
} smag_damp_run_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lon;
  double *stateu_lon;
  double *tendsmag_dudt;
  double *de_lon;
  double *de_lat;
  double *half_cos_lat;
  double *le_lon;
  double *full_cos_lat;
} smag_damp_run_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *tendsmag_dudt;
  double dt;
} smag_damp_run_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lat;
  double *statev_lat;
  double *tendsmag_dvdt;
  double *le_lat;
  double *le_lon;
  double *full_cos_lat;
  double *de_lat;
  double *half_cos_lat;
} smag_damp_run_2_6_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statekmh_lat;
  double *statev_lat;
  double *tendsmag_dvdt;
  double *le_lat;
} smag_damp_run_2_7_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statev_lat;
  double *tendsmag_dvdt;
  double dt;
} smag_damp_run_2_8_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} trickypt_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statephs;
  double *stateph_lev;
  double *hyai;
  double *hybi;
} calc_ph_2_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph_exn_lev;
} calc_ph_2_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
} calc_ph_2_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *statem;
} calc_m_2_52_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *statem_lev;
} calc_m_2_53_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph;
  double *stateph_lev;
  double *statem_lev;
} calc_m_2_54_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateph_lev;
  double *stateph;
  double *statem_lev;
} calc_m_2_55_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lon;
} averageCellToLonEdge_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_lat;
} averageCellToLatEdge_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statem;
  double *statem_vtx;
  double *area_subcell_1;
  double *area_subcell_0;
  double *area_vtx;
} interpCellToVtx_2_13_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *statept;
  double *statem;
} pole_damp_run_2_0_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *stateu_lon;
  double *statev_lat;
  double *stateu;
  double *statev;
} c2a_2_1_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advu0;
  double *advuu;
} accum_uv_cell_2_60_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advv0;
  double *advvv;
} accum_uv_cell_2_61_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldu_lon;
  double *advuu;
} accum_uv_cell_2_62_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldv_lat;
  double *advvv;
} accum_uv_cell_2_63_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
  int nstep;
} accum_uv_cell_2_64_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
  int nstep;
} accum_uv_cell_2_65_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advu0;
} accum_uv_cell_2_66_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advv0;
} accum_uv_cell_2_67_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *state_oldu_lon;
} accum_uv_cell_2_68_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldv_lat;
} accum_uv_cell_2_69_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advcflx;
  double *de_lon;
  double dt;
} accum_uv_cell_2_70_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *advcfly;
  double *de_lat;
  double dt;
} accum_uv_cell_2_71_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advuu;
  double *advvv;
  double *advdivx;
  double *advdivy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
} accum_uv_cell_2_72_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advvv;
  double *state_oldtmpsum;
} accum_uv_cell_2_73_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldtmpsum;
  double *advdivy;
  double *le_lat;
  double *area_cell;
} accum_uv_cell_2_74_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfx_lon;
  double *advmfx;
} accum_mf_cell_2_24_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldmfy_lat;
  double *advmfy;
} accum_mf_cell_2_25_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
  int nstep;
} accum_mf_cell_2_26_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
  int nstep;
} accum_mf_cell_2_27_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfx;
  double *state_oldmfx_lon;
} accum_mf_cell_2_28_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advmfy;
  double *state_oldmfy_lat;
} accum_mf_cell_2_29_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe0;
  double *advm0;
  double *advwe;
  double *advmm;
} accum_we_lev_2_15_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldwe_lev;
  double *state_oldm_lev;
  double *advwe;
  double *advmm;
} accum_we_lev_2_16_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
  double *advwe0;
  double *advm0;
  int nstep;
} accum_we_lev_2_17_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *state_oldwe_lev;
  double *advmm;
  double *state_oldm_lev;
} accum_we_lev_2_18_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advwe;
  double *advmm;
  double *advcflz;
  double dt;
} accum_we_lev_2_19_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_inner_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *state_oldqv;
  double *advuu;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_inner_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advvv;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_inner_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *advqmf_lon;
  double *advdivx;
  double *advqmf_lat;
  double *advdivy;
  double *advqx;
  double *advqy;
  double *le_lon;
  double *area_cell;
  double *le_lat;
  double dt;
} ffsl_calc_tracer_hflx_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqmf_lat;
  double *state_oldtmpsum;
} ffsl_calc_tracer_hflx_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldqv;
  double *state_oldtmpsum;
  double *advdivy;
  double *advqx;
  double *advqy;
  double *le_lat;
  double *area_cell;
  double dt;
} ffsl_calc_tracer_hflx_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqy;
  double *advqx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqly;
  double *advdqy;
  double *advq6y;
} hflx_ppm_outer_2_9_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcflx;
  double *advqy;
  double *advmfx;
  double *advqlx;
  double *advdqx;
  double *advq6x;
  double *advqmf_lon;
} hflx_ppm_outer_2_10_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advcfly;
  double *advmfy;
  double *advqly;
  double *advdqy;
  double *advq6y;
  double *advqmf_lat;
} hflx_ppm_outer_2_11_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advold_m;
  double *state_oldqv;
  double *advqmf_lon;
  double *advqmf_lat;
  double *state_newqv;
  double *le_lon;
  double *le_lat;
  double *area_cell;
} adv_run_2_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advqmf_lat;
  double *state_newtmpsum;
} adv_run_2_1_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *advold_m;
  double *state_oldqv;
  double *state_newtmpsum;
  double *state_newqv;
  double *le_lat;
  double *area_cell;
} adv_run_2_2_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_2_3_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
  double *advqmf_lev;
} adv_run_2_4_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_newqv;
  double *state_oldm;
} adv_run_2_5_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_oldm;
  double *advold_m;
} copy_old_m_2_2_info;



