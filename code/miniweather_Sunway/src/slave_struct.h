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
  double *staterho;
  double *stateu;
  double *statew;
  double *statept;
  double *xp;
  double *zpf;
} StateInitial_0_0_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *staterho;
  double *stateu;
  double *statew;
  double *statept;
  double *state_tmprho;
  double *state_tmpu;
  double *state_tmpw;
  double *state_tmppt;
} StateInitial_0_1_info;


typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_forcingrho;
  double *state_forcingu;
  double *state_forcingw;
  double *state_forcingpt;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *hy_dens_cell;
  double *hy_dens_theta_cell;
  double hv_coef;
} ComputeTendenciesX_0_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *tenddrho;
  double *tenddu;
  double *tenddw;
  double *tenddpt;
} ComputeTendenciesX_0_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_forcingrho;
  double *state_forcingu;
  double *state_forcingw;
  double *state_forcingpt;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *hy_dens_int;
  double *hy_dens_theta_int;
  double *hy_pressure_int;
  double hv_coef;
} ComputeTendenciesZ_0_20_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *tenddw;
  double *state_forcingrho;
  double *fluxfpt;
  double *tenddrho;
  double *tenddu;
  double *tenddpt;
} ComputeTendenciesZ_0_21_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_forcingrho;
  double *state_forcingu;
  double *state_forcingw;
  double *state_forcingpt;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *hy_dens_int;
  double *hy_dens_theta_int;
  double *hy_pressure_int;
  double hv_coef;
} ComputeTendenciesZ_0_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *tenddw;
  double *state_forcingrho;
  double *fluxfpt;
  double *tenddrho;
  double *tenddu;
  double *tenddpt;
} ComputeTendenciesZ_0_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_forcingrho;
  double *state_forcingu;
  double *state_forcingw;
  double *state_forcingpt;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *hy_dens_cell;
  double *hy_dens_theta_cell;
  double hv_coef;
} ComputeTendenciesX_0_22_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *fluxfrho;
  double *fluxfu;
  double *fluxfw;
  double *fluxfpt;
  double *tenddrho;
  double *tenddu;
  double *tenddw;
  double *tenddpt;
} ComputeTendenciesX_0_23_info;

typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *state_initrho;
  double *tenddrho;
  double *state_initu;
  double *tenddu;
  double *state_initw;
  double *tenddw;
  double *state_initpt;
  double *tenddpt;
  double *state_outrho;
  double *state_outu;
  double *state_outw;
  double *state_outpt;
  double dt;
} UpdateState_0_5_info;



typedef struct{
  int lx, ly, lz;
  int ox, oy, oz;
  int hx, hy, hz;
  int bx, by, bz;
  int mx, my, mz;
  double *staterho;
  double *stateu;
  double *statew;
  double *statept;
  double *stateoutdens;
  double *stateoutuwnd;
  double *stateoutwwnd;
  double *stateouttheta;
  double *hy_dens_cell;
  double *hy_dens_theta_cell;
} OutputPrepare_0_0_info;


