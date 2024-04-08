

#100
nlon = 360
nlat = 180
nlev = 32
dt_dyn = 240
dt_adv = 0
ncIn = "gmcore.restart.DSL.100.nc"

#50
# nlon = 720
# nlat = 360
# nlev = 32
# dt_dyn = 120
# dt_adv = 0
# ncIn = "gmcore.restart.DSL.50.nc"

#25
# nlon = 1440
# nlat = 720
# nlev = 32
# dt_dyn = 45
# dt_adv = 0 #90
# ncIn = "gmcore.restart.DSL.25.nc"

#15
# nlon = 2400
# nlat = 1200
# nlev = 32
# dt_dyn = 30
# dt_adv = 30
# ncIn = "gmcore.restart.DSL.15.nc"




planet = 'earth'

num_proc = 20

hydrostatic = True
vert_coord_template = 32
ptop = 2.194e2
p0 = 1.0e5

time_start = 0
time_end = 3600

tangent_wgt_scheme = 1 # 1 = 'classic'
ke_scheme   = 2
ke_cell_wgt = 0.5
pv_scheme   = 2
pv_pole_stokes = True
weno_order_pv = 3
upwind_order_pv = 3
upwind_wgt = 1.0
upwind_wgt_pv = 1

# Filter settings
max_wave_speed       = 300
max_cfl              = 0.5
filter_coef_a        = 1.0
filter_coef_b        = 0.4
filter_coef_c        = 0.2
filter_reset_interval= 0

pgf_scheme = 'lin97'
coriolis_scheme = 1

coarse_pole_mul = 0
coarse_pole_decay = 100

#damp setting
use_smag_damp        = True
div_damp_coef2       = 1.0 / 128.0

smag_damp_coef       = 0.1
