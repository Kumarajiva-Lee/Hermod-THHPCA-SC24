import frontend.gmcore_smallball.mesh as gm
import frontend.gmcore_smallball.global_config as gc
from frontend.pytoy.lang import Float64, Bool, Int32, hybrid   
from frontend.pytoy.lang.dtype_ext import LonLatField 

@hybrid
class HybridStateField:

    u    : LonLatField[Float64]
    u_lon: LonLatField[Float64]
    u_lat: LonLatField[Float64]
    v    : LonLatField[Float64]
    v_lat: LonLatField[Float64]
    v_lon: LonLatField[Float64]

    we_lev : LonLatField[Float64]
    we_lev_lon : LonLatField[Float64]
    we_lev_lat : LonLatField[Float64]

    gz : LonLatField[Float64]
    gz_lev : LonLatField[Float64]

    m : LonLatField[Float64]
    m_vtx : LonLatField[Float64]
    m_lon : LonLatField[Float64]
    m_lat : LonLatField[Float64]
    m_lev : LonLatField[Float64]

    mfx_lon : LonLatField[Float64]
    mfy_lat : LonLatField[Float64]
    mfx_lat : LonLatField[Float64]
    mfy_lon : LonLatField[Float64]

    pv : LonLatField[Float64]
    pv_lon : LonLatField[Float64]
    pv_lat : LonLatField[Float64]

    ke : LonLatField[Float64]

    pt : LonLatField[Float64]
    ptf_lon : LonLatField[Float64]
    ptf_lat : LonLatField[Float64]
    ptf_lev : LonLatField[Float64]

    t : LonLatField[Float64]

    ph : LonLatField[Float64]
    ph_lev : LonLatField[Float64]
    ph_exn_lev : LonLatField[Float64]

    phs : LonLatField[Float64]

    div : LonLatField[Float64]
#     div2 : LonLatField[Float64]
    vor : LonLatField[Float64]

#     #Moist variables
    qv : LonLatField[Float64]
    qm : LonLatField[Float64]

#     #Smagorinsky damping variables
    smag_t : LonLatField[Float64]
    smag_s : LonLatField[Float64]
    kmh : LonLatField[Float64]
    kmh_lon : LonLatField[Float64]
    kmh_lat : LonLatField[Float64]

    #adv variables
    q : LonLatField[Float64]
    qmf_lon : LonLatField[Float64]
    qmf_lat : LonLatField[Float64]
    qmf_lev : LonLatField[Float64]

    #reduce_sum
    tmpsum : LonLatField[Float64]

u = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
u_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
u_lat = LonLatField(Float64,(gm.full_nlat,gm.half_nlat,gm.full_nlev),(True,False,True))
v = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
v_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
v_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))

we_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
we_lev_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.half_nlev),(False,True,False))
we_lev_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.half_nlev),(True,False,False))

gz = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
gz_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

m = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
m_vtx = LonLatField(Float64,(gm.half_nlon,gm.half_nlat,gm.full_nlev),(False,False,True))
m_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
m_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
m_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

mfx_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
mfy_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
mfy_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
mfx_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))

pv = LonLatField(Float64,(gm.half_nlon,gm.half_nlat,gm.full_nlev),(False,False,True))
pv_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
pv_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))

ke = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))

pt = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
ptf_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
ptf_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
ptf_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

t = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
ph = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
ph_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
ph_exn_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

phs = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))

div = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
vor = LonLatField(Float64,(gm.half_nlon,gm.half_nlat,gm.full_nlev),(False,False,True))

# div2 = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))

qv = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qm = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))

smag_t = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
smag_s = LonLatField(Float64,(gm.half_nlon,gm.half_nlat,gm.full_nlev),(False,False,True))
kmh = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
kmh_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
kmh_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))

q = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qmf_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
qmf_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True)) 
qmf_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

tmpsum = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
# state = HybridStateField(u_lon,v_lat,pt,phs)
state = HybridStateField(u, u_lon, u_lat, v, v_lat, v_lon,
                         we_lev, we_lev_lon, we_lev_lat,
                         gz, gz_lev,
                         m, m_vtx, m_lon, m_lat, m_lev,
                         mfx_lon, mfy_lat, mfx_lat, mfy_lon,
                         pv, pv_lon, pv_lat, ke,
                         pt,ptf_lon,ptf_lat,ptf_lev,
                         t, ph, ph_lev, ph_exn_lev, phs,
                         div, vor,
                         qv, qm,
                         smag_t, smag_s, kmh, kmh_lon, kmh_lat,
                         q, qmf_lon, qmf_lat, qmf_lev,
                         tmpsum)
# state = HybridStateField(u, u_lon,u_lat, v, v_lon, v_lat,
#                         we_lev, we_lev_lon, we_lev_lat,
#                         gz, gz_lev,
#                         m, m_vtx, m_lon, m_lat, m_lev,
#                         mfx_lon, mfy_lat, mfx_lat, mfy_lon,
#                         pv, pv_lon, pv_lat, ke,
#                         pt, ptf_lon, ptf_lat, ptf_lev,
#                         t, ph, ph_lev, ph_exn_lev, phs, div, vor, div2,
#                         qv, qm,
#                         smag_t, smag_s, kmh, kmh_lon, kmh_lat)

#ToRemeber: adv[0]为adv_pt , adv[1]为advs[1]
@hybrid
class HybridAdvField:
    old_m : LonLatField[Float64]
    mfx : LonLatField[Float64]
    mfy : LonLatField[Float64]
    mm  : LonLatField[Float64]
    m0 : LonLatField[Float64]
    uu : LonLatField[Float64]
    u0 : LonLatField[Float64]
    vv  : LonLatField[Float64]
    v0 : LonLatField[Float64]
    we : LonLatField[Float64]
    we0 : LonLatField[Float64]
    cflx : LonLatField[Float64]
    cfly : LonLatField[Float64]
    cflz : LonLatField[Float64]
    divx : LonLatField[Float64]
    divy : LonLatField[Float64]
    qlx : LonLatField[Float64]
    qly : LonLatField[Float64]
    dqx : LonLatField[Float64]
    dqy : LonLatField[Float64]
    q6x : LonLatField[Float64]
    q6y : LonLatField[Float64]
    qx : LonLatField[Float64]
    qy : LonLatField[Float64]
    qmf_lon : LonLatField[Float64]
    qmf_lat : LonLatField[Float64]
    qmf_lev : LonLatField[Float64]

old_m = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
mfx = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
mfy = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
mm = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
m0 = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
uu = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
u0 = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
vv = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
v0 = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
we = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
we0 = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
cflx = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
cfly = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
cflz = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))
divx = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
divy = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qlx = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qly = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dqx = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dqy = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
q6x = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
q6y = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qx = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qy = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
qmf_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
qmf_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
qmf_lev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.half_nlev),(True,True,False))

adv = HybridAdvField(old_m, mfx, mfy, mm, m0,
                    uu, u0, vv, v0, we, we0,
                    cflx, cfly, cflz, divx, divy,
                    qlx,qly,dqx,dqy,q6x,q6y,
                    qx,qy,
                    qmf_lon, qmf_lat, qmf_lev)

@hybrid
class HybridAdvPara:
    dynamic : Bool
    nstep : Int32
    uv_step : Int32
    we_step : Int32
    mf_step : Int32

advptPara = HybridAdvPara(True,1,0,0,0)
advmPara = HybridAdvPara(False,1,0,0,0)

@hybrid
class HybridTendField:
    du : LonLatField[Float64]
    dv : LonLatField[Float64]
    dgz : LonLatField[Float64]
    dpt : LonLatField[Float64]
    dphs: LonLatField[Float64]
#     #Tendencies from physics
#     dudt_phys : LonLatField[Float64]
#     dvdt_phys : LonLatField[Float64]
#     dtdt_phys : LonLatField[Float64]
#     dshdt_phys : LonLatField[Float64]

#     #Individual tendencies
    qhv : LonLatField[Float64]
    qhu : LonLatField[Float64]
    dkedlon : LonLatField[Float64]
    dkedlat : LonLatField[Float64]
    dmfdlon : LonLatField[Float64]
    dmfdlat : LonLatField[Float64]
    dptfdlon : LonLatField[Float64]
    dptfdlat : LonLatField[Float64]
    dptfdlev : LonLatField[Float64]
    pgf_lon : LonLatField[Float64]
    pgf_lat : LonLatField[Float64]
    wedudlev : LonLatField[Float64]
    wedvdlev : LonLatField[Float64]
    smag_dptdt : LonLatField[Float64]
    smag_dudt : LonLatField[Float64]
    smag_dvdt : LonLatField[Float64]

du = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
dv = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
dgz = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dpt = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dphs = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))

qhv = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
qhu = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
dkedlon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
dkedlat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
dmfdlon = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dmfdlat = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dptfdlon = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dptfdlat = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
dptfdlev = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
pgf_lon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
pgf_lat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))
wedudlev = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
wedvdlev = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))


smag_dptdt = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,gm.full_nlev),(True,True,True))
smag_dudt  = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,gm.full_nlev),(False,True,True))
smag_dvdt  = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,gm.full_nlev),(True,False,True))

tend = HybridTendField(du, dv, dgz, dpt, dphs, qhv,qhu,
                       dkedlon, dkedlat,
                       dmfdlon, dmfdlat,
                       dptfdlon, dptfdlat, dptfdlev,
                       pgf_lon, pgf_lat,
                       wedudlev, wedvdlev,
                       smag_dptdt, smag_dudt, smag_dvdt
                       )
# tend = HybridTendField(du, dv, dgz, dpt, dphs, qhv, qhu,
#                         dkedlon, dkedlat, dmfdlon, dmfdlat, 
#                         dptfdlon, dptfdlat, dptfdlev,
#                         pgf_lon, pgf_lat, wedudlev, wedvdlev,
#                         smag_dptdt, smag_dudt, smag_dvdt)


@hybrid
class HybridStaticField:
    # landmast : LonLatField[Float64]
    gzs : LonLatField[Float64]
    # zs_std : LonLatField[Float64]
    dzsdlon : LonLatField[Float64]
    dzsdlat : LonLatField[Float64]
    # ref_ps  : LonLatField[Float64]

# landmask = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))
gzs = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))
# zs_std = LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))
dzsdlon = LonLatField(Float64,(gm.half_nlon,gm.full_nlat,1),(False,True,True))
dzsdlat = LonLatField(Float64,(gm.full_nlon,gm.half_nlat,1),(True,False,True))
# ref_ps =  LonLatField(Float64,(gm.full_nlon,gm.full_nlat,1),(True,True,True))

staticv = HybridStaticField(gzs,dzsdlon,dzsdlat)
#static = HybridStaticField(landmask, gzs, zs_std, dzsdlon, dzsdlat, ref_ps)

@hybrid
class HybridTendPara:
    phs: Bool
    pt : Bool
    gz : Bool
    u  : Bool
    v  : Bool

tendPara = HybridTendPara(False, False, False, False, False)

