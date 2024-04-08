from pytoy.lang import Float64, Int32, String, hybrid, space_op, extern_func, const_init  
from pytoy.lang.dtype_ext import LonLatField 
import math
import  gmcore_smallball.const_value as cv
import  gmcore_smallball.global_config as gc
import  gmcore_smallball.sphere_geometry as sg
from  gmcore_smallball.utils import Vector3, MeshVector3, MeshVector6, disect

full_nlon  = gc.nlon
half_nlon  = gc.nlon
full_ids = 0
full_ide = full_nlon - 1
half_ids = 0
half_ide = half_nlon - 1
full_nlat  = gc.nlat
half_nlat  = gc.nlat - 1
full_jds = 0
full_jde = full_nlat - 1
half_jds = 0
half_jde = half_nlat - 1
full_nlev  = gc.nlev
half_nlev  = full_nlev + 1
full_kds = 0
full_kde = full_nlev - 1
half_kds = 0
half_kde = half_nlev - 1

start_lon     = 0
end_lon       = cv.pi2
start_lat     = -cv.pi05
end_lat       = cv.pi05

full_jds_no_pole = full_jds + 1
full_jde_no_pole = full_jde - 1

lon_hw = 2
lat_hw = 2

full_ims = full_ids - lon_hw
full_ime = full_ide + lon_hw
full_jms = full_jds - lat_hw
full_jme = full_jde + lat_hw
half_ims = half_ids - lon_hw
half_ime = half_ide + lon_hw
half_jms = half_jds - lat_hw
half_jme = half_jde + lat_hw
full_kms = full_kds - 2
full_kme = full_kde + 2
half_kms = half_kds - 2
half_kme = half_kde + 2


@hybrid
class HybridMeshField:
    dlat: LonLatField[Float64]
    # full_dlev: LonLatField[Float64]
    # half_dlev: LonLatField[Float64]
    # half_dlev_upper: LonLatField[Float64]
    # half_dlev_lower: LonLatField[Float64]

    full_lon: LonLatField[Float64]
    half_lon: LonLatField[Float64]
    full_lat: LonLatField[Float64]
    half_lat: LonLatField[Float64]
    full_lev: LonLatField[Float64]
    half_lev: LonLatField[Float64]
    full_cos_lon: LonLatField[Float64]
    half_cos_lon: LonLatField[Float64]
    full_sin_lon: LonLatField[Float64]
    half_sin_lon: LonLatField[Float64]
    full_cos_lat: LonLatField[Float64]
    half_cos_lat: LonLatField[Float64]
    full_sin_lat: LonLatField[Float64]
    half_sin_lat: LonLatField[Float64]
    # for output
    full_lon_deg: LonLatField[Float64]
    half_lon_deg: LonLatField[Float64]
    full_lat_deg: LonLatField[Float64]
    half_lat_deg: LonLatField[Float64]
    # Area for weighting
    area_cell: LonLatField[Float64]
    area_lon: LonLatField[Float64]
    area_lon_west: LonLatField[Float64]
    area_lon_east: LonLatField[Float64]
    area_lon_north: LonLatField[Float64]
    area_lon_south: LonLatField[Float64]
    area_lat: LonLatField[Float64]
    area_lat_west: LonLatField[Float64]
    area_lat_east: LonLatField[Float64]
    area_lat_north: LonLatField[Float64]
    area_lat_south: LonLatField[Float64]
    area_vtx: LonLatField[Float64]
    area_subcell_0: LonLatField[Float64]
    area_subcell_1: LonLatField[Float64]
    # Edge length
    de_lon: LonLatField[Float64]
    de_lat: LonLatField[Float64]
    le_lat: LonLatField[Float64]
    le_lon: LonLatField[Float64]
    # Coriolis parameters
    full_f: LonLatField[Float64]
    half_f: LonLatField[Float64]
    # Weight for constructing tangential wind
    full_tangent_wgt_0: LonLatField[Float64]
    full_tangent_wgt_1: LonLatField[Float64]
    half_tangent_wgt_0: LonLatField[Float64] 
    half_tangent_wgt_1: LonLatField[Float64]

    #vert_coord
    hyai : LonLatField[Float64]
    hybi : LonLatField[Float64]
    hyam : LonLatField[Float64]
    hybm : LonLatField[Float64]

    #div_damp
    c_lon : LonLatField[Float64]
    c_lat : LonLatField[Float64]

# full_dlev = LonLatField(Float64, (1,1,full_nlev),  const = True )
# half_dlev = LonLatField(Float64, (1,1,full_nlev),  const = True )
full_lev = LonLatField(Float64, (1,1,full_nlev),  (True,True,True), const = True )
half_lev = LonLatField(Float64, (1,1,half_nlev),  (True,True,False), const = True )

dlat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )

full_lon = LonLatField(Float64, (full_nlon,1,1), (True,True,True), const = True )
half_lon = LonLatField(Float64, (half_nlon,1,1), (False,True,True), const = True )
full_lat = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
half_lat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )

full_cos_lon = LonLatField(Float64, (full_nlon,1,1), (True,True,True), const = True )
half_cos_lon = LonLatField(Float64, (half_nlon,1,1), (False,True,True), const = True )
full_sin_lon = LonLatField(Float64, (full_nlon,1,1), (True,True,True), const = True )
half_sin_lon = LonLatField(Float64, (half_nlon,1,1), (False,True,True), const = True )
full_cos_lat = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
half_cos_lat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
full_sin_lat = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
half_sin_lat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
full_lon_deg = LonLatField(Float64, (full_nlon,1,1),  (True,True,True), const = True )
half_lon_deg = LonLatField(Float64, (half_nlon,1,1), (False,True,True),  const = True )
full_lat_deg = LonLatField(Float64, (1,full_nlat,1),  (True,True,True), const = True )
half_lat_deg = LonLatField(Float64, (1,half_nlat,1),  (True,False,True), const = True )
area_cell    = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lon     = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lon_west = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lon_east = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lon_north = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lon_south = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_lat      = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_lat_west = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_lat_east = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_lat_north = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_lat_south = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_vtx = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
area_subcell_0 = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
area_subcell_1 = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
de_lon = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
de_lat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
le_lat = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
le_lon = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
full_f = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
half_f = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
full_tangent_wgt_0 = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
full_tangent_wgt_1 = LonLatField(Float64, (1,full_nlat,1), (True,True,True), const = True )
half_tangent_wgt_0 = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )
half_tangent_wgt_1 = LonLatField(Float64, (1,half_nlat,1), (True,False,True), const = True )

hyai = LonLatField(Float64, (1,1,half_nlev),  (True,True,False), const = True )
hybi = LonLatField(Float64, (1,1,half_nlev),  (True,True,False), const = True )
hyam = LonLatField(Float64, (1,1,full_nlev),  (True,True,True), const = True )
hybm = LonLatField(Float64, (1,1,full_nlev),  (True,True,True), const = True )

c_lon = LonLatField(Float64, (1,full_nlat,full_nlev),  (True,True,True), const = True )
c_lat = LonLatField(Float64, (1,half_nlat,full_nlev),  (True,False,True), const = True )

mesh = HybridMeshField(dlat,full_lon,half_lon,full_lat,half_lat, full_lev, half_lev,
                       full_cos_lon,half_cos_lon,full_sin_lon,half_sin_lon,
                       full_cos_lat,half_cos_lat,full_sin_lat,half_sin_lat,
                       full_lon_deg, half_lon_deg, full_lat_deg, half_lat_deg,
                       area_cell,
                       area_lon,area_lon_west,area_lon_east,area_lon_north,area_lon_south,
                       area_lat,area_lat_west,area_lat_east,area_lat_north,area_lat_south,
                       area_vtx,area_subcell_0,area_subcell_1,
                       de_lon,de_lat,le_lat,le_lon,
                       full_f,half_f,
                       full_tangent_wgt_0,full_tangent_wgt_1,half_tangent_wgt_0,half_tangent_wgt_1,
                       hyai, hybi, hyam, hybm,
                       c_lon, c_lat)
                       

# mesh = HybridMeshField(dlat,full_lon,half_lon,full_lat,half_lat, full_lev, half_lev,
#                     full_cos_lon,half_cos_lon,full_sin_lon,half_sin_lon,
#                     full_cos_lat,half_cos_lat,full_sin_lat,half_sin_lat,
#                     full_lon_deg,half_lon_deg,full_lat_deg,half_lat_deg,
                    # area_cell,
                    # area_lon,area_lon_west,area_lon_east,area_lon_north,area_lon_south,
                    # area_lat,area_lat_west,area_lat_east,area_lat_north,area_lat_south,
                    # area_vtx,area_subcell,
#                     de_lon,de_lat,le_lat,le_lon,
#                     full_f,half_f,full_tangent_wgt,half_tangent_wgt)

@space_op 
@const_init
def LatLonMeshInit(mesh:HybridMeshField):

    t = 0
    loopinfo = MeshVector6(1, 1, full_nlon,0,0,-1)
    dlon = Float64((end_lon - start_lon) / full_nlon)
    for i in range(-1,full_nlon):
        v_inx = disect(t,loopinfo)
        t = t+1
        mesh.full_lon[i,0,0] = start_lon + (v_inx.z - 1) * dlon
        mesh.half_lon[i,0,0] = mesh.full_lon[i,0,0] + 0.5 * dlon
        mesh.full_lon_deg[i,0,0] = mesh.full_lon[i,0,0] * cv.deg
        mesh.half_lon_deg[i,0,0] = mesh.half_lon[i,0,0] * cv.deg   

    # #Set initial guess latitudes of full merdional grids.
    dlat0 = Float64((end_lat - start_lat) / half_nlat)
    t = 0
    loopinfo = MeshVector6(1, half_nlat, 1,0,1,0)
    for j in range(0,half_nlat):
        v_inx = disect(t,loopinfo)
        t = t+1
        mesh.half_lat[0,j,0] = start_lat + (v_inx.y - 0.5) * dlat0
        if math.fabs(mesh.half_lat[0,j,0]) < 1.0e-12:
            mesh.half_lat[0,j,0] = 0.0

    if gc.coarse_pole_mul != 0:
        #Calculate real dlat which is large at polar region.
        dlat0 = dlon
        for j in range(0,half_nlat):
            mesh.dlat[0,j,0] = dlat0 * (1 + (gc.coarse_pole_mul - 1) * math.exp(-gc.coarse_pole_decay * (abs(mesh.half_lat[0,j,0]) - cv.pi05)**2))
        tmpsum = Float64(0.0)
        for j in range(0,half_nlat):
            tmpsum = tmpsum + mesh.dlat[0,j,0]
        for j in range(0,half_nlat):
            mesh.dlat[0,j,0] = mesh.dlat[0,j,0] * cv.pi / tmpsum
    else:
        for j in range(0,half_nlat):
            mesh.dlat[0,j,0] = dlat0

    #Set latitudes of full merdional grids.
    for j in range(0,1):
        mesh.full_lat[0,j,0] = start_lat
        mesh.full_lat_deg[0,j,0] = start_lat * cv.deg
    for j in range(1,full_nlat-1):
        mesh.full_lat[0,j,0] = mesh.full_lat[0,j-1,0] + mesh.dlat[0,j-1,0]
        if math.fabs(mesh.full_lat[0,j,0]) < 1.0e-12:
            mesh.full_lat[0,j,0] = 0.0
        mesh.full_lat_deg[0,j,0] = mesh.full_lat[0,j,0] * cv.deg
    for j in range(full_nlat-1,full_nlat):
        mesh.full_lat[0,j,0] = end_lat
        mesh.full_lat_deg[0,j,0] = end_lat * cv.deg

    #Set latitudes of half merdional grids.
    for j in range(0,half_nlat):
        if not (mesh.full_lat[0,j,0] == cv.pi05):
            mesh.half_lat[0,j,0] = mesh.full_lat[0,j,0] + 0.5 * mesh.dlat[0,j,0]
            if math.fabs(mesh.half_lat[0,j,0]) < 1e-12:
                mesh.half_lat[0,j,0] = 0.0
            mesh.half_lat_deg[0,j,0] = mesh.half_lat[0,j,0] * cv.deg

    #ToCheck 处理带halo的对称异常麻烦,该段对称赋值只是为了消除极小的误差
    #Ensure the grids are equatorial symmetry.
    # for j in range(0,full_nlat):
    #     if mesh.full_lat[0,j,0] > 0:
    #         mesh.full_lat[0,j,0] = -mesh.full_lat[0,full_nlat - j + 1,0]
    #         mesh.full_lat_deg[0,j,0] = -mesh.full_lat_deg[0,full_nlat - j + 1,0]
    # for j in range(0,half_nlat):
    #     if mesh.half_lat[0,j,0] > 0:
    #         mesh.half_lat[0,j,0] = -mesh.half_lat[0,half_nlat - j + 1,0]
    #         mesh.half_lat_deg[0,j,0] = -mesh.half_lat_deg[0,half_nlat - j + 1,0]

    for i in range(-1,full_nlon):
        mesh.full_cos_lon[i,0,0] = math.cos(mesh.full_lon[i,0,0])
        mesh.full_sin_lon[i,0,0] = math.sin(mesh.full_lon[i,0,0])

    for i in range(-1,half_nlon+1):
        mesh.half_cos_lon[i,0,0] = math.cos(mesh.half_lon[i,0,0])
        mesh.half_sin_lon[i,0,0] = math.sin(mesh.half_lon[i,0,0])

    for j in range(-1,half_nlat):
        if mesh.half_lat[0,j,0] >= -cv.pi05 and mesh.half_lat[0,j,0] <= cv.pi05:
            mesh.half_cos_lat[0,j,0] = math.cos(mesh.half_lat[0,j,0])
            mesh.half_sin_lat[0,j,0] = math.sin(mesh.half_lat[0,j,0])

    for j in range(-1,full_nlat):
        if mesh.full_lat[0,j,0] >= -cv.pi05 and mesh.full_lat[0,j,0] <= cv.pi05:
            mesh.full_cos_lat[0,j,0] = math.cos(mesh.full_lat[0,j,0])
            mesh.full_sin_lat[0,j,0] = math.sin(mesh.full_lat[0,j,0])

    # #Ensure the values of cos_lat and sin_lat are expected at the Poles.
    
    for j in range(0,1):
        mesh.full_cos_lat[0,j,0] =  0.0
        mesh.full_sin_lat[0,j,0] = -1.0
    for j in range(full_jde,full_jde+1):
        mesh.full_cos_lat[0,j,0] =  0.0
        mesh.full_sin_lat[0,j,0] =  1.0

    for j in range(0,1):
        mesh.area_cell[0,j,0] = cv.radius**2 * dlon * (mesh.half_sin_lat[0,j,0] +1.0)
        mesh.area_subcell_1[0,j,0] = cv.radius**2 * 0.5 * dlon * (mesh.half_sin_lat[0,j,0] +1.0)

    for j in range(full_jde,full_jde+1):
        mesh.area_cell[0,j,0] = cv.radius**2 * dlon * (1.0 - mesh.half_sin_lat[0,j-1,0])
        mesh.area_subcell_0[0,j,0] = cv.radius**2 * 0.5 * dlon * (1.0 - mesh.half_sin_lat[0,j-1,0])

    v1 = Vector3(0.0, 0.0, 0.0)
    v2 = Vector3(0.0, 0.0, 0.0)
    v3 = Vector3(0.0, 0.0, 0.0)
    buf = Vector3(0.0, 0.0, 0.0)
    # x = Array((3,),Float64)
    # y = Array((3,),Float64)
    # z = Array((3,),Float64)

    for j in range(1,full_jde - 1 + 1):
        for i in range(0,1):
            mesh.area_cell[0,j,0] = cv.radius**2 * dlon * (mesh.half_sin_lat[0,j,0] - mesh.half_sin_lat[0,j-1,0])
            mesh.area_subcell_0[0,j,0] = cv.radius**2 * 0.5 * dlon * (mesh.full_sin_lat[0,j,0] - mesh.half_sin_lat[0,j-1,0])
            mesh.area_subcell_1[0,j,0] = cv.radius**2 * 0.5 * dlon * (mesh.half_sin_lat[0,j,0] - mesh.full_sin_lat[0,j,0])

            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.full_lon[0], mesh.full_lat[0,j,0])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[j-1])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[0,j,0])
            buf = sg.cartesian_transform(mesh.full_lon[i,0,0], mesh.full_lat[0,j,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j-1,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lon_west[0,j,0] = sg.calc_area(v1,v2,v3)
            mesh.area_lon_east[0,j,0] = mesh.area_lon_west[0,j,0]
            mesh.area_lon[0,j,0] = mesh.area_lon_west[0,j,0] + mesh.area_lon_east[0,j,0]

            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[0,j,0])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.full_lon[0], mesh.full_lat[0,j,0])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[0,j,0])
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i,0,0], mesh.full_lat[0,j,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lon_north[0,j,0] = sg.calc_area_with_last_small_arc(v1,v2,v3)

            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[j-1])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[0,j,0])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.full_lon[0], mesh.full_lat[0,j,0])
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j-1,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i,0,0], mesh.full_lat[0,j,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lon_south[0,j,0] = sg.calc_area_with_last_small_arc(v1,v2,v3)

    for j in range(half_jds,half_jde+1):
        for i in range(0,1):
            mesh.area_vtx[0,j,0] = cv.radius**2 * dlon * (mesh.full_sin_lat[0,j+1,0] - mesh.full_sin_lat[0,j,0])
            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[0,j,0])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[0,j,0])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[j+1])
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j+1,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lat_west[0,j,0] = sg.calc_area(v1,v2,v3)
            mesh.area_lat_east[0,j,0] = mesh.area_lat_west[0,j,0]

    #Reset up or down area to polar sector area.
    for j in range(half_jds,half_jde):
        for i in range(0,1):
            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[j+1])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[0,j,0])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.half_lon[1], mesh.half_lat[0,j,0])
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j+1,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i+1,0,0], mesh.half_lat[0,j,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lat_north[0,j,0] = sg.calc_area_with_last_small_arc(v1,v2,v3)
    for j in range(half_jde,half_jde+1):
        mesh.area_lat_north[0,j,0] = mesh.area_cell[0,j+1,0]

    for j in range(half_jds,half_jds+1):
        mesh.area_lat_south[0,j,0] = mesh.area_cell[0,j,0]
    for j in range(half_jds+1,half_jde+1):
        for i in range(0,1):
            # x[0],y[0],z[0] = sg.cartesian_transform(mesh.full_lon[1], mesh.full_lat[0,j,0])
            # x[1],y[1],z[1] = sg.cartesian_transform(mesh.half_lon[1], mesh.half_lat[0,j,0])
            # x[2],y[2],z[2] = sg.cartesian_transform(mesh.half_lon[0], mesh.half_lat[0,j,0])
            buf = sg.cartesian_transform(mesh.full_lon[i+1,0,0], mesh.full_lat[0,j,0])
            v1.x, v2.x, v3.x = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i+1,0,0], mesh.half_lat[0,j,0])
            v1.y, v2.y, v3.y = buf.x, buf.y, buf.z
            buf = sg.cartesian_transform(mesh.half_lon[i,0,0], mesh.half_lat[0,j,0])
            v1.z, v2.z, v3.z = buf.x, buf.y, buf.z
            mesh.area_lat_south[0,j,0] = sg.calc_area_with_last_small_arc(v1,v2,v3)
        
    for j in range(half_jds,half_jde+1):
        mesh.area_lat[0,j,0] = mesh.area_lat_north[0,j,0] + mesh.area_lat_south[0,j,0]

    for j in range(full_jds_no_pole,full_jde_no_pole+1):
        mesh.de_lon[0,j,0] = cv.radius * mesh.full_cos_lat[0,j,0] * dlon
        mesh.le_lon[0,j,0] = 2.0 * mesh.area_lon[0,j,0] / mesh.de_lon[0,j,0]
    for j in range(full_jds,full_jds+1):
        mesh.le_lon[0,j,0] = 0.0
        mesh.de_lon[0,j,0] = 0.0
    for j in range(full_jde,full_jde+1):
        mesh.le_lon[0,j,0] = 0.0
        mesh.de_lon[0,j,0] = 0.0

    for j in range(half_jds,half_jde+1):
        mesh.le_lat[0,j,0] = cv.radius * mesh.half_cos_lat[0,j,0] * dlon
        mesh.de_lat[0,j,0] = 2.0 * mesh.area_lat[0,j,0] / mesh.le_lat[0,j,0]

    # # TODO: use integer Enum to replace it
    if gc.tangent_wgt_scheme == 1:
        for j in range(full_jds_no_pole,full_jde_no_pole+1):
            mesh.full_tangent_wgt_0[0,j,0] = mesh.le_lat[0,j-1,0] / mesh.de_lon[0,j,0] * 0.25
            mesh.full_tangent_wgt_1[0,j,0] = mesh.le_lat[0,j  ,0] / mesh.de_lon[0,j,0] * 0.25
        for j in range(0,half_nlat):
            mesh.half_tangent_wgt_0[0,j,0] = mesh.le_lon[0,j,0]   / mesh.de_lat[0,j,0] * 0.25
            mesh.half_tangent_wgt_1[0,j,0] = mesh.le_lon[0,j+1,0] / mesh.de_lat[0,j,0] * 0.25

        for j in range(0,full_nlat):
            mesh.full_f[0,j,0] = 2.0 * cv.omega * mesh.full_sin_lat[0,j,0]
        for j in range(0,half_nlat):
            mesh.half_f[0,j,0] = 2.0 * cv.omega * mesh.half_sin_lat[0,j,0]

    #hybrid_coord_init
    for k in range(0,1):
        mesh.hyai[0,0,k] = 0.0000000000000000
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(1,2):
        mesh.hyai[0,0,k] = 0.0027008056640625
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(2,3):
        mesh.hyai[0,0,k] = 0.0076868347823620
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(3,4):
        mesh.hyai[0,0,k] = 0.0158526937798342
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(4,5):
        mesh.hyai[0,0,k] = 0.0276367470939261
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(5,6):
        mesh.hyai[0,0,k] = 0.0414515552099706
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(6,7):
        mesh.hyai[0,0,k] = 0.0562277844155789
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(7,8):
        mesh.hyai[0,0,k] = 0.0724737972365225
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(8,9):
        mesh.hyai[0,0,k] = 0.0893795253309864
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(9,10):
        mesh.hyai[0,0,k] = 0.1063548465048298
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(10,11):
        mesh.hyai[0,0,k] = 0.1254725257849640
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(11,12):
        mesh.hyai[0,0,k] = 0.1479642011314251
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(12,13):
        mesh.hyai[0,0,k] = 0.1744233770857367
        mesh.hybi[0,0,k] = 0.0000000000000000
    for k in range(13,14):
        mesh.hyai[0,0,k] = 0.2054476619883455
        mesh.hybi[0,0,k] = 0.0001056069805219
    for k in range(14,15):
        mesh.hyai[0,0,k] = 0.2362396334209435
        mesh.hybi[0,0,k] = 0.0059397906834944
    for k in range(15,16):
        mesh.hyai[0,0,k] = 0.2617261684535087
        mesh.hybi[0,0,k] = 0.0235374270991938
    for k in range(16,17):
        mesh.hyai[0,0,k] = 0.2783057827853432
        mesh.hybi[0,0,k] = 0.0576439176162314
    for k in range(17,18):
        mesh.hyai[0,0,k] = 0.2816243248017372
        mesh.hybi[0,0,k] = 0.1139578493851603
    for k in range(18,19): 
        mesh.hyai[0,0,k] = 0.2683280168286641
        mesh.hybi[0,0,k] = 0.1932313936022966
    for k in range(19,20):
        mesh.hyai[0,0,k] = 0.2429473354807992
        mesh.hybi[0,0,k] = 0.2808764820852269
    for k in range(20,21):
        mesh.hyai[0,0,k] = 0.2108099294331058
        mesh.hybi[0,0,k] = 0.3715629284436812
    for k in range(21,22):
        mesh.hyai[0,0,k] = 0.1759998110016148
        mesh.hybi[0,0,k] = 0.4612091687348547
    for k in range(22,23): 
        mesh.hyai[0,0,k] = 0.1415334037742616
        mesh.hybi[0,0,k] = 0.5467954414732804
    for k in range(23,24):
        mesh.hyai[0,0,k] = 0.1095123294032484
        mesh.hybi[0,0,k] = 0.6262219272258057
    for k in range(24,25):
        mesh.hyai[0,0,k] = 0.0812761592574764
        mesh.hybi[0,0,k] = 0.6981497395977654 
    for k in range(25,26):
        mesh.hyai[0,0,k] = 0.0575423961814687
        mesh.hybi[0,0,k] = 0.7618575996945027
    for k in range(26,27):
        mesh.hyai[0,0,k] = 0.0385285881906738
        mesh.hybi[0,0,k] = 0.8171268253424389
    for k in range(27,28):
        mesh.hyai[0,0,k] = 0.0240694984015316
        mesh.hybi[0,0,k] = 0.8641240158900223
    for k in range(28,29):
        mesh.hyai[0,0,k] = 0.0137226170696702
        mesh.hybi[0,0,k] = 0.9032994538662682
    for k in range(29,30):
        mesh.hyai[0,0,k] = 0.0068720184783044
        mesh.hybi[0,0,k] = 0.9352560448161740
    for k in range(30,31):
        mesh.hyai[0,0,k] = 0.0027973066609278
        mesh.hybi[0,0,k] = 0.9607227915792738
    for k in range(31,32):
        mesh.hyai[0,0,k] = 0.0007580797514921
        mesh.hybi[0,0,k] = 0.9804356127137464
    for k in range(32,33):
        mesh.hyai[0,0,k] = 0.0000000000000000
        mesh.hybi[0,0,k] = 1.0000000000000000
    
    for k in range(0,full_nlev):
        mesh.hyam[0,0,k] = 0.5 * (mesh.hyai[0,0,k] + mesh.hyai[0,0,k+1])
        mesh.hybm[0,0,k] = 0.5 * (mesh.hybi[0,0,k] + mesh.hybi[0,0,k+1])
    
    for k in range(0,full_nlev):
        mesh.full_lev[0,0,k] = mesh.hyam[0,0,k] + mesh.hybm[0,0,k]
    for k in range(0,half_nlev):
        mesh.half_lev[0,0,k] = mesh.hyai[0,0,k] + mesh.hybi[0,0,k]

    #div_init
    for k in range(full_kds, full_kde+1):
        for j in range(full_jds_no_pole, full_jde_no_pole+1):
            mesh.c_lon[0,j,k] = gc.div_damp_coef2 * max(1.0, 8 * (1 + math.tanh(math.log(gc.ptop / sg.hybrid_coord_calc_ph(mesh.hyam[0,0,k], mesh.hybm[0,0,k], gc.p0))))) * mesh.le_lon[0,j,0] * mesh.de_lon[0,j,0]
    
    for k in range(full_kds, full_kde+1):
        for j in range(half_jds, half_jde+1):
            mesh.c_lat[0,j,k] = gc.div_damp_coef2 * max(1.0, 8 * (1 + math.tanh(math.log(gc.ptop / sg.hybrid_coord_calc_ph(mesh.hyam[0,0,k], mesh.hybm[0,0,k], gc.p0))))) * mesh.le_lat[0,j,0] * mesh.de_lat[0,j,0]

    





