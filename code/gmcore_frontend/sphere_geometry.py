import math

from pytoy.lang import func, Float64
from gmcore_smallball.utils import Vector2, Vector3
import gmcore_smallball.global_config as gc

import gmcore_smallball.const_value as cv

@func
def cartesian_transform(lon:Float64,lat:Float64)->Vector3:
    cos_lat = math.cos(lat)
    x = cv.radius * cos_lat * math.cos(lon)
    y = cv.radius * cos_lat * math.sin(lon)
    z = cv.radius * math.sin(lat)

    # return x,y,z
    return Vector3(x,y,z)

@func
def inverse_cartesian_transform(x:Float64, y:Float64, z:Float64)->Vector2:
    lon = math.atan2(y,x)
    lat = math.asin(z / cv.radius)

    if lon < 0.0:
        lon += cv.pi2

    # return lon,lat
    return Vector2(lon,lat)

@func
def cross_product(a:Vector3, b:Vector3)->Vector3:
    # res = Array(dtype=Float64, length=3)
    # res[1] = x[2]*y[3]-x[3]*y[2]
    # res[2] = x[3]*y[1]-x[1]*y[3]
    # res[3] = x[1]*y[2]-x[2]*y[1]
    res = Vector3(0,0,0)
    res.x = a.y*b.z - a.z*b.y
    res.y = a.z*b.x - a.x*b.z
    res.z = a.x*b.y - a.y*b.x 
    return res

@func
def norm_vector(v:Vector3)->Vector3:
    res = Vector3(0,0,0)
    sum = v.x*v.x + v.y*v.y + v.z*v.z
    sum = math.sqrt(sum)
    if sum != 0:
        res.x = v.x / sum
        res.y = v.y / sum
        res.z = v.z / sum
    return res

@func
def calc_sphere_angle(a:Vector3, b:Vector3, c:Vector3)->Float64:
    # nab = Array((3,),Float64)
    # nbc = Array((3,),Float64)

    # nab = norm_vector(cross_product(a,b))
    # nbc = norm_vector(cross_product(b,c))

    # res = math.acos(- max(min(dot(nab,nbc),1.0),-1.0))

    # #Judge the cyclic direction with respect to point A to handle obtuse angle.
    # if dot(cross_product(nab,nbc),a) < 0.0:
    #     res = cv.pi2 - res
    nab = Vector3(0,0,0)
    nbc = Vector3(0,0,0)
    ta = Vector3(0,0,0)
    tb = Vector3(0,0,0)

    ta = cross_product(a,b)
    tb = cross_product(b,c)
    nab = norm_vector(ta)
    nbc = norm_vector(tb)

    dot = nab.x*nbc.x + nab.y*nbc.y + nab.z*nbc.z
    res = math.acos(- max(min(dot,1.0),-1.0))
    return res
    
@func
def calc_area(x:Vector3, y:Vector3, z:Vector3)->Float64:
    # n = x.length()
    # res = 0.0
    # angel = 0.0
    # vim = Array((3,),Float64)
    # vi  = Array((3,),Float64)
    # vip = Array((3,),Float64)

    # for i in range(0,n):
    #     if i != 0:
    #         im1 = i-1
    #     else:
    #         im1 = n
    #     if i != n:
    #         ip1 = i+1
    #     else:
    #         ip1 = 1
    #     vim = [x[im1],y[im1],z[im1]]
    #     vi  = [x[i],y[i],z[i]]
    #     vip = [x[ip1],y[ip1],z[ip1]]
    #     angle = calc_sphere_angle(vim,vi,vip)
    #     res += angel
    # res = cv.radius**2 * (res - (n - 2) * cv.pi)
    # return res
    n = 3
    res = Float64(0.0)
    angle = Float64(0.0)

    #1
    vi  = Vector3(x.x, y.x, z.x)
    #2
    vip = Vector3(x.y, y.y, z.y)
    #3
    vim = Vector3(x.z, y.z, z.z)

    angle = calc_sphere_angle(vim, vi, vip)
    res += angle

    angle = calc_sphere_angle(vi, vip, vim)
    res += angle

    angle = calc_sphere_angle(vip, vim, vi)
    res += angle

    res = cv.radius**2 * (res - (n - 2) * cv.pi)
    return res

@func 
def calc_area_with_last_small_arc(x:Vector3,y:Vector3,z:Vector3)->Float64:
    xv = Vector3(0,0,0)
    yv = Vector3(0,0,0)
    zv = Vector3(0,0,0)

    res = calc_area(x,y,z)
    # lon0,lat0 = inverse_cartesian_transform(x.x, y.x, z.x)
    # lon1,lat1 = inverse_cartesian_transform(x.y, y.y, z.y)
    # lon2,lat2 = inverse_cartesian_transform(x.z, y.z, z.z)
    lonlat0 = inverse_cartesian_transform(x.x, y.x, z.x)
    lat0 = lonlat0.y

    lonlat1 = inverse_cartesian_transform(x.y, y.y, z.y)
    lon1 = lonlat1.x
    lat1 = lonlat1.y

    lonlat2 = inverse_cartesian_transform(x.z, y.z, z.z)
    lon2 = lonlat2.x

    if lat1 == 0:
        return res
    else:
        if lat0 > lat1:
            dlon = lon2 - lon1
        else:
            dlon = lon1 - lon2
        
        if dlon < 0.0:
            dlon = dlon + cv.pi2

        xv.x = 0.0
        yv.x = 0.0

        if lat0*lat1 >=0 and math.fabs(lat0)>math.fabs(lat1):
            #Point 0 is at the side with the Pole.
            xv.y = x.y
            yv.y = y.y
            zv.y = z.y
            xv.z = x.z
            yv.z = y.z
            zv.z = z.z
        else:
            #Point 0 is at the opposite hemisphere.
            xv.y = x.z
            yv.y = y.z
            zv.y = z.z
            xv.z = x.y
            yv.z = y.y
            zv.z = z.y
        
        if lat1 > 0.0:
            #Small arc is at the North Sphere.
            zv.x = cv.radius
            area1 = cv.radius**2 * dlon * (1.0 - math.sin(lat1))
        else:
            #Small arc is at the South Sphere.
            zv.x = -cv.radius
            area1 = cv.radius**2 * dlon * (math.sin(lat1) + 1.0)
        
        area2 = calc_area(xv,yv,zv)
        area3 = area1 - area2

        if area3 < 0.0 and math.fabs(area3) > 1e-10:
            area3 = 0
        
        if lat0 * lat1 >=0 and math.fabs(lat0) > math.fabs(lat1):
            res += area3
        else:
            res -= area3
        
        return res
    
@func
def hybrid_coord_calc_ph(hyam: Float64, hybm: Float64, phs: Float64) -> Float64:
    return hyam * (gc.p0 - gc.ptop) + hybm * (phs - gc.ptop) + gc.ptop

@func
def hybrid_coord_calc_ph_lev(hyai: Float64, hybi: Float64, phs: Float64) -> Float64:
    return hyai * (gc.p0 - gc.ptop) + hybi * (phs - gc.ptop) + gc.ptop

@func
def hybrid_coord_calc_dphdt_lev(hybi : Float64, dphsdt : Float64) -> Float64:
    return hybi * dphsdt