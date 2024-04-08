import math

from frontend.pytoy.lang import func, Float64

import miniweather_scuderia.const_value as cv
from miniweather_scuderia.utils import Vector2, InitVector

@func
def sample_ellipse_cosine(x: Float64, z: Float64, amp: Float64, x0: Float64, z0: Float64, xrad: Float64, zrad: Float64)->Float64:
    dist = math.sqrt( ((x-x0)/xrad)*((x-x0)/xrad) + ((z-z0)/zrad)*((z-z0)/zrad) ) * cv.pi / 2.0
    if (dist <= cv.pi / 2.0):
        return amp * (math.cos(dist)**2)
    else:
        return 0.0


@func
def hydro_const_theta(z: Float64)->Vector2:
    theta0 = 300.0
    exner0 = 1.0
    t = theta0
    exner = exner0 - cv.grav * z / (cv.cp * theta0)
    p = cv.p0 * (exner**(cv.cp/cv.rd))
    rt = (p / cv.c0)**(1. / cv.gamm)
    r = rt / t
    return Vector2(r,t)

@func
def collision(x: Float64, z:Float64)->InitVector:
    buf = Vector2(0.0, 0.0)
    buf = hydro_const_theta(z)
    ret = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ret.r = 0.0
    ret.t = 0.0
    ret.u = 0.0
    ret.w = 0.0
    ret.t = ret.t + sample_ellipse_cosine(x,z, 20.,cv.xlen/2,2000.0,2000.0,2000.0)
    ret.t = ret.t + sample_ellipse_cosine(x,z,-20.,cv.xlen/2,8000.0,2000.0,2000.0)
    ret.hr = buf.x
    ret.ht = buf.y
    return ret

@func
def thermal(x: Float64, z:Float64)->InitVector:
    buf = Vector2(0.0, 0.0)
    buf = hydro_const_theta(z)
    ret = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ret.r = 0.0
    ret.u = 0.0
    ret.w = 0.0
    ret.t = 0.0 + sample_ellipse_cosine(x, z, 3.0, cv.xlen/2, 2000.0, 2000.0, 2000.0)
    ret.hr = buf.x
    ret.ht = buf.y
    return ret

@func
def hydro_const_bvfreq(z: Float64, bv_freq0: Float64)->Vector2:
    theta0 = 300.0
    exner0 = 1.0
    t = theta0 * math.exp( bv_freq0*bv_freq0 / cv.grav * z )
    exner = exner0 - cv.grav*cv.grav / (cv.cp * bv_freq0*bv_freq0) * (t - theta0) / (t * theta0)
    p = cv.p0 * (exner**(cv.cp/cv.rd))
    rt = (p / cv.c0)**(1.0 / cv.gamm)
    r = rt / t
    return Vector2(r,rt)


# x and z are input coordinates at which to sample
# r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
# hr and ht are output background hydrostatic density and potential temperature at that location
@func
def gravity_waves(x: Float64, z: Float64)->InitVector:
    buf = Vector2(0.0, 0.0)
    bv_freq0 = Float64(0.02)
    buf = hydro_const_bvfreq(z,bv_freq0)
    ret = InitVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ret.r = 0.0
    ret.t = 0.0
    ret.u = 15.0
    ret.w = 0.0
    ret.hr = buf.x
    ret.ht = buf.y
    return ret


