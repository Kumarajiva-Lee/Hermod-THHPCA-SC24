import math

from frontend.pytoy.lang import func, Float64

import gmcore_smallball.const_value as cv
import gmcore_smallball.global_config as gc
import gmcore_smallball.const_value as cv

from gmcore_smallball.utils import Vector2, Vector3, Vector4

@func
def specific_humidity(qv:Float64)->Float64:
    return qv / (1 + qv)

@func
def mixing_ratio(sh:Float64)->Float64:
    return sh / (1 - sh)

@func
def sign(x:Float64, y:Float64)->Float64:
    if y > 0.0:
        return x
    else:
        return -x

@func
def potential_temperature(t:Float64,p:Float64,qv:Float64)->Float64:
    return t * (gc.p0 / p)**cv.Rd_o_cpd * (1 + cv.Rv_o_Rd * qv)

@func
def temperature(pt:Float64,p:Float64,qv:Float64)->Float64:
    return pt * (p / gc.p0)**cv.Rd_o_cpd / (1 + cv.Rv_o_Rd * qv)

@func
def virtual_temperature(t:Float64,sh:Float64)->Float64:
    return t * (1 + (cv.Rv_o_Rd - 1) * sh)

@func
def dry_air_density(pt:Float64,p:Float64)->Float64:
    return gc.p0 / cv.Rd / pt * (p / gc.p0)**cv.cv_o_cpd

@func 
def moist_air_density(t:Float64,p:Float64,qv:Float64)->Float64:
    return p / cv.Rd / virtual_temperature(t, specific_humidity(qv))

@func
def upwind1(dir:Float64,wgt:Float64,f1:Float64,f2:Float64)->Float64:
    # return 0.5 * (f[1] + f[0]) -0.5  * (f[2] - f[1]) * wgt * dir
    return 0.5 * (f2 + f1) - 0.5 * (f2 - f1) * wgt * dir

@func
def upwind3(dir:Float64,wgt:Float64,f1:Float64,f2:Float64,f3:Float64,f4:Float64)->Float64:
    c31 = 7.0 / 12.0
    c32 = -1.0 / 12.0
    c33 = 1.0 / 12.0
    # return c31 * (f[2] + f[1]) + c32 * (f[3] + f[0]) + c33 * (f[3] - f[0] - 3 * (f[2] - f[1])) * wgt * dir
    return c31 * (f3 + f2) + c32 * (f4 + f1) + c33 * (f4 - f1 - 3 * (f3 - f2)) * wgt * dir

@func
def slope(fm1: Float64, f: Float64, fp1: Float64)->Float64:
    df = (fp1 - fm1) * 0.5
    df_min = 2 * (f - min(min(fm1, f),fp1))
    df_max = 2 * (max(max(fm1, f), fp1) - f)
    return sign(min(math.fabs(df),min(df_min,df_max)), df)

@func
def ppm(fm2: Float64, fm1: Float64, f: Float64, fp1: Float64, fp2: Float64)->Vector3:
    dfl = slope(fm2, fm1, f  )
    df  = slope(fm1, f  , fp1)
    dfr = slope(f  , fp1, fp2)
    
    fl = 0.5 * (fm1 + f) + (dfl - df) / 6.0
    fr = 0.5 * (fp1 + f) + (df - dfr) / 6.0

    fl = f - sign(min(math.fabs(df), math.fabs(fl - f)), df)
    fr = f + sign(min(math.fabs(df), math.fabs(fr - f)), df)
    f6 = 6 * f - 3 * (fl + fr)
    df = fr - fl

    return Vector3(fl,df,f6)

