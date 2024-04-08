import math

from frontend.pytoy.lang import func, Float64, Int32

import frontend.mcv_stroud.const_value as cv
from frontend.mcv_stroud.utils import Vector2, Vector3, Vector4

@func 
def matrixA(k:Float64, lamb:Float64, theta:Float64)->Vector4:
    ret = Vector4(0.0,0.0,0.0,0.0)

    if (k <= 4):
        alambda = lamb-(k-1)*cv.pi/2.0 
        atheta  = theta
        a=math.sin(alambda)
        b=math.cos(alambda)
        c=math.sin(atheta) 
        d=math.cos(atheta)
        ret.a11 = d 
        ret.a12 = 0.0
        ret.a21 = -c*d*a/b 
        ret.a22 = b*d*d+c*c/b
    elif (k == 5):
        alambda=lamb 
        atheta=theta
        a= math.sin(alambda)
        b= math.cos(alambda)
        c= math.sin(atheta) 
        d= math.cos(atheta) 
        temp =1.0+a*a*d*d/c/c
        ret.a11=b*c*temp
        ret.a21=-c*c*a*temp
        temp=1.0+b*b*d*d/c/c
        ret.a12=a*c*temp 
        ret.a22=b*c*c*temp
    else:
        alambda=lamb
        atheta=theta
        a=math.sin(alambda) 
        b=math.cos(alambda) 
        c=math.sin(atheta) 
        d=math.cos(atheta) 
        temp=1.+a*a*d*d/c/c
        ret.a11=-b*c*temp
        ret.a21=c*c*a*temp 
        temp=1.+b*b*d*d/c/c
        ret.a12=a*c*temp 
        ret.a22=b*c*c*temp

    return ret

@func
def matrixIG(x:Float64, y:Float64)->Vector4:
    ret = Vector4(0.0,0.0,0.0,0.0)
    rho= math.sqrt(1.0+math.tan(x/cv.r)**2+math.tan(y/cv.r)**2)
    ret.a11 = (1.0+math.tan(y/cv.r)**2)*(rho**2*math.cos(x/cv.r)**2*math.cos(y/cv.r)**2)
    ret.a12 = math.tan(x/cv.r)*math.tan(y/cv.r)*(rho**2*math.cos(x/cv.r)**2*math.cos(y/cv.r)**2)
    ret.a21 = ret.a12
    ret.a22 = (1.0+math.tan(x/cv.r)**2)*(rho**2*math.cos(x/cv.r)**2*math.cos(y/cv.r)**2)

    return ret 

@func
def pprop2sp(k:Int32, x:Float64, y:Float64)->Vector2:
    x1 = x / cv.r
    y1 = y / cv.r
    ret = Vector2(0.0,0.0)
    if (k <= 4):
        ret.x = x1 + (k-1) * cv.pi / 2.0
        ret.y = math.atan2(math.tan(y1)*math.cos(x1),1.0)
    elif ( k == 5):
        a = math.tan(x1)
        b = math.tan(y1)
        ret.x = math.atan2(a, -b)
        ret.y = math.atan2(1.0, math.sqrt(a*a+b*b))
    else:
        a = math.tan(x1)
        b = math.tan(y1)
        ret.x = math.atan2(a, b)
        ret.y = -math.atan2(1.0, math.sqrt(a*a+b*b))
    return ret
        

@func
def covprosp2p(k:Int32, sv1:Float64, sv2:Float64, lamb:Float64, theta:Float64)->Vector2:
    a = Vector4(0.0,0.0,0.0,0.0)
    a = matrixA(k,lamb,theta)

    ret = Vector2(0.0,0.0)

    ret.x=a.a11*sv1+a.a21*sv2
    ret.y=a.a12*sv1+a.a22*sv2

    return ret

@func
def cov2contrav(cov1:Float64, cov2:Float64, x:Float64, y:Float64)->Vector2:
    a = Vector4(0.0,0.0,0.0,0.0)
    a = matrixIG(x,y)
    ret = Vector2(0.0,0.0)
    ret.x = a.a11*cov1 + a.a12*cov2
    ret.y = a.a21*cov1 + a.a22*cov2
    return ret


@func
def computeh(lamb: Float64, theta: Float64)->Float64:
    h = cv.gh0-(cv.r*cv.omega*cv.u0+0.5*cv.u0*cv.u0)*(-math.cos(lamb)*math.cos(theta)*math.sin(0)+math.sin(theta)*math.cos(0))**2
    h = h / cv.gra
    return h

@func
def computejab(x: Float64, y: Float64)->Float64:
    alpha = x /cv.r
    beta = y / cv.r
    jab = 1.0/(1.0+math.tan(alpha)*math.tan(alpha)+math.tan(beta)*math.tan(beta))**1.5/math.cos(alpha)/math.cos(alpha)/math.cos(beta)/math.cos(beta)
    return jab

@func
def computevel(k:Int32, lamb: Float64, theta: Float64)->Vector2:
    alpha=0.0
    us=cv.u0*(math.cos(alpha)*math.cos(theta)+math.sin(alpha)*math.cos(lamb)*math.sin(theta))
    vs=-cv.u0*math.sin(alpha)*math.sin(lamb)

    ret = covprosp2p(k,us,vs,lamb, theta)

    return ret

