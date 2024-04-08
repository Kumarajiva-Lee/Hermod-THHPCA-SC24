from pytoy.lang import hybrid, Float64, Int32, func

@hybrid
class Vector2:
    x: Float64
    y: Float64

@hybrid
class Vector3:
    x: Float64
    y: Float64
    z: Float64

@hybrid
class Vector4:
    e0: Float64
    e1: Float64
    e2: Float64
    e3: Float64

@hybrid
class MeshVector3:
    x: Int32
    y: Int32
    z: Int32

@hybrid
class MeshVector6:
    lx: Int32
    ly: Int32
    lz: Int32
    sx: Int32
    sy: Int32
    sz: Int32

@func
def disect(inx: Int32, loopinfo: MeshVector6) -> MeshVector3:
    k = inx % loopinfo.lz
    tmp = inx / loopinfo.lz
    j = tmp % loopinfo.ly
    tmp = tmp / loopinfo.ly
    i = tmp % loopinfo.lx
    return MeshVector3(i+loopinfo.sx,j+loopinfo.sy,k+loopinfo.sz)

@hybrid
class Vector32:
    e0: Float64
    e1: Float64
    e2: Float64
    e3: Float64
    e4: Float64
    e5: Float64
    e6: Float64
    e7: Float64
    e8: Float64
    e9: Float64
    e10: Float64
    e11: Float64
    e12: Float64
    e13: Float64
    e14: Float64
    e15: Float64
    e16: Float64
    e17: Float64
    e18: Float64
    e19: Float64
    e20: Float64
    e21: Float64
    e22: Float64
    e23: Float64
    e24: Float64
    e25: Float64
    e26: Float64
    e27: Float64
    e28: Float64
    e29: Float64
    e30: Float64
    e31: Float64
    e31: Float64

@hybrid
class Vector33:
    e0: Float64
    e1: Float64
    e2: Float64
    e3: Float64
    e4: Float64
    e5: Float64
    e6: Float64
    e7: Float64
    e8: Float64
    e9: Float64
    e10: Float64
    e11: Float64
    e12: Float64
    e13: Float64
    e14: Float64
    e15: Float64
    e16: Float64
    e17: Float64
    e18: Float64
    e19: Float64
    e20: Float64
    e21: Float64
    e22: Float64
    e23: Float64
    e24: Float64
    e25: Float64
    e26: Float64
    e27: Float64
    e28: Float64
    e29: Float64
    e30: Float64
    e31: Float64
    e32: Float64