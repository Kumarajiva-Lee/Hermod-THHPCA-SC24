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
    r: Float64
    u: Float64
    w: Float64
    t: Float64

@hybrid
class InitVector:
    r: Float64
    u: Float64
    w: Float64
    t: Float64
    hr: Float64
    ht: Float64

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