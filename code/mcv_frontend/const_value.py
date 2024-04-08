import math

from frontend.pytoy.lang import space_op, Float64, Int32

import frontend.mcv_stroud.global_config as gc

r = 6371220.0
gra = 9.80616
omega = 7.292e-5
gh0 = 2.94e4
pi = 2.0* math.asin(1.0)
u0 = 2.0*pi*r/12.0/24.0/3600.0
dx = r * pi / 2.0 / float(gc.nc)
dxi = dx / 2
dy = r * pi / 2.0 / float(gc.nc)
dyi = dy / 2


