import math

from numpy import Inf
import  gmcore_smallball.global_config as gb

pi      = math.atan(1.0) * 4.0
pi2     = pi * 2.0
pi05    = pi * 0.5
deg     = 180.0 / pi
rad     = pi / 180.0
eps     = 1e-24 
inf     = Inf  
karman  = 0.4
if (gb.planet == 'earth'):
    omega      = 2.0 * pi / 86400.0
    radius     = 6371220.0
    g          = 9.80616
    Rd         = 287.04
    Rv         = 461.497
    cpd        = 1004.0
    cvd        = 717.0
    lapse_rate = 0.006
    Rd_o_Rv    = Rd / Rv
    Rv_o_Rd    = Rv / Rd
    p0         = 1.0


Rd_o_g   = Rd / g
Rd_o_cpd = Rd / cpd
cp_o_cvd = cpd / cvd
cv_o_cpd = cvd / cpd


