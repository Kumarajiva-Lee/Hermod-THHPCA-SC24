import miniweather_scuderia.const_value as cv
nx = 100 #400 #200 #100
nz = 50 #200 #100 #50
dx = cv.xlen / nx
dz = cv.zlen / nz
dt = min(dx,dz) / cv.max_speed * cv.cfl

#0 collision
#1 thermal
#2 gravity_waves
#3 density_current
#4 injection
initial = 1


