CC=mpicc
SCC = mpicc
OPT = -O3
CM = -mhost 
CS = -mslave
MOPT = -lm -mftz -mfma -faddress_align=64
SOPT = -msimd -mfma -mftz -lm_slave -fmath-opt -faddress_align=64
CFLAG = -g -DGETTIMEOFDAY $(OPT)  -I/usr/sw/yyzlib/netcdf/include -I/usr/sw/yyzlib/hdf5-1.8.9/include
LDFLAGS = -g -L/usr/sw/yyzlib/netcdf/lib -L/usr/sw/yyzlib/hdf5-1.8.9/lib 
LDLIBS = -lnetcdf -lnetcdff  -lhdf5 -lhdf5_hl -lz -lm -lm_slave 
LINKOPT = -mhybrid


##
PWD_DIR=$(shell pwd)
LIB_DIR=$(PWD_DIR)/LIB
SRC_DIR=$(PWD_DIR)/src
OBJ_DIR=$(PWD_DIR)/obj

#TEST_DIR=$(PWD_DIR)/test

# INC_DIR=$(PWD_DIR)/inc

targets = app



##
export PWD_DIR OBJ_DIR LIB_DIR SRC_DIR CC SCC MOPT CM CS OPT SOPT CFLAG LDFLAGS LDLIBS LINKOPT targets 

##
all:
	make -C $(LIB_DIR)
	make -C $(SRC_DIR)
	make -C $(OBJ_DIR)

##
clean:
	$(RM) -rf $(OBJ_DIR)/*.o
	$(RM) $(targets)
