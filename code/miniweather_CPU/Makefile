##
CC=mpiicc
#OPT = -O3 -ftz -xhost
OPT = -O3 -D A_LOCAL -D ALIGNED -D NOALIAS -ftz -IPF_fltacc -IPF_fma -IPF-fp-relaxed -xhost -restrict -qopt-prefetch -qopt-multi-version-aggressive
#OPT = -O3 -D A_LOCAL -D ALIGNED -D NOALIAS -restrict -xHost -D VEC -fast

# CFLAGS = -Wall -g -debug full -DGETTIMEOFDAY -std=c99 $(OPT)
CFLAG = -Wall -g  -traceback -DGETTIMEOFDAY -std=c99 $(OPT)
LDFLAGS = -Wall -g -lnetcdf

##
PWD_DIR=$(shell pwd)
LIB_DIR=$(PWD_DIR)/LIB
SRC_DIR=$(PWD_DIR)/src
OBJ_DIR=$(PWD_DIR)/obj

#TEST_DIR=$(PWD_DIR)/test

# INC_DIR=$(PWD_DIR)/inc

targets = app



##
export PWD_DIR OBJ_DIR LIB_DIR SRC_DIR CC CFLAG LDFLAGS targets

##
all:
	make -C $(LIB_DIR)
	make -C $(SRC_DIR)
	make -C $(OBJ_DIR)

##
clean:
	$(RM) -rf $(OBJ_DIR)/*.o
	$(RM) $(targets)
