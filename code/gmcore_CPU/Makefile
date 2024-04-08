##
CC=mpiicc
OPT = -O3


# CFLAGS = -Wall -g -debug full -DGETTIMEOFDAY -std=c99 $(OPT)
CFLAG = -Wall -g  -traceback -DGETTIMEOFDAY -std=c99 $(OPT)
LDFLAGS = -Wall -g -lnetcdf

##
PWD_DIR=$(shell pwd)
LIB_DIR=$(PWD_DIR)/LIB
EXTERN_LIB_DIR=$(PWD_DIR)/ExternLIB
SRC_DIR=$(PWD_DIR)/src
OBJ_DIR=$(PWD_DIR)/obj

#TEST_DIR=$(PWD_DIR)/test

# INC_DIR=$(PWD_DIR)/inc

targets = app



##
export PWD_DIR OBJ_DIR EXTERN_LIB_DIR LIB_DIR SRC_DIR CC CFLAG LDFLAGS targets

##
all:
	make -C $(LIB_DIR)
	make -C $(EXTERN_LIB_DIR)
	make -C $(SRC_DIR)
	make -C $(OBJ_DIR)

##
clean:
	$(RM) -rf $(OBJ_DIR)/*.o
	$(RM) $(targets)
