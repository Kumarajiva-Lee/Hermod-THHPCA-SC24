.PHONY:all

OPT = -Dslave_acc

#ifdef slave_acc 
#all:$(OBJ_DIR)/filter.o $(OBJ_DIR)/filter_slave.o
all:$(OBJ_DIR)/filter.o

$(OBJ_DIR)/filter.o:filter.c
	$(CC) -c $^ -o $@ $(CFLAG) $(OPT)

#$(OBJ_DIR)/filter_slave.o:filter_slave.c
#	$(SCC) $(CS) -c $^ -o $@ $(CFLAG) $(SOPT)
