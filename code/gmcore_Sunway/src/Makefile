.PHONY:all

all:$(OBJ_DIR)/master.o $(OBJ_DIR)/physical_variable.o $(OBJ_DIR)/slave.o

$(OBJ_DIR)/master.o:master.c
	$(CC) $(CM) $(MOPT) -c $^ -o $@ $(CFLAG) 

$(OBJ_DIR)/slave.o:slave.c
	$(SCC) $(CS) -c $^ -o $@ $(CFLAG) $(SOPT)

$(OBJ_DIR)/physical_variable.o:physical_variable.c
	$(CC) -c $^ -o $@ $(CFLAG)  




# .PHONY:all

# all:$(OBJ_DIR)/main.o $(OBJ_DIR)/physical_variable.o

# $(OBJ_DIR)/main.o:main.c
# 	$(CC) -c $^ -o $@ $(CFLAG)

# $(OBJ_DIR)/physical_variable.o:physical_variable.c
# 	$(CC) -c $^ -o $@ $(CFLAG)
