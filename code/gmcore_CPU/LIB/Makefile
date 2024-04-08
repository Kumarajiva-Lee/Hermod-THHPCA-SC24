.PHONY:all

all:$(OBJ_DIR)/Process.o $(OBJ_DIR)/Communicator.o $(OBJ_DIR)/Memory.o $(OBJ_DIR)/NCop.o $(OBJ_DIR)/Time.o $(OBJ_DIR)/Diagnose.o

$(OBJ_DIR)/Process.o:Process.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/Communicator.o:Communicator.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/Memory.o:Memory.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/NCop.o:NCop.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/Time.o:Time.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/Diagnose.o:Diagnose.c
	$(CC) -c $^ -o $@ $(CFLAG)