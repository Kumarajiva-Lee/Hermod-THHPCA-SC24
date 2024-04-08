.PHONY:all

all:$(OBJ_DIR)/main.o $(OBJ_DIR)/physical_variable.o

$(OBJ_DIR)/main.o:main.c
	$(CC) -c $^ -o $@ $(CFLAG)

$(OBJ_DIR)/physical_variable.o:physical_variable.c
	$(CC) -c $^ -o $@ $(CFLAG)
