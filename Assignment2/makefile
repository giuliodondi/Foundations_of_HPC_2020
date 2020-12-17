COMPILER = gcc
EXE_NAME = blur.x
FLAGS = -Wall -Wextra  -O0
INCLUDE = -I src/include/

COMMON_OBJ = src/common/pgm_lib.o src/common/blur_common_lib.o src/common/kernel_lib.o
SERIAL_OBJ = src/blur_serial_main.o 
OPENMP_OBJ =
MPI_OBJ =

MOD= *.mod

default: serial

serial: $(COMMON_OBJ)  $(SERIAL_OBJ)
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm
	
openmp: $(COMMON_OBJ) $(OPENMP_OBJ) 
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm

mpi: $(COMMON_OBJ) $(MPI_OBJ)
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm

clean:
	rm -rf $(SERIAL_OBJ) $(COMMON_OBJ) $(OPENMP_OBJ) $(MPI_OBJ) $(SERIAL_EXE_NAME) $(MOD)
	
.PHONY: clean default



	

	
$(MPI_EXE_NAME): $(MPI_OBJ) $(COMMON_OBJ)
	$(COMPILER) $^ -o $@ $(FLAGS) -lm
	
%.o: %.c
	$(COMPILER) -c $< -o $@ $(FLAGS) $(INCLUDE)
	

main.o: main.c
node.o: src/node.c