COMPILER = gcc
MPI_COMPILER = mpicc 
EXE_NAME = blur.x
FLAGS = -Wall -Wextra  -march=native -O3
INCLUDE = -I src/include/

COMMON_OBJ = src/common/blur_common_lib.o  src/common/pgm_lib.o src/common/kernel_lib.o 
SERIAL_OBJ = src/serial/blur_serial_main.o src/serial/pgm_blur_serial.o
SERIALTEST_OBJ = src/serial/blur_serial_main_perftest.o src/serial/pgm_blur_serial.o
OPENMP_OBJ =
MPI_OBJ = src/mpi/blur_mpi_main.o src/mpi/blur_mpi_lib.o




MOD= *.mod

default: serial

serial: $(COMMON_OBJ)  $(SERIAL_OBJ)
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm
	
serialtest: $(COMMON_OBJ)  $(SERIALTEST_OBJ)
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm
	
	
openmp: $(COMMON_OBJ) $(OPENMP_OBJ) 
	$(COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm

mpi: $(COMMON_OBJ) $(MPI_OBJ)
	$(MPI_COMPILER) $^ -o $(EXE_NAME) $(FLAGS) -lm

clean:
	rm -rf $(SERIAL_OBJ) $(SERIALTEST_OBJ) $(COMMON_OBJ) $(OPENMP_OBJ) $(MPI_OBJ) $(EXE_NAME) $(MOD)
	
.PHONY: clean default



	

	
$(MPI_EXE_NAME): $(MPI_OBJ) $(COMMON_OBJ)
	$(COMPILER) $^ -o $@ $(FLAGS) -lm
	
src/mpi/blur_mpi_main.o : src/mpi/blur_mpi_main.c
	$(MPI_COMPILER) -c $< -o $@ $(FLAGS) $(INCLUDE)
	
%.o: %.c
	$(COMPILER) -c $< -o $@ $(FLAGS) $(INCLUDE)


	

main.o: main.c
node.o: src/node.c