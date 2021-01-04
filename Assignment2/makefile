COMPILER = gcc
MPI_COMPILER = mpicc 
OMP_COMPILER = gcc -fopenmp 
GEN_EXE_NAME = generate_pgm.x
SERIAL_EXE_NAME = blur_serial.x
MPI_EXE_NAME = blur_mpi.x
OMP_EXE_NAME = blur_omp.x
FLAGS = -Wall -Wextra  -march=native -O3 -g3
INCLUDE = -I src/include/
CFLAGS=

CHK_OBJ = src/check_pgm.o  src/common/pgm_lib.o
GEN_OBJ = src/generate_pgm.o src/common/pgm_lib.o
COMMON_OBJ = src/common/blur_common_lib.o  src/common/pgm_lib.o src/common/kernel_lib.o 
SERIAL_OBJ = src/blur_serial_main.o src/common/blur_pgm_lib.o
SERIALTEST_OBJ = src/blur_serial_main_perftest.o src/common/blur_pgm_lib.o
OMP_OBJ = src/blur_omp_main.o src/common/blur_pgm_halo_lib.o src/common/img_cell.o
MPI_OBJ = src/blur_mpi_main.o src/common/blur_pgm_halo_lib.o src/common/img_cell.o




MOD= *.mod

default: serial

all : serial mpi omp

gen:  $(GEN_OBJ)
	$(COMPILER) $^ -o $(GEN_EXE_NAME) $(FLAGS) -lm

serial: $(COMMON_OBJ)  $(SERIAL_OBJ)
	$(COMPILER) $^ -o $(SERIAL_EXE_NAME) $(FLAGS) -lm
	
serialtest: $(COMMON_OBJ)  $(SERIALTEST_OBJ)
	$(COMPILER) $^ -o $(SERIAL_EXE_NAME) $(FLAGS) -lm
	
	
omp: $(COMMON_OBJ) $(OMP_OBJ) 
	$(OMP_COMPILER) $^ -o $(OMP_EXE_NAME) $(FLAGS) -lm

mpi: $(COMMON_OBJ) $(MPI_OBJ)
	$(MPI_COMPILER) $^ -o $(MPI_EXE_NAME) $(FLAGS) -lm
	
	
check: $(CHK_OBJ)
	$(COMPILER) $^ -o check_pgm.x $(INCLUDE) $(FLAGS) -lm

clean:
	rm -rf $(SERIAL_OBJ) $(SERIALTEST_OBJ) $(COMMON_OBJ) $(OMP_OBJ) $(MPI_OBJ) $(CHK_OBJ)
	rm -rf $(SERIAL_EXE_NAME) $(OMP_EXE_NAME) $(MPI_EXE_NAME) $(GEN_EXE_NAME) check_pgm.x 
	rm -rf $(MOD)
	
.PHONY: clean default



src/blur_omp_main.o : src/blur_omp_main.c
	$(OMP_COMPILER) -c $< -o $@ $(FLAGS) $(CFLAGS) $(INCLUDE)
	
src/blur_mpi_main.o : src/blur_mpi_main.c
	$(MPI_COMPILER) -c $< -o $@ $(FLAGS) $(CFLAGS) $(INCLUDE)
	
%.o: %.c
	$(COMPILER) -c $< -o $@ $(FLAGS) $(CFLAGS) $(INCLUDE)


	

main.o: main.c
node.o: src/node.c