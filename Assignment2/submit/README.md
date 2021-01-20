# Assignment 2 - Parallel Blurring algorithm

## Compilation

The makefile offers three targets: 

- **serial**, which uses gcc to produce a serialised version of the code, will produce the execuable **blur.x**
- **omp**, which uses gcc -fopenmp to produce a multithreaded version of the code **blur.omp.x**
- **mpi**, which uses mpicc to produce a multi-process version of the code **blur.mpi.x**

They may all be compiled with make all.
Several compilation flags may be supplied to the make command to customise the program:

- **CFLAGS= '-DTIME'** will enable the timing functionality, times for sections of the code will be measured and printed to terminal at the end
- **CFLAGS= '-DINFO'** will toggle debug information on the status of the program, in the parallelised code this will also output information on the grid arrangements of processes/threads and the size of the image subdivisions

Flags may be supplied to select the blurring function to be used:

- **CFLAGS= '-DBL_COPY'** will duplicate the image buffer to store the blurred pixels during the operation (serial only)
- **CFLAGS= '-DBL_LINEBUF'** will use a smaller buffer option to reduce memory impact (serial only)
- **CFLAGS= '-DBL_UNROL2'** will apply x2 unrolling to the convolution loops
- **CFLAGS= '-DBL_UNROL4'** will apply x4 unrolling to the convolution loops
- **CFLAGS= '-DBL_UNROL8'** will apply x8 unrolling to the convolution loops

if no flag is supplied, x4 unrolling is selected by default.

## Running

The standard parameters are:

**exe_name -input input_image .pgm -kernel-type type -kernel-size size_w [ size_h ]**
      **\[ -output output_image.pgm ]  \[-kernel-weight weight ] \[ -kernel-file filename ]**
