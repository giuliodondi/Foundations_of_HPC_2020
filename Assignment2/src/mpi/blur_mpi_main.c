#include <pgm.h>
#include <kernel_t.h>
#include <img_cell.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include <mpi.h>
#define USE MPI

//blurring function headers
void get_cell_1D( int nprocs, int proc_id, img_cell* proc_cell, int width, int height, char img_bytes, int halowidth) ;
void pgm_blur_serial(  pgm* input_img , kernel_t* k);
void pgm_blur_halo(  pgm* input_img , kernel_t* k,  const char* halos);



int main( int argc, char **argv ) 
{ 
	//mpi common variables
	int nprocs, proc_id;
	int collecttag = 123;

	//initialise the MPI communicators
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	
	//variables common to all processes;
	double start_t, elapsed;
	pgm  local_image;
	pgm  original_image ;
	kernel_t kernel_ptr;
	char img_bytes;
	img_cell proc_cell;
	int halowidth;
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";

	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel_ptr) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	
	
	
	
	if (read_pgm( &original_image , infile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}

	img_bytes  = (1 + (original_image.maxval > 255));
	halowidth = (kernel_ptr.size - 1)/2;

	
	
	
	//build the cells accountign for the halo above and below the cells
	get_cell_1D( nprocs, proc_id, &proc_cell, original_image.width, original_image.height, img_bytes, halowidth);
	
	printf("Process %d will process %d rows, %d cols of the image.\n", proc_id, proc_cell.height, proc_cell.width );
	
	//initialise the working image
	//get the pointer to the beginnig of the memory section
	local_image.width = proc_cell.width;
	local_image.height = proc_cell.height;
	local_image.maxval = original_image.maxval;
	local_image.data = &original_image.data[proc_cell.idx];
		
	start_t = MPI_Wtime();
	pgm_blur_halo( &local_image, &kernel_ptr , proc_cell.halos);
	elapsed = MPI_Wtime() - start_t;

	printf("Wall time for process %d : %f s.\n",proc_id,elapsed);
	
	
		
	if (proc_id == 0) {

		//retrieve data directly in the right place
		//update the now-dummy cell variable with the information about the process we're receiving from
		for (int p=1;p < nprocs; ++p) {
			get_cell_1D( nprocs, p, &proc_cell, original_image.width, original_image.height, img_bytes, 0);
			MPI_Recv(&original_image.data[proc_cell.idx] , proc_cell.size , MPI_UINT8_T, p, collecttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			printf("Master has received from process %d.\n", p);
		}
		
		
		//write the image
		if ( write_pgm( &original_image, outfile)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
	   	printf("Master has written output file \"%s\".\n",outfile);

	}
	else {
		
		//prepare to send back the data without the halos
		//update the cells properties removing the halos since we don't send them back
		get_cell_1D( nprocs, proc_id, &proc_cell, original_image.width, original_image.height, img_bytes, 0);

		MPI_Send(&original_image.data[proc_cell.idx], proc_cell.size, MPI_UINT8_T, 0, collecttag, MPI_COMM_WORLD);
		
	}
		
	
	

	clear_pgm( &original_image);
	//the local image points to the memory location of the original image which is free
	local_image.data = NULL;
	delete_kernel(&kernel_ptr);
	clear_pgm( &local_image);
	MPI_Finalize();
	return 0;
	

} 