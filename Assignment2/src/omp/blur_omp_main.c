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

#include <omp.h>


//blurring function headers
void get_cell_1D( int nprocs, int proc_id, img_cell* proc_cell, int width, int height, char img_bytes, int halowidth) ;
int trim_halo_1D( img_cell* proc_cell, char img_bytes, int halowidth );
char read_write_cell_1D( pgm* original_img, pgm* local_img, img_cell* proc_cell, char img_bytes, int halowidth, char* mode);
void pgm_blur_serial(  pgm* input_img , kernel_t* k);
void pgm_blur_halo(  pgm* input_img , kernel_t* k,  const char* halos);



int main( int argc, char **argv ) 
{ 
	
	
		
	pgm  original_image ;
	kernel_t kernel;
	char img_bytes;
	int halowidth;
	double global_start_t, global_elapsed;
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";

	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	
	
	
	
	if (read_pgm( &original_image , infile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}

	img_bytes  = (1 + (original_image.maxval > 255));
	halowidth = (kernel.size - 1)/2;

	global_start_t = omp_get_wtime();
	#pragma omp parallel shared( kernel, original_image, halowidth, img_bytes)
	{
		
		int nprocs = omp_get_num_threads();
		int proc_id =  omp_get_thread_num();
		
		//create a local copy of the kernel
		kernel_t local_kernel;
		
		
		if (copy_kernel(&local_kernel, &kernel)== -1 ) {
			printf("Thread %d couldn't allocate memory for the local kernel.\n",proc_id);
			//return -1;	
		} 
		
		img_cell proc_cell;
		pgm  local_image;
		double start_t, elapsed;
		
		//build the cells accountign for the halo above and below the cells
		get_cell_1D( nprocs, proc_id, &proc_cell, original_image.width, original_image.height, img_bytes, halowidth);

		printf("Thread %d will process %d rows, %d cols of the image.\n", proc_id, proc_cell.height, proc_cell.width );

		//initialise the working image
		//get the pointer to the beginnig of the memory section
		local_image.width = proc_cell.width;
		local_image.height = proc_cell.height;
		local_image.maxval = original_image.maxval;
		
		if (read_write_cell_1D(&original_image, &local_image, &proc_cell, img_bytes, halowidth,"r")== -1 ) {
			printf("Thread %d couldn't allocate memory for the image buffer.\n",proc_id);
			clear_pgm( &local_image);
			if (proc_id==0) {
				clear_pgm( &original_image);
				delete_kernel(&local_kernel);
			}
			//return 0;
		}
		
		//allocate local buffer for the image portion
		
		start_t = omp_get_wtime();
		pgm_blur_halo( &local_image, &local_kernel , proc_cell.halos);
		elapsed = omp_get_wtime() - start_t;
		
		printf("Wall time for Thread %d : %f s.\n",proc_id,elapsed);
		
		//collect the results now
		//there should be no need to put a barrier since we never write to the same place
		
		//update the cell parameters setting zero halosize
		get_cell_1D( nprocs, proc_id, &proc_cell, original_image.width, original_image.height, img_bytes, 0);
		
		read_write_cell_1D(&original_image, &local_image, &proc_cell, img_bytes, halowidth, "w");
		
		clear_pgm( &local_image);
		delete_kernel(&local_kernel);
	}

	global_elapsed = omp_get_wtime() - global_start_t;
	printf("Wall time for parallel region: %f s.\n",global_elapsed);
	
	//write the image
	if ( write_pgm( &original_image, outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	printf("Master has written output file \"%s\".\n",outfile);
	

	clear_pgm( &original_image);
	//the local image points to the memory location of the original image which is free
	delete_kernel(&kernel);
	return 0;
	

} 