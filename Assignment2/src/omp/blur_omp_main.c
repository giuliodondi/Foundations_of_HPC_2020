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
void get_cell_1D(const int nprocs, const int proc_id, img_cell* proc_cell, const pgm* image, const unsigned int* halowidth);
int trim_halo_1D( const img_cell* proc_cell, const char img_bytes, const unsigned int* halowidth );
char read_write_cell_1D( pgm* original_img, pgm* local_img, img_cell* proc_cell, const unsigned int* halowidth, const char* mode);
void pgm_blur_halo(  pgm* input_img , kernel_t* k,  const char* halos);




int main( int argc, char **argv ) 
{ 
	#ifdef TIME
	double header_time, read_time, write_time, global_t, avg_buf_read_t=0, avg_blur_time_t=0, avg_buf_write_t=0 ;
	global_t = omp_get_wtime();
	#endif
		
	pgm  original_image = new_pgm();
	kernel_t kernel;
	unsigned int halowidth0[2] = {0,0};
	
	long int header_offs=0;
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";

	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	
	
	
	header_time = omp_get_wtime();
	//read the file header
	if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
	header_time = omp_get_wtime() - header_time;
	#endif

	

	//allocate memory for the image
	if (allocate_pgm_memory( &original_image)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}


	
	#pragma omp parallel  shared( kernel, original_image) reduction(+: avg_buf_read_t, avg_blur_time_t, avg_buf_write_t)
	{

		int nprocs = omp_get_num_threads();
		int proc_id =  omp_get_thread_num();
		
		img_cell proc_cell;
		pgm  local_image = new_pgm();
		double buf_read_time, blur_time, buf_write_time;
		
		
		
		//touch the image buffer
		//first we need information on the local image cells without halos
		//build the cells accountign for the halo above and below the cells
		get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, halowidth0);

		
		//touch the memory
		{
			//wortk out the beginning of the buffer in the original image
			register int img_idx;// = ( original_img->size[0]*proc_cell->idx[1]+ proc_cell->idx[0])*original_img->pix_bytes;
			for (int j=proc_cell.idx[0]; j< (proc_cell.idx[0] + proc_cell.size[0]); ++j) {
			for (int i=proc_cell.idx[1]; i< (proc_cell.idx[1] + proc_cell.size[1]); ++i) {
				
					//printf("proc %d i j (%d , %d)\n",proc_id,i,j);
					img_idx =  (original_image.size[0]*i + j ) *original_image.pix_bytes;
					original_image.data[ img_idx ] =0;
				}
			}
		}
		
		
		
		//only the master reads the image
		#pragma omp master 
		{
			#ifdef TIME
			read_time = omp_get_wtime();
			#endif
			//read the image data
			if (read_pgm_data( &original_image , infile, &header_offs)== -1 ) {
				printf("Aborting.\n");
				clear_pgm( &original_image);
				delete_kernel(&kernel);
				//return -1;
			}
			#ifdef TIME
			read_time = omp_get_wtime() - read_time;
			#endif
			
			
			 #ifdef INFO
				printf("Input file \"%s\" has been read.\n",infile);
				printf("The image is %d x %d.\n",original_image.size[0],original_image.size[1]);
			#endif
		
		}
		
		#pragma omp barrier
	
		
		
		//build the cells accountign for the halo above and below the cells
		get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, kernel.halfsize);
		
		//initialise the working image
		//get the pointer to the beginnig of the memory section
		local_image.size[0] = proc_cell.size[0];
		local_image.size[1] = proc_cell.size[1];
		local_image.maxval = original_image.maxval;
		local_image.pix_bytes = original_image.pix_bytes;
		
		#ifdef TIME
		buf_read_time = omp_get_wtime();
		#endif
		//allocate local buffer for the image portion
		if (read_write_cell_1D(&original_image, &local_image, &proc_cell, kernel.halfsize,"r")== -1 ) {
			printf("Thread %d couldn't allocate memory for the image buffer.\n",proc_id);
			clear_pgm( &local_image);
			if (proc_id==0) {
				clear_pgm( &original_image);
			}
			//return 0;
		}
		#ifdef TIME
		avg_buf_read_t += omp_get_wtime() - buf_read_time;
		#endif

		
		#ifdef INFO
		printf("Thread %d will process %d rows, %d cols of the image.\n", proc_id, proc_cell.size[1], proc_cell.size[0] );
		#endif
		
		//create a local copy of the kernel
		kernel_t local_kernel;
		
		if (copy_kernel(&local_kernel, &kernel)== -1 ) {
			printf("Thread %d couldn't allocate memory for the local kernel.\n",proc_id);
			//return -1;	
		} 
		
		
		#ifdef TIME
		blur_time = omp_get_wtime();
		#endif
		
		pgm_blur_halo( &local_image, &local_kernel , proc_cell.halos);
		
		#ifdef TIME
		avg_blur_time_t += omp_get_wtime() - blur_time;
		#endif
		
		
		
		//update the cell parameters setting zero halosize
		get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, halowidth0);
		
		//collect the results now
		#ifdef TIME
		buf_write_time = omp_get_wtime();
		#endif
		read_write_cell_1D(&original_image, &local_image, &proc_cell, kernel.halfsize, "w");
		#ifdef TIME
		avg_buf_write_t += omp_get_wtime() - buf_write_time;
		#endif
		
		clear_pgm( &local_image);
		delete_kernel(&local_kernel);
	}

	
	

	header_time = omp_get_wtime();
	//write the file header
	if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	header_time = omp_get_wtime() - header_time;
	
	
	write_time = omp_get_wtime();
	//write the image data
	if (write_pgm_data( &original_image , outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	write_time = omp_get_wtime() - write_time;
	
	
	#ifdef INFO
	printf("Master has written output file \"%s\".\n",outfile);
	#endif

	clear_pgm( &original_image);
	//the local image points to the memory location of the original image which is free
	delete_kernel(&kernel);
	
	#ifdef TIME
	global_t = omp_get_wtime() - global_t;
	
	int nprocs = atoi(getenv ("OMP_NUM_THREADS"));
	
	avg_buf_read_t = avg_buf_read_t/nprocs;
	avg_blur_time_t = avg_blur_time_t/nprocs;
	avg_buf_write_t = avg_buf_write_t/nprocs;
	
	printf("Header read time 	: %f s.\n",header_time);
	printf("Data read time		: %f s.\n",read_time);
	printf("Avg buf read time 	: %f s.\n",avg_buf_read_t);
	printf("Avg blur time 		: %f s.\n",avg_blur_time_t);
	printf("Avg buf write time 	: %f s.\n",avg_buf_write_t);
	printf("Header write time 	: %f s.\n",header_time);
	printf("Data write time 	: %f s.\n",write_time);
	printf("Total time 			: %f s.\n",global_t);
	#endif
	return 0;
	

} 