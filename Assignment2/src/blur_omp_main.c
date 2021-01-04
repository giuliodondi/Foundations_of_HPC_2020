#include <pgm.h>
#include <kernel_t.h>
#include <p_grid.h>
#include <img_cell.h>
#include <common_headers.h>
#include <blur_pgm.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include <omp.h>





int main( int argc, char **argv ) 
{ 
	#ifdef TIME
	double avg_buf_read_t=0, avg_blur_time_t=0, avg_buf_write_t=0 ;
	double avg_buf_read_t2=0, avg_blur_time_t2=0, avg_buf_write_t2=0 ;
	double header_time, read_time, write_time, total_t;
	total_t = omp_get_wtime();
	#endif
		
	pgm  original_image = new_pgm();
	kernel_t kernel;
	
	
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
	
	
	#ifdef TIME
	header_time = omp_get_wtime();
	#endif
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


	#ifdef TIME
	#pragma omp parallel  shared( kernel, original_image) reduction(+: avg_buf_read_t, avg_blur_time_t, avg_buf_write_t, avg_buf_read_t2, avg_blur_time_t2, avg_buf_write_t2)
	#else
	#pragma omp parallel  shared( kernel, original_image)
	#endif
	
	{
		int nprocs = omp_get_num_threads();
		int proc_id =  omp_get_thread_num();
		
		p_grid grid;
		build_grid(&grid,nprocs);
		//grid.size[0]= 1 ;
		//grid.size[1] = nprocs;
		
		#ifdef INFO
		if (proc_id==0) {
			printf("Threads arranged on a grid %d x %d\n",grid.size[0],grid.size[1]);
		}
		#endif
		
		
		img_cell cell_halo, cell_nohalo;
		pgm  local_image = new_pgm();
		unsigned int halowidth0[2] = {0,0};
		kernel_t local_kernel;
		
		#ifdef TIME
			double buf_read_time, blur_time, buf_write_time;
		#endif
		
		//create a local copy of the kernel
		if (copy_kernel(&local_kernel, &kernel)== -1 ) {
			printf("Thread %d couldn't allocate memory for the local kernel.\n",proc_id);
			//return -1;	
		} 
		
		
		//initialise cells with halo and without halo
		memcpy( cell_halo.coords , get_grid_coords(  &grid, proc_id ) , 2*sizeof(unsigned int) );
		memcpy( cell_nohalo.coords , cell_halo.coords , 2*sizeof(unsigned int) );
		
		get_cell_grid( &grid, &cell_halo, &original_image, local_kernel.halfsize);
		get_cell_grid( &grid, &cell_nohalo, &original_image, halowidth0);

		
		#ifdef INFO
		printf("\nCell %d has coords ( %d , %d ) \n", proc_id, cell_halo.coords[0] , cell_halo.coords[1] );
		printf("Cell %d  has size %d x %d , starts at img line %d col %d .\n", proc_id, cell_halo.size[0], cell_halo.size[1], cell_halo.idx[1], cell_halo.idx[0] );
		printf("Cell %d halos : (%d %d %d %d).\n", proc_id, cell_halo.halos[0], cell_halo.halos[1], cell_halo.halos[2], cell_halo.halos[3] );
		printf("\n");
		#endif
		
		
		
		//touch the image buffer with the nohalo parameters
		{
			//wortk out the beginning of the buffer in the original image
			register int img_idx;
			for (size_t i=cell_nohalo.idx[1]; i< (cell_nohalo.idx[1] + cell_nohalo.size[1]); ++i) {
			for (size_t j=cell_nohalo.idx[0]; j< (cell_nohalo.idx[0] + cell_nohalo.size[0]); ++j) {				
					//printf("proc %d i j (%d , %d)\n",proc_id,i,j);
					img_idx =  (original_image.size[0]*i + j ) *original_image.pix_bytes;
					original_image.data[ img_idx ] =0;
				}
			}
		}
		
		#pragma omp barrier
		
		//only the master reads the image file
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
		
		//otherwise the children processes use data which isn't loaded yet
		#pragma omp barrier
		
		
		//initialise the working image
		//get the pointer to the beginnig of the memory section
		local_image.size[0] = cell_halo.size[0];
		local_image.size[1] = cell_halo.size[1];
		local_image.maxval = original_image.maxval;
		local_image.pix_bytes = original_image.pix_bytes;
		
		
		
		local_image.data = (uint8_t*)malloc( cell_halo.size_*sizeof(uint8_t) );
		if ( ! local_image.data) {
			printf("Error allocating memory for a cell.\n");
			clear_pgm( &local_image);
			if (proc_id==0) {
				clear_pgm( &original_image);
			}
			//return -1;
		}
		
		
		#ifdef TIME
		buf_read_time = omp_get_wtime();
		#endif
		
		read_img_buffer( &original_image , &local_image, &cell_halo);

			
		#ifdef TIME
		buf_read_time = omp_get_wtime() - buf_read_time;
		avg_buf_read_t += buf_read_time;
		avg_buf_read_t2 += buf_read_time*buf_read_time;
		#endif


		
		#ifdef TIME
		blur_time = omp_get_wtime();
		#endif
		
		blur_halo_func_manager( &local_image, &local_kernel , cell_halo.halos);
		
		#ifdef TIME
		blur_time = omp_get_wtime() - blur_time;
		avg_blur_time_t += blur_time;
		avg_blur_time_t2 += blur_time*blur_time;
		#endif
		

		
		#pragma omp barrier
		//collect the results back into the image buffer
		
		#ifdef TIME
		buf_write_time = omp_get_wtime();
		#endif
		
		
		write_img_buffer( &original_image , &local_image, &cell_halo, &cell_nohalo);
		
		
		#ifdef TIME
		buf_write_time = omp_get_wtime() - buf_write_time;
		avg_buf_write_t += buf_write_time;
		avg_buf_write_t2 += buf_write_time*buf_write_time;
		#endif
		
		#pragma omp barrier
		
		clear_pgm( &local_image);
		delete_kernel(&local_kernel);
	}
	delete_kernel(&kernel);
	
	
	#ifdef TIME
	header_time = omp_get_wtime();
	#endif
	//write the file header
	if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
	header_time = omp_get_wtime() - header_time;
	write_time = omp_get_wtime();
	#endif
	//write the image data
	if (write_pgm_data( &original_image , outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
	write_time = omp_get_wtime() - write_time;
	#endif
	
	
	#ifdef INFO
	printf("Master has written output file \"%s\".\n",outfile);
	#endif

	clear_pgm( &original_image);
	//the local image points to the memory location of the original image which is free
	
	#ifdef TIME
	total_t = omp_get_wtime() - total_t;
	
	int nprocs = atoi(getenv ("OMP_NUM_THREADS"));
	
	avg_buf_read_t = avg_buf_read_t/nprocs;
	avg_blur_time_t = avg_blur_time_t/nprocs;
	avg_buf_write_t = avg_buf_write_t/nprocs;
	
	if ( nprocs==1 ) {
		avg_buf_read_t2 = 0;
		avg_blur_time_t2 = 0;
		avg_buf_write_t2 = 0;
	} else {
		avg_buf_read_t2 = sqrt(avg_buf_read_t2/nprocs - avg_buf_read_t*avg_buf_read_t);
		avg_blur_time_t2 = sqrt(avg_blur_time_t2/nprocs - avg_blur_time_t*avg_blur_time_t);
		avg_buf_write_t2 = sqrt(avg_buf_write_t2/nprocs - avg_buf_write_t*avg_buf_write_t);
	}
	
	printf("Header read time 	: %f s.\n",header_time);
	printf("Data read time		: %f s.\n",read_time);
	printf("Avg buf read time 	: %f +- %f s.\n",avg_buf_read_t, avg_buf_read_t2 );
	printf("Avg blur time 		: %f +- %f s.\n",avg_blur_time_t, avg_blur_time_t2 );
	printf("Avg buf write time 	: %f +- %f s.\n",avg_buf_write_t, avg_buf_write_t2 );
	printf("Header write time 	: %f s.\n",header_time);
	printf("Data write time 	: %f s.\n",write_time);
	printf("Total time 			: %f s.\n",total_t);
	#endif
	return 0;
	

} 