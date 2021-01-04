#include <pgm.h>
#include <kernel_t.h>
#include <img_cell.h>
#include <common_headers.h>
#include <blur_pgm.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include <mpi.h>
#define USE MPI



int main( int argc, char **argv ) 
{ 
	
	#ifdef TIME
	double read_time=0, blur_time=0, write_time=0 ;
	double read_time2=0, blur_time2=0, write_time2=0 ;
	double header_read_time, header_write_time,  total_t, total_t2;
	total_t = MPI_Wtime();
	#endif
	
	//mpi common variables
	int nprocs, proc_id;
	MPI_Status status;
	MPI_Offset my_data_offs;

	//initialise the MPI communicators
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	
	p_grid grid;
	build_grid(&grid,nprocs);

	#ifdef INFO
	if (proc_id==0) {
		printf("Threads arranged on a grid %d x %d\n",grid.size[0],grid.size[1]);
	}
	#endif


	

	pgm  local_image = new_pgm();
	pgm  original_image = new_pgm();
	kernel_t kernel;
	img_cell cell_halo, cell_nohalo;
	long int header_offs=0;
	unsigned int halowidth0[2] = {0,0};
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "";
	
	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		delete_kernel(&kernel);
		return -1;
		MPI_Finalize();
	}
	
	
	//only master reads the file header
	if (proc_id == 0) {
		#ifdef TIME
		header_read_time = MPI_Wtime();
		#endif
		//read the file header
		if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			delete_kernel(&kernel);
			MPI_Finalize();
			return -1;
		}
		#ifdef TIME
		header_read_time = MPI_Wtime() - header_read_time;
		#endif
		 #ifdef INFO
			printf("Input file \"%s\" has been opened.\n",infile);
			printf("The image is %d x %d.\n",original_image.size[0],original_image.size[1]);
		#endif
	}
	
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	
	//send info on the incoming image metadata
	MPI_Bcast(&original_image.size, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.maxval, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.pix_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	

	
	
	//initialise cells with halo and without halo
	memcpy( cell_halo.coords , get_grid_coords(  &grid, proc_id ) , 2*sizeof(unsigned int) );
	memcpy( cell_nohalo.coords , cell_halo.coords , 2*sizeof(unsigned int) );

	get_cell_grid( &grid, &cell_halo, &original_image, kernel.halfsize);
	get_cell_grid( &grid, &cell_nohalo, &original_image, halowidth0);
	
	
	#ifdef INFO
	printf("\nCell %d has coords ( %d , %d ) \n", proc_id, cell_halo.coords[0] , cell_halo.coords[1] );
	printf("Cell %d  has size %d x %d , starts at img line %d col %d .\n", proc_id, cell_halo.size[0], cell_halo.size[1], cell_halo.idx[1], cell_halo.idx[0] );
	printf("Cell %d halos : (%d %d %d %d).\n", proc_id, cell_halo.halos[0], cell_halo.halos[1], cell_halo.halos[2], cell_halo.halos[3] );
	printf("\n");
	#endif

	//initialise the local image
	local_image.size[0] = cell_halo.size[0];
	local_image.size[1] = cell_halo.size[1];
	local_image.maxval = original_image.maxval;
	local_image.pix_bytes = original_image.pix_bytes;
	//allocte memory for the local image
	local_image.data = (uint8_t*)malloc( cell_halo.size_*sizeof(uint8_t) );
	if ( ! local_image.data) {
		printf("Error allocating memory for a cell.\n");
		delete_kernel(&kernel);
		clear_pgm( &local_image);
		if (proc_id==0) {
			clear_pgm( &original_image);
		}
		MPI_Finalize();
		return -1;
	}
	
	
	
	
	//custom mpi datatype for the pixel
	MPI_Datatype pixel;
	
	if (original_image.pix_bytes == 2) {
		pixel= MPI_UINT16_T;
	} else {
		pixel = MPI_UINT8_T;
	}
	MPI_Type_commit ( &pixel );
	
	
	//parallel file read
	my_data_offs = header_offs;
	
	//create subarray describing the data with halos within the original image
	int dtype_size[2];
	int dtype_subsize[2];
	int dtype_idx[2];

	
	dtype_size[0] = original_image.size[1];
	dtype_size[1] = original_image.size[0];
	dtype_subsize[0] = cell_halo.size[1];
	dtype_subsize[1] = cell_halo.size[0];
	dtype_idx[0] = cell_halo.idx[1];
	dtype_idx[1] = cell_halo.idx[0];
	
	MPI_Datatype img_subarr_halo;	
	MPI_Type_create_subarray(2, dtype_size, dtype_subsize, dtype_idx ,MPI_ORDER_C, pixel, &img_subarr_halo);
	MPI_Type_commit(&img_subarr_halo);
	
	
	#ifdef TIME
	read_time = MPI_Wtime();
	#endif
	
	MPI_File in_file;
	MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY,MPI_INFO_NULL, &in_file);
	MPI_File_set_view(in_file, my_data_offs, pixel , img_subarr_halo, "native", MPI_INFO_NULL);
	MPI_File_read(in_file, &local_image.data[0], cell_halo.size_, pixel,  &status);
	MPI_File_close(&in_file);
	
	endian_swap(&local_image);
	
	#ifdef TIME
	read_time = MPI_Wtime() - read_time;
	read_time2 += read_time*read_time;
	#endif
	
	#ifdef INFO
		printf("Process %d has read the buffer.\n",proc_id);
	#endif
	
	 MPI_Type_free(&img_subarr_halo);
	
	//blurring
	#ifdef TIME
	blur_time = MPI_Wtime();
	#endif
	
	
	blur_halo_func_manager( &local_image, &kernel , cell_halo.halos);
	
	
	#ifdef TIME
	blur_time = MPI_Wtime() - blur_time;
	blur_time2 += blur_time*blur_time;
	#endif
	
	delete_kernel(&kernel);
	
	
	if (proc_id == 0) {
		#ifdef TIME
		header_write_time = MPI_Wtime();
		#endif
		//write the file header
		if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			clear_pgm( &local_image);
			delete_kernel(&kernel);
			MPI_Finalize();
			return -1;
		}
		#ifdef TIME
		header_write_time = MPI_Wtime() - header_write_time;
		#endif
	}
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	
	
	my_data_offs = header_offs;
	
	//create subarray describing the data without halos within the original image
	dtype_size[0] = original_image.size[1];
	dtype_size[1] = original_image.size[0];
	dtype_subsize[0] = cell_nohalo.size[1];
	dtype_subsize[1] = cell_nohalo.size[0];
	dtype_idx[0] = cell_nohalo.idx[1];
	dtype_idx[1] = cell_nohalo.idx[0];
	
	MPI_Datatype img_subarr_nohalo;
	MPI_Type_create_subarray(2, dtype_size, dtype_subsize, dtype_idx ,MPI_ORDER_C, pixel, &img_subarr_nohalo);
	MPI_Type_commit(&img_subarr_nohalo);
	
	//create subarray describing the data without halos within the local image buffer (which contains halos)
	dtype_size[0] = cell_halo.size[1];
	dtype_size[1] = cell_halo.size[0];
	dtype_subsize[0] = cell_nohalo.size[1];
	dtype_subsize[1] = cell_nohalo.size[0];
	dtype_idx[0] = cell_halo.halos[1];
	dtype_idx[1] = cell_halo.halos[0];
	
	MPI_Datatype cell_subarr_nohalo;
	MPI_Type_create_subarray(2, dtype_size, dtype_subsize, dtype_idx ,MPI_ORDER_C, pixel, &cell_subarr_nohalo);
	MPI_Type_commit(&cell_subarr_nohalo);
	
	//parallel write
	#ifdef TIME
	write_time = MPI_Wtime();
	#endif
	
	endian_swap(&local_image);
	
	
	
	MPI_File out_file;
	MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &out_file);
	MPI_File_set_view(out_file, my_data_offs, pixel , img_subarr_nohalo, "native", MPI_INFO_NULL);
	MPI_File_write(out_file, &local_image.data[0], 1 , cell_subarr_nohalo,  &status);
	MPI_File_close(&out_file);
	

	#ifdef TIME
	write_time = MPI_Wtime() - write_time;
	write_time2 += write_time*write_time;
	#endif
	
	
	#ifdef INFO
	if (proc_id==0) {
		printf("Output file \"%s\" has been written.\n",outfile);
	}
	#endif
	
	
	MPI_Type_free(&img_subarr_nohalo);
	MPI_Type_free(&cell_subarr_nohalo);
	//MPI_Type_free(&pixel);
	
	clear_pgm( &original_image);
	clear_pgm( &local_image);
	
	#ifdef TIME
	total_t = MPI_Wtime() - total_t;
	total_t2 = total_t*total_t;

	
	double time_arr[] = {
					read_time,
					read_time2,
					blur_time,
					blur_time2,
					write_time,
					write_time2,
					total_t,
					total_t2
	};
	
	double avg_time_arr[] = {0,0,0,0,0,0,0,0};
	
	MPI_Reduce(&time_arr, &avg_time_arr, 8, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
	
	if (proc_id==0) {
		for (int i=0; i<8; ++i) {
			avg_time_arr[i] /= nprocs;	
		}
		
		for (int i=1; i<8; i+=2) {
			if ( nprocs==1 ) {
				avg_time_arr[i] = 0;
			} else {
				avg_time_arr[i] = sqrt(avg_time_arr[i] - avg_time_arr[i-1]*avg_time_arr[i-1]);
			}
		}
		
		printf("Header read time 	: %f s.\n",header_read_time);
		printf("Avg read time 		: %f +- %f s.\n",avg_time_arr[0], avg_time_arr[1] );
		printf("Avg blur time 		: %f +- %f s.\n",avg_time_arr[2], avg_time_arr[3] );
		printf("Header write time 	: %f s.\n",header_write_time);
		printf("Avg write time 		: %f +- %f s.\n",avg_time_arr[4], avg_time_arr[5] );
		printf("Avg total time 		: %f +- %f s.\n",avg_time_arr[6], avg_time_arr[7] );
	}
	
	
	#endif
	
	 
	MPI_Finalize();
	return 0;
	

} 
