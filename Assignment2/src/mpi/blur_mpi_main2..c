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
void get_cell_1D( int nprocs, int proc_id, img_cell* proc_cell, pgm* image, int halowidth) ;
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
	MPI_Status status;
	
	//variables common to all processes;
	double start_t, elapsed;
	pgm  local_image;
	pgm  original_image ;
	kernel_t kernel_ptr;
	img_cell proc_cell;
	int halowidth;
	long int header_offs=0;
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";

	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel_ptr) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	
	
	if (proc_id == 0) {
		//read the file header
		if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
	}
	
	//send info on the incoming image data
	MPI_Bcast(&original_image.size, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.maxval, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.pix_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	

	if (read_pgm_data( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	//send out the header offset
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	

	halowidth = (kernel_ptr.size - 1)/2;
	
	
	//to use mpi subarrays we need a cart network
	//build a column of processes for the moment
	//if we want to implement the process grid, do it here
	int mpi_dims[2] = {1,nprocs};
	int mpi_periods[2] = {0,0};
	MPI_Comm MPI_COMM_CART;
	
	MPI_Cart_create(
		MPI_COMM_WORLD,	
		2,	
		mpi_dims,	
		mpi_periods,
		0,
		&MPI_COMM_CART
	);	
	
	MPI_Comm_rank(MPI_COMM_CART, &proc_id);
	int proc_coords[2];
	MPI_Cart_coords(MPI_COMM_CART, proc_id, 2, proc_coords);
	
	
	
	
	//build the cells accounting for the halo above and below the cells
	get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, halowidth);
	
	printf("\nI am process %d, with coordinates (%d, %d), \nI will process %d rows, %d cols of the image\n", 
		   proc_id, proc_coords[0], proc_coords[1], proc_cell.size[1], proc_cell.size[0] );
	
	//initialise the working image
	//get the pointer to the beginnig of the memory section
	local_image.size[0] = proc_cell.size[0];
	local_image.size[1] = proc_cell.size[1];
	local_image.maxval = original_image.maxval;
	local_image.pix_bytes = original_image.pix_bytes;
	
	local_image.data = (uint8_t*)malloc( proc_cell.size_*sizeof(uint8_t) );
	
	//initialise data type for the pixel size and the cell data with halos
	MPI_Datatype mpi_cell_halo;
	/*
	MPI_Datatype mpi_pixel_type;
	if (original_image.maxval>255) {
		pixel_type = MPI_UINT16_T;
	} else {
		pixel_type = MPI_UINT8_T;
	}
	*/
	
	printf("From process %d, the image cell is %d x %d and located at (%d,%d) in the image.\n",
		   proc_id, proc_cell.size[0], proc_cell.size[1] , proc_cell.idx[0], proc_cell.idx[1] );
	
	
	MPI_Type_create_subarray (
		2, 
		original_image.size,
		proc_cell.size, 
		proc_cell.idx, 
		MPI_ORDER_C,
		MPI_UINT8_T,
		&mpi_cell_halo
	);
	MPI_Type_commit(&mpi_cell_halo);
	
	MPI_File in_file;
	MPI_File_open(MPI_COMM_CART, infile, MPI_MODE_RDONLY,MPI_INFO_NULL, &in_file);
	
	//set the view past the header
	MPI_File_set_view(in_file, header_offs ,MPI_UINT8_T, mpi_cell_halo, "native", MPI_INFO_NULL);
	//read into local buffer
	MPI_File_read(in_file, local_image.data , proc_cell.size_, MPI_UINT8_T, &status);
	
	MPI_File_close(&in_file);
	MPI_Type_free(&mpi_cell_halo);
	
	
	
	start_t = MPI_Wtime();
	pgm_blur_halo( &local_image, &kernel_ptr , proc_cell.halos);
	elapsed = MPI_Wtime() - start_t;

	printf("Wall time for process %d : %f s.\n",proc_id,elapsed);
	
	
		
	if (proc_id == 0) {

		//retrieve data directly in the right place
		//update the now-dummy cell variable with the information about the process we're receiving from
		for (int p=1;p < nprocs; ++p) {
			get_cell_1D( nprocs, p, &proc_cell, &original_image, 0);
			int orig_img_idx = proc_cell.idx[0]*proc_cell.size[0];
			MPI_Recv(&original_image.data[orig_img_idx] , proc_cell.size_ , MPI_UINT8_T, p, collecttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			printf("Master has received from process %d.\n", p);
		}
		
		
		//write the file header
		if (read_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}

		//write the image data
		if (read_pgm_data( &original_image , outfile, &header_offs)== -1 ) {
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
		get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, 0);
		int local_img_idx = trim_halo_1D( proc_cell, local_image.pix_bytes, halowidth);
		
		MPI_Send(&local_image.data[local_img_idx], proc_cell.size_, MPI_UINT8_T, 0, collecttag, MPI_COMM_WORLD);
		
	}
		
	
	
	
	/*
	//write the file header
	if (proc_id == 0) {
		if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			clear_pgm( &local_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
	}
	//re-send out the header offset
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	
	//update the cells properties removing the halos since we don't write them back
	get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, 0);
	
	printf("From process %d, the image cell is %d x %d and located at (%d,%d) in the image.\n",
		   proc_id, proc_cell.size[0], proc_cell.size[1] , proc_cell.idx[0], proc_cell.idx[1] );
	
	MPI_Datatype mpi_cell_nohalo;
	MPI_Type_create_subarray (
		2, 
		original_image.size,
		proc_cell.size, 
		proc_cell.idx, 
		MPI_ORDER_C,
		MPI_UINT8_T,
		&mpi_cell_nohalo
	);
	MPI_Type_commit(&mpi_cell_nohalo);
	
	MPI_File out_file;
	MPI_File_open(MPI_COMM_CART, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY ,MPI_INFO_NULL, &out_file);
	
	//set the view past the header
	MPI_File_set_view(out_file, header_offs ,MPI_UINT8_T, mpi_cell_nohalo, "native", MPI_INFO_NULL);
	
	
	//write into local buffer
	int local_img_idx = trim_halo_1D( proc_cell, local_image.pix_bytes, halowidth);
	//printf("\nProcess %d is writing its image data.\n the buffer is %d elements and starts at  ", proc_id, proc_cell.size_ );
	MPI_File_write(out_file, &local_image.data[local_img_idx] , proc_cell.size_, MPI_UINT8_T, &status);
		
	MPI_File_close(&in_file);
	MPI_Type_free(&mpi_cell_nohalo);

	*/

	clear_pgm( &original_image);
	//the local image points to the memory location of the original image which is free
	local_image.data = NULL;
	delete_kernel(&kernel_ptr);
	clear_pgm( &local_image);
	MPI_Finalize();
	return 0;
	

} 