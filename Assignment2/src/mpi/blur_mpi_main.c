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
void get_cell_1D(const int nprocs, const int proc_id, img_cell* proc_cell, const pgm* image, const unsigned int* halowidth);
int trim_halo_1D( const img_cell* proc_cell, const char img_bytes, const unsigned int* halowidth );
void pgm_blur_halo(  pgm* input_img , kernel_t* k,  const char* halos);


int main( int argc, char **argv ) 
{ 
	//mpi common variables
	int nprocs, proc_id;
	MPI_Status status;
	MPI_Offset my_data_offs;

	//initialise the MPI communicators
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	
	
	//variables common to all processes;
	double header_time, read_time, blur_time, write_time;
	pgm  local_image = new_pgm();
	pgm  original_image = new_pgm();
	kernel_t kernel;
	img_cell proc_cell;
	long int header_offs=0;
	unsigned int halowidth0[2] = {0,0};
	
	
	//read command-line arguments and initialise the variables
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";

	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		delete_kernel(&kernel);
		return -1;
		MPI_Finalize();
	}
	
	
	//only master reads the file header
	if (proc_id == 0) {
		header_time = MPI_Wtime();
		//read the file header
		if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			clear_pgm( &local_image);
			delete_kernel(&kernel);
			MPI_Finalize();
			return -1;
		}
		header_time = MPI_Wtime() - header_time;
		printf("Header read time for process %d : %f s.\n",proc_id,header_time);
	}
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	
	//send info on the incoming image metadata
	MPI_Bcast(&original_image.size, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.maxval, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&original_image.pix_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	
	//build the cells accounting for the halo above and below the cells
	get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, kernel.halfsize);
	
	printf("Process %d will process %d rows, %d cols of the image.\n", proc_id, proc_cell.size[1], proc_cell.size[0] );
	
	//initialise the local image
	local_image.size[0] = proc_cell.size[0];
	local_image.size[1] = proc_cell.size[1];
	local_image.maxval = original_image.maxval;
	local_image.pix_bytes = original_image.pix_bytes;
	//allocte data for the local image
	local_image.data = (uint8_t*)malloc( proc_cell.size_*sizeof(uint8_t) );
	
	
	//parallel file read
	read_time = MPI_Wtime();
	//wortk out the beginning of the buffer in the original image
	my_data_offs = header_offs + ( original_image.size[0]*proc_cell.idx[1]+ proc_cell.idx[0])*original_image.pix_bytes;
	
	MPI_File in_file;
	MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY,MPI_INFO_NULL, &in_file);
	MPI_File_read_at(in_file,  my_data_offs, &local_image.data[0]  , proc_cell.size_ , MPI_UINT8_T, &status);
	MPI_File_close(&in_file);
	
	endian_swap(&local_image);
	
	read_time = MPI_Wtime() - read_time;
	
	printf("Data read time for process %d : %f s.\n",proc_id,read_time);
	
	
	//blurring
	blur_time = MPI_Wtime();
	pgm_blur_halo( &local_image, &kernel , proc_cell.halos);
	blur_time = MPI_Wtime() - blur_time;

	printf("Blur time for process %d : %f s.\n",proc_id,blur_time);
	
	
	if (proc_id == 0) {
		header_time = MPI_Wtime();
		//write the file header
		if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			clear_pgm( &local_image);
			delete_kernel(&kernel);
			MPI_Finalize();
			return -1;
		}
		header_time = MPI_Wtime() - header_time;
		printf("Header write time for process %d : %f s.\n",proc_id,header_time);
	}
	MPI_Bcast(&header_offs, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
	
	//update cells setting halo width to zero
	get_cell_1D( nprocs, proc_id, &proc_cell, &original_image, halowidth0);
	my_data_offs = header_offs + ( original_image.size[0]*proc_cell.idx[1]+ proc_cell.idx[0])*original_image.pix_bytes;
	
	//parallel write
	write_time = MPI_Wtime();
	endian_swap(&local_image);
	
	MPI_File out_file;
	MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &out_file);

	MPI_File_write_at(out_file,  my_data_offs, &local_image.data[trim_halo_1D( &proc_cell, local_image.pix_bytes, kernel.halfsize)]  , proc_cell.size_ , MPI_UINT8_T, &status);
	
	MPI_File_close(&out_file);
	
	write_time = MPI_Wtime() - write_time;
	
	printf("Data write time for process %d : %f s.\n",proc_id,write_time);
	
	delete_kernel(&kernel);
	clear_pgm( &original_image);
	clear_pgm( &local_image);
	MPI_Finalize();
	return 0;
	

} 


