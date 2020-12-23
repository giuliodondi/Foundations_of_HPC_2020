#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include <mpi.h>
#define USE MPI

//serial function headers
void pgm_blur_serial(  pgm* input_img , kernel_t* k);
void pgm_blur_halo(  pgm* input_img , kernel_t* k,  const char* halos);





int main( int argc, char **argv ) 
{ 

	int nprocs, irank;
	double start_t, elapsed;

	//first initialise the MPI communicators
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);
	
	if ( argc <=1) {
		fprintf (stderr , "\nToo few arguments. Usage : mpirun -np n %s -input input.pgm -ker \n", argv[0] ) ;
		MPI_Finalize() ;
		exit(-1) ;
	  }

	 // master process
  	int master = 0;
  	int sendtag = 123;
	int recvtag = 321;
	
	//variables common to all processes;
	kernel_t kernel_ptr;
	pgm  local_image;
	
	if (irank == master) {
	
		//read command-line arguments and initialise the variables
		pgm  original_image ;
		kernel_t kernel_ptr;
		
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


		printf("Master has read input file \"%s\".\n",infile);
		printf("The image is %d x %d.\n",original_image.width,original_image.height);
		
		
		//first broadcast the kernel
		MPI_Bcast(&kernel_ptr.size, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(kernel_ptr.ker, kernel_ptr.size*kernel_ptr.size, MPI_DOUBLE, master, MPI_COMM_WORLD);
		MPI_Bcast(kernel_ptr.kernorm, kernel_ptr.size*kernel_ptr.size, MPI_DOUBLE, master, MPI_COMM_WORLD);
		
		
		//now figure out the number of image lines to send to each process
		//master will take care of the top and bottom halo aeras plus an arbitrary number of lines
		//we'll try to give master about half as much work as the other processors
		
		int halosize = (kernel_ptr.size - 1)/2;
		int masterlines = 100;
		int childlines = ceil( (original_image.height - 2*halosize -  masterlines)/(nprocs - 1) );
		int tmp = original_image.height - childlines*(nprocs - 1);
		while (tmp < 2*(halosize) ) {
			++masterlines;
			childlines  = ceil( (original_image.height - 2*halosize - masterlines)/(nprocs - 1) );
			tmp = original_image.height - childlines*(nprocs - 1);
		}
		//add the top and bottom halos
		childlines += 2*halosize;
		
		
		printf("\nMaster will process %d lines of the image.\n", masterlines + 2*halosize);
		printf("\nEach children process will receive %d lines of the image.\n", childlines);
		
		//send info on the incoming image data
		MPI_Bcast(&original_image.width, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&childlines, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&original_image.maxval, 1, MPI_INT, master, MPI_COMM_WORLD);
		char img_bytes = (1 + (original_image.maxval > 255));
		
		//allocate the comm buffer
		unsigned int buf_size =  original_image.width * childlines*img_bytes;
		uint8_t* commbuf = (uint8_t*)malloc( buf_size*sizeof(uint8_t) );
		if ( ! commbuf) {
			printf("Master couldn't allocate memory for the comm buffer.\n");
			free(commbuf);
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
		
		//send data to the processes
		//startline starts at the beginning of the first buffer
		//each time it increments to the start of the next buffer
		//accounting for the overlapping halos
		int buf_startline = masterlines;
		//this contains
		int buf_idx;
		for (int p=1;p < nprocs; ++p) {
			buf_idx = buf_startline*original_image.width*img_bytes;
			memcpy( commbuf , &original_image.data[buf_idx], buf_size*sizeof(uint8_t));
			printf("Master is sending to process %d.\n", p);
			MPI_Send(commbuf, buf_size, MPI_UINT8_T, p, sendtag, MPI_COMM_WORLD);
			buf_startline += childlines - 2*halosize;
		}
		

		
		//create the working pgm image for master
		//first it will store the top portion, then the bottom
		//we already know everything about where the data starts and ends
		pgm  local_image;
		local_image.width = original_image.width;
		local_image.height = masterlines + 2*halosize;
		local_image.maxval = original_image.maxval;
		
		
		unsigned int masterbuf_size =  local_image.width * local_image.height*img_bytes;
		local_image.data = (uint8_t*)malloc( masterbuf_size*sizeof(uint8_t) );
		if ( ! local_image.data) {
			printf("Master couldn't allocate memory for the image.\n");
			clear_pgm(&local_image);
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
		memcpy( local_image.data , &original_image.data[0], masterbuf_size*sizeof(uint8_t));
		
		//blur this portion nd write it back
		
		
		
		char halos[] = {0,0,0,1};
		start_t = MPI_Wtime();
		pgm_blur_halo( &local_image, &kernel_ptr , halos);
		elapsed = MPI_Wtime() - start_t;
		memcpy( &original_image.data[0] , &local_image.data[0], masterbuf_size*sizeof(uint8_t));
		
		
		
		
		//now the bottom portion
		local_image.height = 2*halosize;
		free(local_image.data);
		masterbuf_size =  local_image.width * local_image.height*img_bytes;
		local_image.data = (uint8_t*)malloc( masterbuf_size*sizeof(uint8_t) );
		if ( ! local_image.data) {
			printf("Master couldn't allocate memory for the image.\n");
			clear_pgm(&local_image);
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
		
		
			int orig_image_size = original_image.width*original_image.height*img_bytes;
			//int wk_img_size = local_image.width*local_image.height;
			memcpy( &local_image.data[0] , &original_image.data[orig_image_size - masterbuf_size], masterbuf_size*sizeof(uint8_t));
		
			//blur this portion nd write it back
			 halos[1] = 1 ;
			 halos[3] = 0;
			start_t = MPI_Wtime();
			pgm_blur_halo( &local_image, &kernel_ptr , halos);
			elapsed += MPI_Wtime() - start_t;
			
			memcpy( &original_image.data[orig_image_size - masterbuf_size] , &local_image.data[0], masterbuf_size*sizeof(uint8_t));
		
		
		printf("Wall time for Master : %f s.\n",elapsed);
		
		
		//receive the children portions of the image
		//the comm buffer is the one we used to send them out, but without
		//the top and bottom halos
		buf_size -= 2*halosize*original_image.width*img_bytes;
		
		//this is where the first data buffer starts, includes the halo
		buf_startline = masterlines + halosize;
		//send data to the processes
		for (int p=1;p < nprocs; ++p) {
			buf_idx = buf_startline*original_image.width*img_bytes;
			MPI_Recv(commbuf, buf_size , MPI_UINT8_T, p, recvtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			printf("Master has received from process %d.\n", p);
			memcpy( &original_image.data[buf_idx] , commbuf, buf_size*sizeof(uint8_t));
			buf_startline += childlines - 2*halosize;	
		}
		free(commbuf);
		
		
		
		//write the image
		if ( write_pgm( &original_image, outfile)== -1 ) {
			printf("Aborting.\n");
			clear_pgm( &original_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
	   	printf("Master has written output file \"%s\".\n",outfile);
		
		
		clear_pgm( &original_image);	
	}
	else {

		//get the kernel from Master
		
		MPI_Bcast(&kernel_ptr.size, 1, MPI_INT, master, MPI_COMM_WORLD);
		kernel_ptr.ker = (double*)malloc( kernel_ptr.size*kernel_ptr.size *sizeof(double));
		kernel_ptr.kernorm = (double*)malloc( kernel_ptr.size*kernel_ptr.size *sizeof(double));
		if (!(kernel_ptr.ker && kernel_ptr.kernorm)) {
			printf("Process %d couldn't allocate memory for the kernel.\n",	irank);
			delete_kernel(&kernel_ptr);
			return -1;
		}
		MPI_Bcast(kernel_ptr.ker, kernel_ptr.size*kernel_ptr.size, MPI_DOUBLE, master, MPI_COMM_WORLD);
		MPI_Bcast(kernel_ptr.kernorm, kernel_ptr.size*kernel_ptr.size, MPI_DOUBLE, master, MPI_COMM_WORLD);
		
		printf("Process %d received the kernel from Master.\n",	irank);
		
		//prepare the working pgm image
		MPI_Bcast(&local_image.width, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&local_image.height, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&local_image.maxval, 1, MPI_INT, master, MPI_COMM_WORLD);
		char img_bytes = (1 + (local_image.maxval > 255));

		//allocate memory for the image
		unsigned int img_size =  local_image.width * local_image.height*img_bytes;
		local_image.data = (uint8_t*)malloc( img_size*sizeof(uint8_t) );
		if ( ! local_image.data) {
			printf("Process %d couldn't allocate memory for the image.\n",	irank);
			clear_pgm(&local_image);
			delete_kernel(&kernel_ptr);
			return -1;
		}
		
		
		//receive the data
		MPI_Recv(local_image.data, img_size , MPI_UINT8_T, master, sendtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Process %d received a %d x %d image.\n",irank,local_image.width,local_image.height);
		
		char halos[] = {0,1,0,1};

		start_t = MPI_Wtime();
		pgm_blur_halo( &local_image, &kernel_ptr , halos);
		elapsed = MPI_Wtime() - start_t;
		printf("Wall time for process %d : %f s.\n",irank,elapsed);
		
		//send back the data
		//we need to trim off the halos at the top and bottom
		int halosize = (kernel_ptr.size - 1)/2;
		img_size -= 2*halosize*local_image.width*img_bytes;
		MPI_Send(&local_image.data[halosize*local_image.width*img_bytes], img_size, MPI_UINT8_T, master, recvtag, MPI_COMM_WORLD);
		
	}
		
	
	

	
	delete_kernel(&kernel_ptr);
	clear_pgm( &local_image);
	MPI_Finalize();
	return 0;
		
	

} 