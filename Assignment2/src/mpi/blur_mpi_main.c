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

	//mpi common variables
	int nprocs, irank;
	double start_t, elapsed;
  	int master = 0;
  	int sendtag = 123;
	int recvtag = 321;

	//first initialise the MPI communicators
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);
	
	if ( argc <=1) {
		fprintf (stderr , "\nToo few arguments. Usage : mpirun -np n %s -input input.pgm -ker \n", argv[0] ) ;
		MPI_Finalize() ;
		exit(-1) ;
	  }

	
	
	//variables common to all processes;
	kernel_t kernel_ptr;
	pgm  local_image;
	int halowidth,tmp;
	char img_bytes;
	
	if (irank == master) {
		
		//variables only for master 
		uint8_t* commbuf;
		size_t childbuf_size, masterbuf_size ;
		int masterlines, childlines, buf_idx;
		
	
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
		
		halowidth = (kernel_ptr.size - 1)/2;
		//now figure out the number of image lines to send to each process
		//master will take care of the top and bottom halo aeras plus an arbitrary number of lines
		//we'll try to give master about half as much work as the other processors
		//but master must take care at least of the top and bottom halo regions
		
		masterlines = max(2*halowidth,floor((float)original_image.height/(2*nprocs - 1)));
		childlines = floor(original_image.height - masterlines)/(nprocs - 1) ;
	
		while ( ( original_image.height - masterlines - childlines*(nprocs - 1) )>0 ) {
			++masterlines ;	
			childlines = floor(original_image.height - masterlines)/(nprocs - 1) ;
		}
		
		
		
		printf("\nMaster will process %d lines of the image.\n", masterlines + 2*halowidth);
		printf("Each children process will receive %d lines of the image.\n", childlines);
		
		//add the top and bottom halos
		childlines += 2*halowidth;
		//subtract them from the master lines so we keep them separate
		masterlines -= 2*halowidth;
		
		//send info on the incoming image data
		MPI_Bcast(&original_image.width, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&childlines, 1, MPI_INT, master, MPI_COMM_WORLD);
		MPI_Bcast(&original_image.maxval, 1, MPI_INT, master, MPI_COMM_WORLD);
		img_bytes = (1 + (original_image.maxval > 255));
		
		//allocate the comm buffer
		childbuf_size =  original_image.width * childlines*img_bytes;
		commbuf = (uint8_t*)malloc( childbuf_size*sizeof(uint8_t) );
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
		tmp = masterlines;
		for (int p=1;p < nprocs; ++p) {
			buf_idx = tmp*original_image.width*img_bytes;
			printf("Master is sending to process %d.\n", p);
			MPI_Send(&original_image.data[buf_idx], childbuf_size, MPI_UINT8_T, p, sendtag, MPI_COMM_WORLD);
			tmp += childlines - 2*halowidth;
		}
		

		
		//create the working pgm image for master
		//first it will store the top portion, then the bottom
		//we already know everything about where the data starts and ends
		local_image.width = original_image.width;
		local_image.height = masterlines + 2*halowidth;
		local_image.maxval = original_image.maxval;
		
		
		masterbuf_size =  local_image.width * local_image.height*img_bytes;
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
		local_image.height = 2*halowidth + 1;
		
		masterbuf_size =  local_image.width * local_image.height*img_bytes;
		
		tmp = original_image.width*original_image.height*img_bytes - masterbuf_size;
		memcpy( &local_image.data[0] , &original_image.data[ tmp ], masterbuf_size*sizeof(uint8_t));

		//blur this portion nd write it back
		halos[1] = 1 ;
		halos[3] = 0;
		start_t = MPI_Wtime();
		pgm_blur_halo( &local_image, &kernel_ptr , halos);
		elapsed += MPI_Wtime() - start_t;

		memcpy( &original_image.data[tmp] , &local_image.data[0], masterbuf_size*sizeof(uint8_t));

		
		printf("Wall time for Master : %f s.\n",elapsed);
		
		
		//receive the children portions of the image
		//the comm buffer is the one we used to send them out, but without
		//the top and bottom halos
		//decrease the size by the halo memory size
		childbuf_size -= 2*halowidth*original_image.width*img_bytes;
		
		//this is where the first data buffer starts, includes the halo
		tmp = masterlines + halowidth;
		for (int p=1;p < nprocs; ++p) {
			buf_idx = tmp*original_image.width*img_bytes;
			MPI_Recv(&original_image.data[buf_idx] , childbuf_size , MPI_UINT8_T, p, recvtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			printf("Master has received from process %d.\n", p);
			tmp += childlines - 2*halowidth;	
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
		img_bytes = (1 + (local_image.maxval > 255));

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
		
		//send back the data without the halos
		halowidth = (kernel_ptr.size - 1)/2;
		//this is the index right after the halo
		tmp = halowidth*local_image.width*img_bytes;
		//trim off the top and bottom halos from the total size
		img_size -= 2*tmp;
		MPI_Send(&local_image.data[tmp], img_size, MPI_UINT8_T, master, recvtag, MPI_COMM_WORLD);
		
	}
		
	
	

	
	delete_kernel(&kernel_ptr);
	clear_pgm( &local_image);
	MPI_Finalize();
	return 0;
		
	

} 