#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <blur_pgm.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>



int main( int argc, char **argv ) 
{ 
	#ifdef TIME
	clock_t blur_t, header_time, read_time, write_time, total_t;
	total_t = clock();
    #endif 
	
	//read command-line arguments and initialise the variables
	
	
	char infile[80] = "";
	char outfile[80] = "";
	
	pgm  image = new_pgm();
	kernel_t kernel;
	long int header_offs=0;
	

	//read command line parameters 
	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		return -1;
	}
	
	
	/*
	printf("ker size : %d x %d \n", kernel.size[0],kernel.size[1] );
	printf("ker hsize : %d x %d \n", kernel.halfsize[0],kernel.halfsize[1] );
	double norm=0;
	for (size_t i=0; i<kernel.size[0]*kernel.size[1]; ++i) {
		norm += kernel.ker[i];
		printf("ker %ld %f %f\n",i,kernel.ker[i], kernel.kernorm[i]);	
	}
	printf ("ker norm %f\n.",norm);
	return 0;
	*/
	
	#ifdef TIME
	header_time = clock();
	#endif
	//read the file header
	if (read_pgm_header( &image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
	header_time = clock() - header_time;
	#endif
	
	//allocate memory for the image
	if (allocate_pgm_memory( &image)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel);
		return -1;
	}
	
	#ifdef TIME
	read_time = clock();
	#endif
	//read the image data
	if (read_pgm_data( &image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
		read_time = clock() - read_time;
		#endif
	endian_swap(&image);
	
    #ifdef INFO
    	printf("Input file \"%s\" has been read.\n",infile);
		printf("The image is %d x %d.\n",image.size[0],image.size[1]);
	#endif
	
	#ifdef TIME
	blur_t = clock();
	#endif
	//pgm_blur_copy( &image, &kernel );
	//pgm_blur_linebuf( &image, &kernel );
	
	blur_func_manager( &image, &kernel );
	
	#ifdef TIME
	blur_t = clock() - blur_t;
	#endif
	
	#ifdef TIME
	header_time = clock();
	#endif
    //write the file header
	if (write_pgm_header( &image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef INFO
	printf("The header is %li bytes.\n",header_offs);
	#endif
	
	#ifdef TIME
	header_time = clock() - header_time;
	write_time = clock();
	#endif
	
	endian_swap(&image);
	
	//write the image data
	if (write_pgm_data( &image , outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel);
		return -1;
	}
	#ifdef TIME
	write_time = clock() - write_time;
	#endif
		
	#ifdef INFO
    printf("Output file \"%s\" has been written.\n",outfile);
	#endif
	
	#ifdef TIME
	total_t = clock()-total_t;
	printf("Header read time 	: %f s.\n",(double)(header_time)/ CLOCKS_PER_SEC );
	printf("Data read time		: %f s.\n",(double)(read_time)/ CLOCKS_PER_SEC );
	printf("Blur time   		: %f s.\n", (double)(blur_t) / CLOCKS_PER_SEC );
	printf("Header write time 	: %f s.\n",(double)(header_time)/ CLOCKS_PER_SEC );
	printf("Data write time 	: %f s.\n",(double)(write_time)/ CLOCKS_PER_SEC );
	printf("Total time 			: %f s.\n", (double)(total_t) / CLOCKS_PER_SEC );
	#endif

	clear_pgm( &image);
	delete_kernel(&kernel);
    return 0;
} 