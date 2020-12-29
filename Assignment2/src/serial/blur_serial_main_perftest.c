#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


#define NITER 10

//serial function headers
void pgm_blur_copy(  pgm* input_img , kernel_t* k);
void pgm_blur_linebuf(  pgm* input_img , kernel_t* k);
void pgm_blur_linebuf_unrol(  pgm* input_img , kernel_t* k);


int main( int argc, char **argv ) 
{ 

	

    
	//read command-line arguments and initialise the variables
	
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";
	
	pgm  original_image = new_pgm();
	kernel_t kernel;
	long int header_offs=0;	

	
	
	//read command line parameters 
	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	

	//read the file header
	if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	
	//allocate memory for the image
	if (allocate_pgm_memory( &original_image)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	
	//read the image data
	if (read_pgm_data( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	
	
    
    printf("Input file \"%s\" has been read.\n",infile);
	printf("The image is %d x %d.\n",original_image.size[0],original_image.size[1]);

	
	
	pgm  image = new_pgm();
	
	
	printf("running the copy algorithm %d times.\n", NITER);
	
	double avg_time = 0;
	for (int n = 0; n< NITER; ++n) {
		copy_pgm( &original_image, &image) ;
		
		clock_t begin = clock();
   		pgm_blur_copy( &image, &kernel );
		clock_t end = clock();
		avg_time += (double)(end - begin) / CLOCKS_PER_SEC ;
		
		
	}
	
	printf("Average runtime : %f s \n", avg_time/NITER );
	
	


	
	
	

	printf("running the line-buffered algorithm %d times.\n", NITER);
	
	avg_time = 0;
	for (int n = 0; n< NITER; ++n) {
		copy_pgm( &original_image, &image) ;
		clock_t begin = clock();
   		pgm_blur_linebuf( &image, &kernel );
		clock_t end = clock();
		avg_time += (double)(end - begin) / CLOCKS_PER_SEC ;
	}
	
	printf("Average runtime : %f s \n", avg_time/NITER );
	
	

	printf("running the line-buffered + unrolling algorithm %d times.\n", NITER);
	
	avg_time = 0;
	for (int n = 0; n< NITER; ++n) {
		copy_pgm( &original_image, &image) ;
		clock_t begin = clock();
   		pgm_blur_linebuf_unrol( &image, &kernel );
		clock_t end = clock();
		avg_time += (double)(end - begin) / CLOCKS_PER_SEC ;
	}
	
	printf("Average runtime : %f s \n", avg_time/NITER );
	
	

        
   
	//write the file header
	if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
	printf("the header is %li bytes.\n",header_offs);
	
	//write the image data
	if (write_pgm_data( &original_image , outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel);
		return -1;
	}
    printf("Output file \"%s\" has been written.\n",outfile);
	
	

	clear_pgm( &image);
	delete_kernel(&kernel);
    return 0;
} 