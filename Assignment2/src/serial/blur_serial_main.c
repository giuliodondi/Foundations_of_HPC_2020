#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

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
	kernel_t kernel_ptr;
	long int header_offs=0;
	

	//read command line parameters 
	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel_ptr) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	

	//read the file header
	if (read_pgm_header( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	
	//read the image data
	if (read_pgm_data( &original_image , infile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	
	
	
    
    printf("Input file \"%s\" has been read.\n",infile);
	printf("The image is %d x %d.\n",original_image.size[0],original_image.size[1]);

	
	
	clock_t begin = clock();
    

	//pgm_blur_copy( &original_image, &kernel_ptr );
	//pgm_blur_linebuf( &original_image, &kernel_ptr );
	pgm_blur_linebuf_unrol( &original_image, &kernel_ptr );
	
	clock_t end = clock();
	printf("Elapsed: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC );

	
    //write the file header
	if (write_pgm_header( &original_image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	printf("the header is %li bytes.\n",header_offs);
	
	//write the image data
	if (write_pgm_data( &original_image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &original_image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	
	
		
	
    printf("Output file \"%s\" has been written.\n",outfile);
	
	

	clear_pgm( &original_image);
	delete_kernel(&kernel_ptr);
    return 0;
} 