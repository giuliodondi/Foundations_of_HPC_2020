#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <blur_serial_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>





int main( int argc, char **argv ) 
{ 

	
    
	//read command-line arguments and initialise the variables
	
	
	char infile[80] = "";
	char outfile[80] = "output.pgm";
	
	pgm  image ;
	kernel_t kernel_ptr;
	

	
	if (read_params_initialise_kernel(argc, argv, infile, outfile, &kernel_ptr) == -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
	

	if (read_pgm( &image , infile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel_ptr);
		return -1;
	}

	
    
    printf("Input file \"%s\" has been read.\n",infile);
	printf("The image is %d x %d.\n",image.width,image.height);
	
	
	clock_t begin = clock();
    

	//pgm_blur_copy( &image, &kernel_ptr );
	//pgm_blur_linebuf( &image, &kernel_ptr );
	pgm_blur_linebuf_unrol( &image, &kernel_ptr );
	//pgm_blur_linebuf_unrol3( &image, &kernel_ptr );
	
	clock_t end = clock();
	printf("Elapsed: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC );

        
   
	if ( write_pgm( &image, outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		delete_kernel(&kernel_ptr);
		return -1;
	}
    printf("Output file \"%s\" has been written.\n",outfile);
	
	

	clear_pgm( &image);
	delete_kernel(&kernel_ptr);
    return 0;
} 