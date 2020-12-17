#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <kernel.h>


#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535

#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif


int main( int argc, char **argv ) 
{ 
	
	 int xsize      = 0;
    int ysize      = 0;
    int maxval     = 0;

    // print information about endianism
    //printf("this machine is %s\n", (I_M_LITTLE_ENDIAN)?"little endian":"big endian");

    // you can use also the system-defined macro LITTLE_ENDIAN
    //printf("2nd check: this machine definitely is %s\n", (LITTLE_ENDIAN)?"little endian":"big endian");
    
	//read command-line arguments and initialise the variables
	char infile[80] = "";
	char outfile[80] = "";
	
	void* image_ptr ;
	kernel kernel_ptr;
	
	
	if (read_params_initialise_kernel(argc, argv, &infile, &outfile, &kernel_ptr) == -1 ) {
		printf("Aborting.\n");
		return -1;
	}
																			  
	
	printf("kernel size : %d\n", kernel_ptr.size);
	
	for (int i = 0; i < (kernel_ptr.size*kernel_ptr.size); ++i)
		printf("kernel[%d] = %f\n", i, kernel_ptr.ker[i]);

	
		
	

    read_pgm_image( &image_ptr, &maxval, &xsize, &ysize, infile);
    printf("The image has been read\n");

    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
      swap_image( image_ptr, xsize, ysize, maxval);

    // do something on the image (for instance, blur it)
    // ...
    //
    
    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
      swap_image( image_ptr, xsize, ysize, maxval);
        
    write_pgm_image( image_ptr, maxval, xsize, ysize, outfile);
    printf("The image has been written back\n");
	
	

    free(image_ptr);
	delete_kernel(&kernel_ptr);
    return 0;
} 