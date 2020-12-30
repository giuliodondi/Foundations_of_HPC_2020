#include <pgm.h>

#include <stdlib.h>
#include <stdio.h> 

#define XWIDTH 512
#define YWIDTH 512
#define MAXVAL 65534

#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif

char  generate_gradient( pgm* gradient );


int main( int argc, char **argv ) 
{ 
    int xsize      = XWIDTH;
    int ysize      = YWIDTH;
    int maxval     = MAXVAL;
    
    long int header_offs=0;

 
    if ( argc > 2 ) {
        xsize = atoi( argv[1] );
        ysize = atoi( argv[2] );
       
        if ( argc > 3 ) {
             maxval = atoi( argv[3] ) % 65536;
        }
    }
    
   
    
    pgm gradient = new_pgm();
    gradient.size[0] = xsize;
    gradient.size[1] = ysize;
    gradient.maxval = maxval;
    gradient.pix_bytes = 1 + (maxval > 255);
    
     printf("%d %d %d %d \n", xsize, ysize, maxval, gradient.pix_bytes);

    // ---------------------------------------------------
    //
    
    if (generate_gradient( &gradient )== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &gradient);
		return -1;
	}
    printf("The gradient has been generated\n");

    //write the file header
	if (write_pgm_header( &gradient , "gradient.pgm", &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &gradient);
		return -1;
	}
	printf("the header is %li bytes.\n",header_offs);
	
	//write the image data
	if (write_pgm_data( &gradient , "gradient.pgm")== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &gradient);
		return -1;
	}
    
   clear_pgm( &gradient);
    return 0;
} 


char  generate_gradient( pgm* gradient )
/*
 * just and example about how to generate a vertical gradient
 * maxval is either 255 or 65536, xsize and ysize are the
 * x and y dimensions of the image to be generated.
 * The memory region that will contain the image is returned
 * by the function as a void *

 */
{
    
    int xsize      = gradient->size[0];
    int ysize      = gradient->size[1];
    int maxval     = gradient->maxval;
    int minval      = 0; 
    double delta    = (maxval - minval)/ ( (double)ysize );
    
    
     if(delta < 1 )
    delta = 1;
    

    if (allocate_pgm_memory( gradient) == -1) {
        return -1;
    }
    
   

    if( maxval < 256 )
    // generate a gradient with 1 byte of color depth
    {
        uint8_t* cImage = (uint8_t*)gradient->data; 
        uint8_t _maxval = (uint8_t)maxval;
        int idx = 0;
        for ( int yy = 0; yy < ysize; yy++ )
        {
            uint8_t value = minval + (int)(yy*delta);
            for( int xx = 0; xx < xsize; xx++ )
                cImage[idx++] = (value > _maxval)?_maxval:value;
        }
    }
    else
    // generate a gradient with 2 bytes of color depth
    {
        uint16_t* sImage = (uint16_t*)gradient->data; 
        uint16_t _maxval = (unsigned short int)maxval;
        int idx = 0;
        for ( int yy = 0; yy < ysize; yy++ )
        {
            uint16_t value  = (uint16_t) (minval+ yy*delta);

            for( int xx = 0; xx < xsize; xx++ )
                sImage[idx++] = (value > maxval)?_maxval:value;
        }
    }

    return 0;
}

