#include <pgm.h>

#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include <math.h>
#include <time.h>


void generate_image(pgm* image, int type) {
	
	register double pixel;
	register size_t w = image->size[0];
	register size_t h = image->size[1];

	
	if( type == 0) {
		
		double scale = ((double)image->maxval )/((double) RAND_MAX + 1 ) ;
			
		srand(time(NULL));	
		
		uint8_t* img = (uint8_t*)image->data; 
		
		if (image->maxval <= 255 ) {
			for ( size_t i = 0; i< h; ++i) {
				for (size_t j=0; j< w; ++j) {
					pixel = ((double)rand() )*scale;
					img[i*w + j] = (uint8_t)pixel;

				}
				
			}

		} else {
			uint16_t* img = (uint16_t*)image->data; 
			for ( size_t i = 0; i< h; ++i) {
				for (size_t j=0; j< w; ++j) {
					pixel = ((double)rand() )*scale;
					img[i*w + j] = (uint16_t)pixel;

				}
				
			}
		}
	}
	
	
}


int main( int argc, char **argv ) {
	
	int xsize,ysize,type,maxval;
	char outfile[80] = "gen_image.pgm";
	long int header_offs=0;
	
	if ( argc< 4 ) {
		printf("usage: %s image-type max-brightness x-size y-size [out_name.pgm]\n",argv[0]);
		return 0;
	} else {
		type=atoi(argv[1]);
		maxval=atoi(argv[2]);
		xsize=atoi(argv[3]);
		ysize=atoi(argv[4]);
		
		
		if (argc > 5 ) {
			strcpy(outfile,argv[5]);
		}
	}
	
		
	pgm image = new_pgm();
	image.size[0]=xsize;
	image.size[1]=ysize;
	image.maxval=maxval;
	image.pix_bytes =  1 + (image.maxval > 255);
	
	//allocate memory for the image
	if (allocate_pgm_memory( &image)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		return -1;
	}

	
	
	generate_image(&image, type);
	
	printf("Done Generating\n");
	

	
	//write the file header
	if (write_pgm_header( &image , outfile, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		return -1;
	}

	//write the image data
	if (write_pgm_data( &image , outfile)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image);
		return -1;
	}
	
	 printf("Output file \"%s\" has been written.\n",outfile);

	clear_pgm( &image);
	return 0;
}