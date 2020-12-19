#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>





void pgm_blur_copy(  pgm* input_img , kernel_t k) {
	
	
	register const int xdim = input_img->width ;
	register const int ydim = input_img->height ;
	register const int xdim1 = input_img->width - 1;
	register const int ydim1 = input_img->height  - 1;

	
	register const int ker_s = k.size ;
	register const int ker_hsize = ( k.size - 1)/2 ;
	double* kernel = k.ker;
	
	
	uint16_t* image = (uint16_t*)input_img->data; 
	uint16_t* out_image = (uint16_t*)malloc( xdim*ydim*2*sizeof(uint8_t) );
	

	
	register int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	register uint8_t normx=0, normy=0;

	for (int i=0; i<ydim; ++i) {	
		
		normy = ( (i<ker_hsize) || (i>=(ydim - ker_hsize) ) );
		
		for (int j=0; j<xdim; ++j) {
			
			if (!normy)
				normx = ( (j<ker_hsize) || (j>=(xdim - ker_hsize) ) );
			
			//printf("i j : %d %d\n",i,j);
			//printf("normx normy : %d %d\n",normx , normy);
			//general code to stay within the boundaries
			//if we're closer to the edges than the kernel half size
			
			offs_l = -min( ker_hsize, j );
			offs_u = -min( ker_hsize, i );
			offs_r = min( ker_hsize, xdim1 - j );
			offs_d = min( ker_hsize, ydim1 - i );
			
			
			accum=0;
			for (int k = offs_u; k<= offs_d; ++k) {
				for (int t =  offs_l; t<= offs_r; ++t) {
					accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ]*( image[(i + k)*xdim + j + t] ) ;
				}
				
			}
			//vignette kernel renormalization loop
			if (normy || normx ) {
				normc=0;
				for (int k = offs_u; k<= offs_d; ++k) {
					for (int t =  offs_l; t<= offs_r; ++t) {
						normc += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ] ;
					}
				
				}
				accum = accum/normc;
			}
			out_image[i*xdim + j] = (uint16_t)accum;
		}
	}
	
	input_img->data = (uint8_t*)out_image;
	free(image);
	printf("Blurring complete.\n");
	return ;

	
}





void pgm_blur_linebuf(  pgm* input_img , kernel_t k) {
	
	
	register const int xdim = input_img->width ;
	register const int ydim = input_img->height ;
	register const int xdim1 = input_img->width - 1;
	register const int ydim1 = input_img->height  - 1;

	
	register const int ker_s = k.size ;
	register const int ker_hsize = ( k.size - 1)/2 ;
	double* kernel = k.ker;
	
	
	uint16_t* image = (uint16_t*)input_img->data; 

	
	register const int buffer_lines = (ker_hsize + 1);
	register int line_idx; 
	
	int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*(2)*sizeof(int8_t) );
	
	register int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	register uint8_t normx=0, normy=0;

	
	
	for (int i=0; i<ydim; ++i) {
		line_idx = xdim*(i%buffer_lines);
		
		if (i>=buffer_lines) {
			//printf("line written : %d\n",i-2);
			//printf("line_idx : %d\n",i%buffer_lines);
			memcpy( &image[(i-buffer_lines)*xdim] , (uint8_t*)(&linebuf[line_idx]), xdim*(2)*sizeof(int8_t) );
		}
		
		//are we close to the top/bottom edges?
		normy = ( (i<ker_hsize) || (i>=(ydim - ker_hsize) ) );
		
		for (int j=0; j<xdim; ++j) {
			
			//are we close to the left/right edges?
			normx = ( (j<ker_hsize) || (j>=(xdim - ker_hsize) ) );
			
			//general code to stay within the boundaries of the image
			//the offsets # are the pixels between the current pixel and the edge
			//they provide the boundaries of the kernel scanning loop
			//image pixels and kernel value indices are accordingly adjusted
			
			offs_l = -min( ker_hsize, j );
			offs_u = -min( ker_hsize, i );
			offs_r = min( ker_hsize, xdim1 - j );
			offs_d = min( ker_hsize, ydim1 - i );
			
			
			accum=0;
			for (int k = offs_u; k<= offs_d; ++k) {
				for (int t =  offs_l; t<= offs_r; ++t) {
					accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ]*( image[(i + k)*xdim + j + t] ) ;
				}
				
			}
			//vignette kernel renormalization loop
			
			if (normy || normx ) {
					normc=0;
					for (int k = offs_u; k<= offs_d; ++k) {
						for (int t =  offs_l; t<= offs_r; ++t) {
							normc += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ] ;
						}

					}
				accum = accum/normc;
			}
			linebuf[line_idx + j] = (uint16_t)accum;
		}
	}
	
	
	for (int i=0; i<buffer_lines; ++i) {
		//printf("line written : %d\n",ydim-buffer_lines + i);
		//printf("line_idx : %d\n",i);
		memcpy( &image[(ydim-buffer_lines + i)*xdim] , (uint8_t*)(&linebuf[i]), xdim*(2)*sizeof(int8_t) );	
	}

	free(linebuf);
	printf("Blurring complete.\n");
	return ;

	
}