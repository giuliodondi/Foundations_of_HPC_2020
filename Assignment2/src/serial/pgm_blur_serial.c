#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>





void pgm_blur_copy(  pgm* input_img , const kernel_t* k) {
	
	
	register const int xdim = input_img->width ;
	register const int ydim = input_img->height ;
	register const int xdim1 = input_img->width - 1;
	register const int ydim1 = input_img->height  - 1;

	
	register const int ker_s = k->size ;
	register const int ker_hsize = ( k->size - 1)/2 ;
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	
	uint16_t* image = (uint16_t*)input_img->data; 
	
	//store the blurred pixels in a buffer
	//it will essentially be a new image that will rpelace the input one at the end
	uint16_t* out_image = (uint16_t*)malloc( xdim*ydim*sizeof(uint16_t) );
	

	
	register int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	register uint8_t normflagx,normflagy;

	for (int i=0; i<ydim; ++i) {	
		
		normflagy = ( (i<ker_hsize) || (i>=(ydim - ker_hsize) ) );
		
		for (int j=0; j<xdim; ++j) {
			
			normflagx = normflagy +  ( (j<ker_hsize) || (j>=(xdim - ker_hsize) ) );
			
			/*
			general code to stay within the boundaries of the image
			the offsets # are the pixels between the current pixel and the edge
			they provide the boundaries of the kernel scanning loop
			image pixels and kernel indices are accordingly adjusted
			
			ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
			
			overlaying the kernel over the top left corner of the image 
			the offsets become = l 0, u 0, r +1, d +1
			*/
			
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
			//vignette renormalisation
			//we use the normalisation matrix pre-computed and use the offsets variables
			//to extract the right entry, which is the normalisation value
			//if the kernel is fully within the borders the indices would point to the central value
			//which is always 1, we skip these calculations
			
			if (normflagx ) {
				normc=ker_norm[ker_s*(ker_hsize + offs_u + offs_d) +  ker_hsize + offs_l + offs_r];
				accum = accum*normc;
			}
			out_image[i*xdim + j] = (uint16_t)accum;
		}
	}
	
	input_img->data = (uint8_t*)out_image;
	free(image);
	printf("Blurring complete.\n");
	return ;

	
}





void pgm_blur_linebuf(  pgm* input_img , const kernel_t* k) {
	
	
	register const int xdim = input_img->width ;
	register const int ydim = input_img->height ;
	register const int xdim1 = input_img->width - 1;
	register const int ydim1 = input_img->height  - 1;

	
	register const int ker_s = k->size ;
	register const int ker_hsize = ( k->size - 1)/2 ;
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	
	uint16_t* image = (uint16_t*)input_img->data; 
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
	register const int buffer_lines = (ker_hsize + 1);
	register int line_idx; 
	int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
	
	register int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	register uint8_t normflagx,normflagy;

	
	
	for (int i=0; i<ydim; ++i) {
		line_idx = xdim*(i%buffer_lines);
		
		if (i>=buffer_lines) {
			memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
		}
		
		//are we close to the left/right edges?
		normflagy = ( (i<ker_hsize) || (i>=(ydim - ker_hsize) ) );
		
		for (int j=0; j<xdim; ++j) {
			
			
			//are we close to the top/bottom edges?
			normflagx = normflagy + ( (j<ker_hsize) || (j>=(xdim - ker_hsize) ) );
			
			/*
			general code to stay within the boundaries of the image
			the offsets # are the pixels between the current pixel and the edge
			they provide the boundaries of the kernel scanning loop
			image pixels and kernel indices are accordingly adjusted
			
			ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
			
			overlaying the kernel over the top left corner of the image 
			the offsets become = l 0, u 0, r +1, d +1
			*/
			
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
			
			//vignette renormalisation
			//we use the normalisation matrix pre-computed and use the offsets variables
			//to extract the right entry, which is the normalisation value
			//if the kernel is fully within the borders the indices would point to the central value
			//which is always 1, we skip these calculations
			
			if (normflagx ) {
				printf("i j f: %d %d %d\n",i,j,normflagx);
				normc=ker_norm[ker_s*(ker_hsize + offs_u + offs_d) +  ker_hsize + offs_l + offs_r];
				accum = accum*normc;
			}
			
			linebuf[line_idx + j] = (uint16_t)accum;
		}
	}
	
	
	//the buffer is full of lines still to write
	//but the next buffer line to write will be at an arbitrary index
	
	//this is the buffer line that contains the next line of data
	int j = (ydim - buffer_lines)%buffer_lines ;
	
	//write the data in the last lines of the image
	//when we reach the end of the buffer, wrap around as the last lines to write are at the top
	for (int i = ydim - buffer_lines; i<ydim; ++i) {
		memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
		j = (j + 1)%buffer_lines;
	}
	
	

	free(linebuf);
	printf("Blurring complete.\n");
	return ;

	
}




void pgm_blur_linebuf_unrol(  pgm* input_img , const kernel_t* k) {
	
	
	register const int xdim = input_img->width ;
	register const int ydim = input_img->height ;
	register const int xdim1 = input_img->width - 1;
	register const int ydim1 = input_img->height  - 1;

	
	register const int ker_s = k->size ;
	register const int ker_hsize = ( k->size - 1)/2 ;
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	
	uint16_t* image = (uint16_t*)input_img->data; 
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
	register const int buffer_lines = (ker_hsize + 1);
	register int line_idx; 
	int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
	
	register int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	register uint8_t normflagx,normflagy;

	
	
	for (int i=0; i<ydim; ++i) {
		line_idx = xdim*(i%buffer_lines);
		
		if (i>=buffer_lines) {
			memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
		}
		
		//are we close to the left/right edges?
		normflagy = ( (i<ker_hsize) || (i>=(ydim - ker_hsize) ) );
		
		for (int j=0; j<xdim; ++j) {
			
			
			//are we close to the top/bottom edges?
			normflagx = normflagy + ( (j<ker_hsize) || (j>=(xdim - ker_hsize) ) );
			
			/*
			general code to stay within the boundaries of the image
			the offsets # are the pixels between the current pixel and the edge
			they provide the boundaries of the kernel scanning loop
			image pixels and kernel indices are accordingly adjusted
			
			ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
			
			overlaying the kernel over the top left corner of the image 
			the offsets become = l 0, u 0, r +1, d +1
			*/
			
			offs_l = -min( ker_hsize, j );
			offs_u = -min( ker_hsize, i );
			offs_r = min( ker_hsize, xdim1 - j );
			offs_d = min( ker_hsize, ydim1 - i );
			
			
			accum=0;
			if (( offs_r - offs_l)%2) {
				//this branch selects the odd-sized convolutions which can only happen if the normflag is enabled
				//so this branch is never entered when we're far from the borders
				for (int k = offs_u; k<= offs_d; ++k) {
					accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + offs_l  ]*( image[(i + k)*xdim + j + offs_l] )
							+ kernel[ ker_s*(ker_hsize + k) +  ker_hsize + offs_l + 1  ]*( image[(i + k)*xdim + j + offs_l + 1] );
					for (int t =  offs_l + 2; t<= offs_r; t+=2) {
						//printf("k t : %d %d \n",k,t);
						//getchar();
						accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ]*( image[(i + k)*xdim + j + t] ) 
								+ kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t +1 ]*( image[(i + k)*xdim + j + t + 1] );
					}

				}
			}
			else {
				for (int k = offs_u; k<= offs_d; ++k) {
					accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + offs_l  ]*( image[(i + k)*xdim + j + offs_l] );
					for (int t =  offs_l + 1; t<= offs_r; t+=2) {
						//printf("k t : %d %d \n",k,t);
						//getchar();
						accum += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ]*( image[(i + k)*xdim + j + t] ) 
								+ kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t +1 ]*( image[(i + k)*xdim + j + t + 1] );
					}

				}	
			}
			
			//vignette renormalisation
			//we use the normalisation matrix pre-computed and use the offsets variables
			//to extract the right entry, which is the normalisation value
			//if the kernel is fully within the borders the indices would point to the central value
			//which is always 1, we skip these calculations
			
			if (normflagx ) {
				printf("i j f: %d %d %d\n",i,j,normflagx);
				normc=ker_norm[ker_s*(ker_hsize + offs_u + offs_d) +  ker_hsize + offs_l + offs_r];
				accum = accum*normc;
			}
			
			linebuf[line_idx + j] = (uint16_t)accum;
		}
	}
	
	
	//the buffer is full of lines still to write
	//but the next buffer line to write will be at an arbitrary index
	
	//this is the buffer line that contains the next line of data
	int j = (ydim - buffer_lines)%buffer_lines ;
	
	//write the data in the last lines of the image
	//when we reach the end of the buffer, wrap around as the last lines to write are at the top
	for (int i = ydim - buffer_lines; i<ydim; ++i) {
		memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
		j = (j + 1)%buffer_lines;
	}
	
	

	free(linebuf);
	printf("Blurring complete.\n");
	return ;

	
}


