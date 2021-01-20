#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>


/*
Blurring functions, take as input a struct containing the image data and a kernel struct 
containing the kernel and normalization matrices

the hsize_h, hsize_v variables indicate the # of kernel elements to either side of the central element 
in the horizontal and vertical directions.
e.g. for a 5x3 kernel hsize_h=2 , hsize_v = 1


All functions perform a convolution pixel by pixel scanning the image data in row-major order
an inner loop performs the convolution scan for every pixel

The borders are handled by means of four offs_x indices which indicates the number of kernel elements to all
four sides of the central element which should be considered. 
if we're close to the border, these will be less than hsize_h and hsize_v since only a portion of the kernel
should be used.
e.g. for a 3x3 kernel overlaid on the top left corner of the image, these indices are :
 offs_l = 0, offs_u = 0, offs_r = +1, offs_d = +1

These indices determine the boundaries of the convolution loop and the image pixels to select.


the border vignette effect is handled by first figuring out when the border effect takes place using the 
normflagx , normflagy variables and then selecting a pre-computed appropriate normalisation constant

the ker_norm matrix stores the the normalization constants for each possible sub-portion of the kernel
when this matrix was created, the offs_x indices were used to encode the location of the right constant
within the ner_norm matrix and now are used to select the right value

All functions buffer at least a portion of the processed pixels before writing them back

In all functions the cases of an 8 or 16-bit encoded pgm are handled separately as the buffers need a different 
data type


*/



//this function buffers the entire image and copies it over in the end
void pgm_blur_copy(  pgm* input_img , const kernel_t* k) {
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;
	

	if (input_img->pix_bytes == 2) {
		
		uint16_t* image = (uint16_t*)input_img->data; 
	
		uint16_t* out_image = (uint16_t*)malloc( xdim*ydim*sizeof(uint16_t) );

		for (size_t i=0; i<ydim; ++i) {	

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				for (int u = offs_u; u<= offs_d; ++u) {
					for (int t =  offs_l; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] ) ;
					}

				}

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				out_image[i*xdim + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}
		input_img->data = (uint8_t*)out_image;
		free(image);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 

		uint8_t* out_image = (uint8_t*)malloc( xdim*ydim*sizeof(uint8_t) );

		for (size_t i=0; i<ydim; ++i) {	

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				for (int u = offs_u; u<= offs_d; ++u) {
					for (int t =  offs_l; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] ) ;
					}

				}
	
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				out_image[i*xdim + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		input_img->data = (uint8_t*)out_image;
		free(image);
	
	}

}



/*
this function buffers only a small portion of the image
the buffer is as wide as the image and has half as many lines as the kernel
including the middle line

we roll through the buffer and fill the lines with the blurred pixels

only when the buffer is full we start back at line 0 and start memcpy its contents 
starting from line 0 of the original image
we're safe that we only overwrite when the data is no longer needed

at the end of the loop all the calculations are done but the buffer is full of data still to be written
but the next buffer line to write will be at an arbitrary index that depends on which buffer line 
was last written to the image in the loop

we find the right buffer line to write next, memcpy line by line until we reach the end of the buffer
then wrap around the top and write the last lines  

*/

void pgm_blur_linebuf(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				for (int u = offs_u; u<= offs_d; ++u) {
					for (int t =  offs_l; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] ) ;
					}

				}
				
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}

		size_t j = (ydim - buffer_lines)%buffer_lines ;

		
		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*xdim*sizeof(int8_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int8_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);



				accum=0;
				for (int u = offs_u; u<= offs_d; ++u) {
					for (int t =  offs_l; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] ) ;
					}

				}
			
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}


//this function implements some x2 unrolling of the convolution loop
//in the border regions the kernel sub-portion may not be odd-sized
//so the even case is handled separetely
void pgm_blur_linebuf_unrolx2(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				if (( offs_r - offs_l)%2) {
					//this branch selects the even-sized convolutions
					//this branch is never entered when we're far from the borders
					for (int u = offs_u; u<= offs_d; ++u) {
						for (int t =  offs_l; t<= offs_r; t+=2) {
							accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
									+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] ) ;
						}

					}
				}
				else {
					for (int u = offs_u; u<= offs_d; ++u) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + offs_l  ]*( image[(i + u)*xdim + j + offs_l] ) ;
						for (int t =  offs_l + 1; t<= offs_r; t+=2) {
							accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
									+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] ) ;
						}

					}
				}
				
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}

		size_t j = (ydim - buffer_lines)%buffer_lines ;

		
		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*xdim*sizeof(int8_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int8_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				if (( offs_r - offs_l)%2) {
					//this branch selects the even-sized convolutions
					//this branch is never entered when we're far from the borders
					for (int u = offs_u; u<= offs_d; ++u) {
						for (int t =  offs_l; t<= offs_r; t+=2) {
							accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
									+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] ) ;
						}

					}
				}
				else {
					for (int u = offs_u; u<= offs_d; ++u) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + offs_l  ]*( image[(i + u)*xdim + j + offs_l] ) ;
						for (int t =  offs_l + 1; t<= offs_r; t+=2) {
							accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
									+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] ) ;
						}

					}
				}
			
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}




	
	


//this function implements x4 unrolling
//a single loop takes care of odd and even kernel sizes
void pgm_blur_linebuf_unrolx4(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 3) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						t+=4;
					}
					for (; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*xdim*sizeof(int8_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int8_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 3) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						t+=4;
					}
					for (; t<= offs_r; ++t) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}
				
				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}




//this function implements x8 unrolling
//a single loop takes care of odd and even kernel sizes
void pgm_blur_linebuf_unrolx8(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum1, accum2, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*xdim*sizeof(int16_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int16_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum1=0;
				accum2=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 7) {
						accum1 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						
						accum2 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 4 ]*( image[(i + u)*xdim + j + t + 4] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 5 ]*( image[(i + u)*xdim + j + t + 5] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 6 ]*( image[(i + u)*xdim + j + t + 6] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 7 ]*( image[(i + u)*xdim + j + t + 7] ) ;
						t+=8;
					}
					for (; t<= offs_r; ++t) {
						accum1 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}
				
				accum1 += accum2;

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum1 = accum1*normc;
				}
				linebuf[line_idx + j] = (uint16_t)min(input_img->maxval,accum1);
			}
		}
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int16_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*xdim*sizeof(int8_t) );
		
		for (size_t i=0; i<ydim; ++i) {
			line_idx = xdim*(i%buffer_lines);

			if (i>=buffer_lines) {
				memcpy( &image[(i-buffer_lines)*xdim] , (&linebuf[line_idx]), xdim*sizeof(int8_t) );
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum1=0;
				accum2=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 7) {
						accum1 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						
						accum2 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 4 ]*( image[(i + u)*xdim + j + t + 4] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 5 ]*( image[(i + u)*xdim + j + t + 5] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 6 ]*( image[(i + u)*xdim + j + t + 6] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 7 ]*( image[(i + u)*xdim + j + t + 7] ) ;
						t+=8;
					}
					for (; t<= offs_r; ++t) {
						accum1 += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}
				
				accum1 += accum2;

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum1 = accum1*normc;
				}
				linebuf[line_idx + j] = (uint8_t)min(input_img->maxval,accum1);
			}
		}
		
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}




