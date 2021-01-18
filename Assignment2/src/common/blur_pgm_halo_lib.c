#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>


/*
Haloed Blurring functions, take as input a struct containing the image data and a kernel struct 
containing the kernel and normalization matrices
the input halos array contains the width in pixels of the left, up, right, down halo layers

the halos define the boundaries of the outer loop as we only blur the internal region
	

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

only a small portion of the image is buffered
the buffer is as wide as the internal image portion and has half as many lines as the kernel
including the middle line

we roll through the buffer and fill the lines with the blurred pixels

only when the buffer is full we start back at line 0 and start memcpy its contents 
starting from line 0 of the original image
we're safe that we only overwrite when the data is no longer needed


the lines actualy procesed may not fill entirely the buffer as the kernel size and the image slice dimensions
are unrelated. so after the loop we either simply copy the first lines of the buffer as needed or 
if the buffer was filled at least once, we proceed as in the un-haloed case

we find the right buffer line to write next, memcpy line by line until we reach the end of the buffer
then wrap around the top and write the last lines  

The cases of an 8 or 16-bit encoded pgm are handled separately as the buffers need a different 
data type



*/




//this function implements some x2 unrolling of the convolution loop
//in the border regions the kernel sub-portion may not be odd-sized
//so the even case is handled separetely
void pgm_blur_halo_unrolx2(  pgm* input_img , const kernel_t* k,  const int* halos) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;
	

	double* kernel = k->ker;
	double* ker_norm = k->kernorm;

	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	size_t bound_l = halos[0];
	size_t bound_u = halos[1];
	size_t bound_r = xdim - halos[2];
	size_t bound_d = ydim - halos[3];
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	int linebuf_length = (bound_r - bound_l);

	
	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*linebuf_length*sizeof(int16_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);


			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int16_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				if (( offs_r - offs_l)%2) {
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
				linebuf[line_idx + j - bound_l] = (uint16_t)min(input_img->maxval,accum);
			}
		}
		//the buffer is full of lines still to write, 2 different cases
		if ((bound_d - bound_u)<=buffer_lines) {
			//not a single line of the buffer has been written since the "image lines" do not completely fill the buffer
			//therefore we write the buffer lines from the beginning until the last valid line
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				++j;
			}
		}
		else {
			//the buffer has overflowed so every line must be written
			//but the next buffer line to write depends where the loop ended

			//this is the buffer line that contains the next line of data
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			//write the data in the last lines of the image
			//when we reach the end of the buffer, wrap around as the last lines to write are at the top
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*linebuf_length*sizeof(int8_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);


			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int8_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				if (( offs_r - offs_l)%2) {
					for (int u = offs_u; u<= offs_d; ++u) {
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + offs_l  ]*( image[(i + u)*xdim + j + offs_l] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + offs_l + 1 ]*( image[(i + u)*xdim + j + offs_l + 1] );
						for (int t =  offs_l + 2; t<= offs_r; t+=2) {
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
				linebuf[line_idx + j - bound_l] = (int8_t)min(input_img->maxval,accum);;
			}
		}
		if ((bound_d - bound_u)<=buffer_lines) {
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				++j;
			}
		}
		else {
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
}



//this function implements x4 unrolling
//a single loop takes care of odd and even kernel sizes
void pgm_blur_halo_unrolx4(  pgm* input_img , const kernel_t* k,  const int* halos) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;
	

	double* kernel = k->ker;
	double* ker_norm = k->kernorm;

	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	size_t bound_l = halos[0];
	size_t bound_u = halos[1];
	size_t bound_r = xdim - halos[2];
	size_t bound_d = ydim - halos[3];
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	int linebuf_length = (bound_r - bound_l);

	int offs_l, offs_r, offs_u, offs_d;
	register double accum, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*linebuf_length*sizeof(int16_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);
			

			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int16_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (uint16_t)min(input_img->maxval,accum);;
			}
		}
		if ((bound_d - bound_u)<=buffer_lines) {
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				++j;
			}
		}
		else {
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*linebuf_length*sizeof(int8_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);


			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int8_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (int8_t)min(input_img->maxval,accum);;
			}
		}
		if ((bound_d - bound_u)<=buffer_lines) {
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				++j;
			}
		}
		else {
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
}


//this function implements x4 unrolling
//a single loop takes care of odd and even kernel sizes
void pgm_blur_halo_unrolx8(  pgm* input_img , const kernel_t* k,  const int* halos) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;
	

	double* kernel = k->ker;
	double* ker_norm = k->kernorm;

	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	size_t bound_l = halos[0];
	size_t bound_u = halos[1];
	size_t bound_r = xdim - halos[2];
	size_t bound_d = ydim - halos[3];
	
	
	const size_t buffer_lines = (hsize_v + 1);
	size_t line_idx; 
	int linebuf_length = (bound_r - bound_l);

	int offs_l, offs_r, offs_u, offs_d;
	register double accum1, accum2, normc;
	uint8_t normflagx,normflagy;

	if (input_img->pix_bytes == 2) {
		uint16_t* image = (uint16_t*)input_img->data; 
		int16_t* linebuf = (int16_t*)malloc( buffer_lines*linebuf_length*sizeof(int16_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);
			

			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int16_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (uint16_t)min(input_img->maxval,accum1);
			}
		}
		if ((bound_d - bound_u)<=buffer_lines) {
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				++j;
			}
		}
		else {
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
	else {
		uint8_t* image = (uint8_t*)input_img->data; 
		int8_t* linebuf = (int8_t*)malloc( buffer_lines*linebuf_length*sizeof(int8_t) );

		for (size_t i=bound_u; i<bound_d; ++i) {
			line_idx = linebuf_length*((i - bound_u)%buffer_lines);


			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int8_t));
			}

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (int8_t)min(input_img->maxval,accum1);
			}
		}
		if ((bound_d - bound_u)<=buffer_lines) {
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				++j;
			}
		}
		else {
			size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
			for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
}

//wrapper to switch between the unrolled versions for both omp and mpi
void  blur_halo_func_manager( pgm* input_img , const kernel_t* k, const int* halos) {
	
	//pgm_blur_halo_unrolx2(input_img,k,halos);
	pgm_blur_halo_unrolx4(input_img,k,halos);
	//pgm_blur_halo_unrolx8(input_img,k,halos);
	
}

