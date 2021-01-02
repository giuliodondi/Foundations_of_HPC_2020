#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>


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
	
	
	
	//code duplication in exchange for if checks in the loops
	if (input_img->pix_bytes == 2) {
		
		uint16_t* image = (uint16_t*)input_img->data; 
	
		//store the blurred pixels in a buffer
		//it will essentially be a new image that will rpelace the input one at the end

		uint16_t* out_image = (uint16_t*)malloc( xdim*ydim*sizeof(uint16_t) );

		for (size_t i=0; i<ydim; ++i) {	

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

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
	
		//store the blurred pixels in a buffer
		//it will essentially be a new image that will rpelace the input one at the end

		uint8_t* out_image = (uint8_t*)malloc( xdim*ydim*sizeof(uint8_t) );

		for (size_t i=0; i<ydim; ++i) {	

			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

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




void pgm_blur_linebuf(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}

void pgm_blur_linebuf_unrolx2(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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

				
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

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

				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}


void pgm_blur_linebuf_unrolx4(  pgm* input_img , const kernel_t* k) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;

	
	double* kernel = k->ker;
	double* ker_norm = k->kernorm;
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 3) {
						//printf("t %d\n",t);
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						t+=4;
					}
					//printf("aa\n");
					for (; t<= offs_r; ++t) {
						//printf("t %d\n",t);
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}
				
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint16_t) min(input_img->maxval,accum);
			}
		}
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=0; j<xdim; ++j) {

				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h) ) );

				/*
				general code to stay within the boundaries of the image
				the offsets # are the pixels between the current pixel and the edge
				they provide the boundaries of the kernel scanning loop
				image pixels and kernel indices are accordingly adjusted
				ex. for a full 3x3 kernel the offses are = l -1, u - 1, r +1, d +1
				overlaying the kernel over the top left corner of the image 
				the offsets become = l 0, u 0, r +1, d +1
				*/

				offs_l = -min( hsize_h, j );
				offs_u = -min( hsize_v, i );
				offs_r = min( hsize_h, xdim - j - 1);
				offs_d = min( hsize_v, ydim - i - 1);

				accum=0;
				int t;
				for (int u = offs_u; u<= offs_d; ++u) {
					t = offs_l;
					while (t <= offs_r - 3) {
						//printf("t %d\n",t);
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 1 ]*( image[(i + u)*xdim + j + t + 1] )
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 2 ]*( image[(i + u)*xdim + j + t + 2] ) 
								+ kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t + 3 ]*( image[(i + u)*xdim + j + t + 3] ) ;
						t+=4;
					}
					//printf("aa\n");
					for (; t<= offs_r; ++t) {
						//printf("t %d\n",t);
						accum += kernel[ k->size[0]*(hsize_v + u) +  hsize_h + t  ]*( image[(i + u)*xdim + j + t] );
					}
				}
				
				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j] = (uint8_t) min(input_img->maxval,accum);
			}
		}
		
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (ydim - buffer_lines)%buffer_lines ;

		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
		for (size_t i = ydim - buffer_lines; i<ydim; ++i) {
			memcpy( &image[i*xdim] , (&linebuf[xdim*j]), xdim*sizeof(int8_t) );	
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}

}



void  blur_func_manager( pgm* input_img , const kernel_t* k ) {
	if (k->halfsize[0]>=4) {
		pgm_blur_linebuf_unrolx4(input_img,k);
	}
	else if (k->halfsize[0]>=2) {
		pgm_blur_linebuf_unrolx2(input_img,k);
	}
	
}

