#include <pgm.h>
#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>




void pgm_blur_halo_unrolx2(  pgm* input_img , const kernel_t* k,  const int* halos) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;
	

	double* kernel = k->ker;
	double* ker_norm = k->kernorm;

	
	/* the halos array contains four elements corresponding to left, up, right, down edges
		the values are 0 or 1 which indicates whether the edge in question contains a halo
		(=1) or if it corresponds to the actual image edge (=0)
		if there is a halo we offset the correspnding loop boundary by the kernel halfsize
		
	*/
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	size_t bound_l = halos[0];
	size_t bound_u = halos[1];
	size_t bound_r = xdim - halos[2];
	size_t bound_d = ydim - halos[3];
	
	
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				//are we close to the top/bottom edges?
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (uint16_t)accum;
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				//are we close to the top/bottom edges?
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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

				//vignette renormalisation
				//we use the normalisation matrix pre-computed and use the offsets variables
				//to extract the right entry, which is the normalisation value
				//if the kernel is fully within the borders the indices would point to the central value
				//which is always 1, we skip these calculations

				if (normflagx ) {
					normc=ker_norm[k->size[0]*(hsize_v + offs_u + offs_d) +  hsize_h + offs_l + offs_r];
					accum = accum*normc;
				}
				linebuf[line_idx + j - bound_l] = (int8_t)accum;
			}
		}
		//the buffer is full of lines still to write, 2 different cases
		if ((bound_d - bound_u)<=buffer_lines) {
			//not a single line of the buffer has been written since the "image lines" do not completely fill the buffer
			//therefore we write the buffer lines from the beginning until the last valid line
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
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
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
}


void pgm_blur_halo_unrolx4(  pgm* input_img , const kernel_t* k,  const int* halos) {
	
	
	register const size_t xdim = input_img->size[0] ;
	register const size_t ydim = input_img->size[1] ;
	

	double* kernel = k->ker;
	double* ker_norm = k->kernorm;

	
	/* the halos array contains four elements corresponding to left, up, right, down edges
		the values are 0 or 1 which indicates whether the edge in question contains a halo
		(=1) or if it corresponds to the actual image edge (=0)
		if there is a halo we offset the correspnding loop boundary by the kernel halfsize
		
	*/
	
	register size_t hsize_h = k->halfsize[0];
	register size_t hsize_v = k->halfsize[1];
	
	
	size_t bound_l = halos[0];
	size_t bound_u = halos[1];
	size_t bound_r = xdim - halos[2];
	size_t bound_d = ydim - halos[3];
	
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite image lines no longer needed
	*/
	
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
			
			//printf("p %d line %d img %d\n",omp_get_thread_num(),(i - bound_u)%buffer_lines,i);


			if (i>=bound_u + buffer_lines) {
				memcpy( &image[(i - buffer_lines)*xdim + bound_l] , (&linebuf[line_idx]), linebuf_length*sizeof(int16_t));
			}

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				//are we close to the top/bottom edges?
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
				linebuf[line_idx + j - bound_l] = (uint16_t)accum;
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

			//are we close to the left/right edges?
			normflagy = ( (i<hsize_v) || (i>=(ydim - hsize_v) ) );

			for (size_t j=bound_l; j<bound_r; ++j) {
				//are we close to the top/bottom edges?
				normflagx = normflagy +  ( (j<hsize_h) || (j>=(xdim - hsize_h ) ) );

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
					//getchar();
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
				linebuf[line_idx + j - bound_l] = (int8_t)accum;
			}
		}
		//the buffer is full of lines still to write, 2 different cases
		if ((bound_d - bound_u)<=buffer_lines) {
			//not a single line of the buffer has been written since the "image lines" do not completely fill the buffer
			//therefore we write the buffer lines from the beginning until the last valid line
			size_t j = 0;
			for (size_t i = bound_u; i<bound_d; ++i) {
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
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
				memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
				j = (j + 1)%buffer_lines;
			}
		}
		free(linebuf);
	}
}


void  blur_halo_func_manager( pgm* input_img , const kernel_t* k, const int* halos) {
	if (k->halfsize[0]>=4) {
		pgm_blur_halo_unrolx4(input_img,k,halos);
	}
	else if (k->halfsize[0]>=2) {
		pgm_blur_halo_unrolx2(input_img,k,halos);
	}
	
}

