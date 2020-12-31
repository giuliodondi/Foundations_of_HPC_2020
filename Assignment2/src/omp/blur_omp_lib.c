#include <pgm.h>
#include <kernel_t.h>
#include <img_cell.h>
#include <common_headers.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>


//one-dimensional image splitting
void get_cell_1D(const int nprocs, const int proc_id, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth) {
	
	int childlines = ceil(((float)image->size[1])/((float)nprocs) );
	int masterlines = image->size[1] - (nprocs - 1)*childlines;
	int lastchildlines = childlines;
	if (masterlines<0) {
		childlines -=1;
		lastchildlines = image->size[1] - (nprocs - 1)*childlines;
		masterlines = childlines;
	}
	
	#ifdef INFO
	if (proc_id==0) {
		printf("Master: %d , Child: %d , Lastchild %d , Halowidth: %d\n",masterlines,childlines,lastchildlines,kerhwidth[1]);
	}
	#endif

	proc_cell->size[0] = image->size[0];
	
	if (proc_id == 0) {
		proc_cell->idx[0]=0;
		proc_cell->idx[1]=0;
		
		proc_cell->size[1] = masterlines;
		proc_cell->halos[0] = 0;
		proc_cell->halos[1] = 0;
		proc_cell->halos[2] = 0;
		proc_cell->halos[3] = 0;
		
		//if there is just one worker it takes care of the entire image
		//no halo needs to be included
		if (nprocs>1 ) {
			proc_cell->size[1] += kerhwidth[1];
			proc_cell->halos[3] = kerhwidth[1];
		} 
		
	} 	else {
		
		proc_cell->idx[0]=0;
		proc_cell->idx[1]=masterlines + (proc_id - 1)*childlines ;
		
		
		proc_cell->halos[0] = 0;
		proc_cell->halos[2] = 0;
		
		//check if the cell is closer to the top edge of the image by the kernel half height
		//if so the halo is the number of lines between the top edge and the cell
		//otherwise we can take the entire kernel half height as the halo
		if (  proc_cell->idx[1] < kerhwidth[1]) {
			proc_cell->halos[1] = proc_cell->idx[1];
			proc_cell->idx[1] = 0;
		}
		else {
			proc_cell->idx[1]-=kerhwidth[1];
			proc_cell->halos[1] = kerhwidth[1];	
		}
		
		proc_cell->size[1] = proc_cell->halos[1];
			
		//similar check for the lower halo
		//but not for the very last cell
		if (proc_id < (nprocs - 1)) {
			proc_cell->size[1] += childlines;
			unsigned int cell_lastline = proc_cell->idx[1] + childlines;

			if ( image->size[1] < cell_lastline + kerhwidth[1] ) {
				proc_cell->halos[3] =  image->size[1] - cell_lastline;
			} else {
				proc_cell->halos[3] =  kerhwidth[1];
			}


			

		} else {
			proc_cell->size[1] += lastchildlines; 
			proc_cell->halos[3] = 0;
		}
		
		proc_cell->size[1] += proc_cell->halos[3];
		
		
	}
	proc_cell->size_ = proc_cell->size[0]*proc_cell->size[1]*image->pix_bytes;
	
}




//returns the idx in the local cell buffer of the first "real" image pixel (i.e. not halo)

int trim_halo_1D( const img_cell* proc_cell, const char img_bytes) {
	//width of the actual pixels plus left and right halo if present
	int w = (proc_cell->halos[0] + proc_cell->halos[2]) + proc_cell->size[0];
	//skip the top halo rows if present and add the left halo on the first "actual" row if present
	return (w*proc_cell->halos[1] + proc_cell->halos[0])*img_bytes;
}






void pgm_blur_halo(  pgm* input_img , const kernel_t* k,  const unsigned int* halos) {
	
	
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
	
	
	register size_t bound_l = halos[0];
	register size_t bound_u = halos[1];
	register size_t bound_r = xdim - halos[2];
	register size_t bound_d = ydim - halos[3];
	
	/*
	register size_t bound_l = halos[0]*hsize_h;
	register size_t bound_u = halos[1]*hsize_v;
	register size_t bound_r = xdim - halos[2]*hsize_h;
	register size_t bound_d = ydim - halos[3]*hsize_v;
	*/
	
	
	/*
	store the blurred pixels in a buffer
	the buffer is as wide as the image and has half as many lines as the kernel
	including the middle line
	we roll through the buffer and fill the lines with the blurred pixels
	only when the buffer is full we start back at line 0 and memcpy it at line 0 of the image data
	we're safe that we only overwrite when the data is no longer needed
	
	at the end all the calculations are done but the buffer is full of data still to be written
	*/
	
	register const size_t buffer_lines = (hsize_v + 1);
	register size_t line_idx; 
	register int linebuf_length = (bound_r - bound_l);

	
	register int offs_l, offs_r, offs_u, offs_d;
	double normc;
	register double accum;
	register uint8_t normflagx,normflagy;

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
				linebuf[line_idx + j - bound_l] = (uint16_t)accum;
			}
		}
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
		for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {

			memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int16_t));
			j = (j + 1)%buffer_lines;
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
		//the buffer is full of lines still to write
		//but the next buffer line to write will be at an arbitrary index

		//this is the buffer line that contains the next line of data
		size_t j = (bound_d - bound_u - buffer_lines)%buffer_lines ;
		//write the data in the last lines of the image
		//when we reach the end of the buffer, wrap around as the last lines to write are at the top
		for (size_t i = bound_d - buffer_lines; i<bound_d ; ++i) {

			memcpy( &image[i*xdim + bound_l] , (&linebuf[linebuf_length*j]), linebuf_length*sizeof(int8_t));
			j = (j + 1)%buffer_lines;
		}
		free(linebuf);
	}
}



