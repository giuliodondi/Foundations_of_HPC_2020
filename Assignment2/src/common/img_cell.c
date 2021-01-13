#include <pgm.h>
#include <img_cell.h>

#include <pgm.h>
#include <p_grid.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>


void split_dimension(const int dim, const int nprocs, img_cell* proc_cell, const pgm* image, const int* kerhwidth ) {

	
	float baselines = ((float)image->size[dim])/((float)nprocs);
	
	int h1 = (int)floor(baselines);
	int h2 = (int)ceil(baselines);
	
	int n1=0, n2=0, tot_size;
	
	for (int i=1; i<nprocs; ++i) {
		n1=i;
		n2 =(nprocs - i) ;
		tot_size =n1*h1 + n2*h2 ;
		if (tot_size==image->size[dim]) {
			break;	
		}
	}
	
	
	#ifdef INFO
	if ( (proc_cell->coords[0]==0) && (proc_cell->coords[1]==0) ){
		printf("Dim: %d, Size 1 : %d x %d , Size 2 : %d x %d , Halowidth: %d\n",dim,n1, h1,n2, h2,kerhwidth[dim]);
	}
	#endif

	//first initialise the base cell size and index along the specified dimension
	if (proc_cell->coords[dim] < n1 ) {
		proc_cell->idx[dim]= proc_cell->coords[dim]*h1 ;
		proc_cell->size[dim]=h1;
	}
	else {
		proc_cell->idx[dim]= n1*h1 + (  proc_cell->coords[dim] - n1 )*h2  ;
		proc_cell->size[dim]=h2;
	}
	
	//now work out if halos are necessary, their size and update index and size
	
	//we check if the inclusion of the full kernel width overhangs the image edge
	//but ONLY for the cells that are not the frst and last along the dimension
	//this is valid if the kernel width is << image size along the dimension
	
	if (proc_cell->coords[dim] == 0 ) {
		//the first cell has no "top" halo
		proc_cell->halos[dim]=0;
		
		//if there is just one worker along the dimension
		//no halo needs to be included
		if (nprocs>1 ) {
			proc_cell->halos[dim + 2] = kerhwidth[dim];
		} else {
			proc_cell->halos[dim+2]=0;	
		}
		
	} else if (proc_cell->coords[dim] == nprocs - 1 ) {
		proc_cell->idx[dim]-=kerhwidth[dim];
		proc_cell->halos[dim]=kerhwidth[dim];
		
		proc_cell->halos[dim+2]=0;	
		
	} 

	else {
		
		//check if the cell is closer to the top edge of the image by the kernel half height
		//if so the halo is the number of lines between the top edge and the cell
		//otherwise we can take the entire kernel half height as the halo
		if (  proc_cell->idx[dim] < kerhwidth[dim]) {
			proc_cell->halos[dim] = proc_cell->idx[dim];
			proc_cell->idx[dim] = 0;
		}
		else {
			proc_cell->idx[dim]-=kerhwidth[dim];
			proc_cell->halos[dim] = kerhwidth[dim];	
		}
		
		//similar check for the lower halo
		int cell_lastline = proc_cell->idx[dim] + proc_cell->size[dim];
		if ( image->size[dim] < cell_lastline + kerhwidth[dim] ) {
			proc_cell->halos[dim+2] =  image->size[dim] - cell_lastline;
		} else {
			proc_cell->halos[dim+2] =  kerhwidth[dim];
		}
	}
	
	proc_cell->size[dim] += proc_cell->halos[dim] + proc_cell->halos[dim + 2];
	
}

//one-dimensional image splitting
void get_cell_1D(p_grid* grid, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth) {
	
	
	//dimension x 
	proc_cell->size[0] = image->size[0];
	proc_cell->idx[0] = 0;
	proc_cell->halos[0] = 0;
	proc_cell->halos[2] = 0;
	
	//dimension y
	split_dimension(1, grid->size[1], proc_cell, image, (int*)kerhwidth );
	
	proc_cell->size_ = proc_cell->size[0]*proc_cell->size[1]*image->pix_bytes;
}


// image splitting on a grid
void get_cell_grid(p_grid* grid, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth) {
	
	
	//dimension x 
	split_dimension(0, grid->size[0], proc_cell, image, (int*)kerhwidth );
	
	//dimension y
	split_dimension(1, grid->size[1], proc_cell, image, (int*)kerhwidth );
	//proc_cell->size[1] = image->size[1];
	//proc_cell->idx[1] = 0;
	//proc_cell->halos[1] = 0;
	//proc_cell->halos[3] = 0;
	
	proc_cell->size_ = proc_cell->size[0]*proc_cell->size[1]*image->pix_bytes;
}




//returns the idx in the local cell buffer of the first "real" image pixel (i.e. not halo)

int trim_halo( const img_cell* proc_cell, const char img_bytes) {
	//width of the actual pixels plus left and right halo if present
	int w = (proc_cell->halos[0] + proc_cell->halos[2]) + proc_cell->size[0];
	//skip the top halo rows if present and add the left halo on the first "actual" row if present
	return (w*proc_cell->halos[1] + proc_cell->halos[0])*img_bytes;
}


void read_img_buffer( pgm* image , pgm* local_image, img_cell* cell_halo) {
	
	//work out the beginning of the buffer in the original image
	register int img_idx = img_idx_convert(image, cell_halo->idx);
	register int cell_idx=0 ;
	
	int cell_line_size_halo=cell_halo->size[0]*image->pix_bytes;
	
	for (int i=cell_halo->idx[1]; i< (cell_halo->idx[1] + cell_halo->size[1]); ++i) {	
		
			memcpy( &(local_image->data[ cell_idx ]) , &(image->data[ img_idx] ) , cell_line_size_halo*sizeof(uint8_t) );
		
			cell_idx += cell_line_size_halo;
			img_idx += (image->size[0])*image->pix_bytes;
	}
	
}

void write_img_buffer( pgm* image , pgm* local_image, img_cell* cell_halo, img_cell* cell_nohalo) {
	
	
	//wortk out the beginning of the buffer in the original image
	register int img_idx = img_idx_convert(image, cell_nohalo->idx);
	register int cell_idx = (cell_halo->size[0]*cell_halo->halos[1] + cell_halo->halos[0])*image->pix_bytes;
	
	int cell_line_size_nohalo=cell_nohalo->size[0]*image->pix_bytes;
	int cell_line_size_halo=cell_halo->size[0]*image->pix_bytes;
	
	for (int i=cell_nohalo->idx[1]; i< (cell_nohalo->idx[1] + cell_nohalo->size[1]); ++i) {			
			memcpy(  &(image->data[ img_idx] ) , &(local_image->data[ cell_idx ]) , cell_line_size_nohalo*sizeof(uint8_t) );
		
			cell_idx += cell_line_size_halo;
			img_idx += image->size[0]*image->pix_bytes;
	}
	
}