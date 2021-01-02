#include <pgm.h>
#include <img_cell.h>

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