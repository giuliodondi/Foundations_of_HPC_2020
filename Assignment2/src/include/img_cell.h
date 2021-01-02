#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#ifndef IMG_CELL_H
#define IMG_CELL_H

typedef struct {
	unsigned int size[2];
	unsigned int idx[2];
	unsigned int size_;
	//int idx;
	unsigned int halos[4];
} img_cell;


#endif


void get_cell_1D(const int nprocs, const int proc_id, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth);
int trim_halo_1D( const img_cell* proc_cell, const char img_bytes);



