//#include <stdlib.h>
//#include <string.h>
//#include <stdio.h>
//#include <stdint.h>
//#include <math.h>
#include <pgm.h>
#include <p_grid.h>

#ifndef IMG_CELL_H
#define IMG_CELL_H

typedef struct {
	unsigned int coords[2];
	unsigned int size[2];
	unsigned int idx[2];
	unsigned int size_;
	unsigned int halos[4];
} img_cell;


#endif

void split_dimension(const int dim, const unsigned int nprocs, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth );
void get_cell_grid(p_grid* grid, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth);
void get_cell_1D(p_grid* grid, img_cell* proc_cell, const pgm* image, const unsigned int* kerhwidth);
int trim_halo( const img_cell* proc_cell, const char img_bytes);
void read_img_buffer( pgm* image , pgm* local_image, img_cell* cell_halo);
void write_img_buffer( pgm* image , pgm* local_image, img_cell* cell_halo, img_cell* cell_nohalo);


