
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