
#ifndef IMG_CELL_H
#define IMG_CELL_H

typedef struct {
	int width;
	int height;
	int idx;
	int size;
	char halos[4];
} img_cell;


#endif