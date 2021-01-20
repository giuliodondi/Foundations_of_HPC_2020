

#ifndef P_GRID_H
#define P_GRID_H

typedef struct {
	int size[2];
} p_grid;


#endif

void build_grid( p_grid* grid, int p );
int* get_grid_coords( p_grid* grid, int proc_id ); 