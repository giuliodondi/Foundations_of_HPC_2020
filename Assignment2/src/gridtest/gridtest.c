#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef P_GRID_H
#define P_GRID_H

typedef struct {
	unsigned int size[2];
} p_grid;


#endif

void build_grid( p_grid* grid, int p ) {
	
	int p_x=1, p_y=p, b=(int)sqrt(p) + 1;
	
	for (int i=2; i < b; ++i) {
		
		if ( p%i==0) {
			
			p_x = i;
			p_y = p/i;
		}
	}
	
	grid->size[0] = p_x;
	grid->size[1] = p_y;
	
}


int* get_grid_coords( p_grid* grid, int proc_id ) {
	static int grid_coords[2];

	
	for (int i=0; i<grid->size[1]; ++i) {
		if (i*grid->size[0] <= proc_id) {
			grid_coords[1]=i;
		} else {
			break;
		}
	}
	
	grid_coords[0]= proc_id -  grid_coords[1]*grid->size[0];
	
	return grid_coords;
}

int main( int argc, char **argv ) {

	if (argc<3) {
		printf("error\n");
		return 0;
	}
	
	int p = atoi(argv[1]);
	int p_id = atoi(argv[2]);
	
	printf("%d %d\n",p, p_id);
	
	p_grid g;
	
	build_grid(&g,p);
	
	int* coords = get_grid_coords( &g, p_id );
	
	printf("grid %d x %d\n",g.size[0],g.size[1]);
	printf("coords %d x %d\n",coords[0],coords[1]);

}