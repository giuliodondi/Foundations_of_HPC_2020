#include <p_grid.h>
#include <math.h>


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
	
	unsigned int pid = (unsigned int) proc_id  ;
	
	if (grid->size[0]>1) {
		grid_coords[0] = pid%grid->size[0];
	} else {
		grid_coords[0] = 0;
	}
	grid_coords[1] = floor(pid/grid->size[0]);
		
	return grid_coords;
}