#include <stdlib.h>
#include <malloc.h>
#include <kernel.h>
#include <math.h>



void average_kernel(kernel* k, int kernel_size) {
	
	int kernel_size2 = kernel_size*kernel_size;
	
	float* kernel_matrix = (float*)calloc( kernel_size2 , sizeof(float));
	
	for (int i=0; i<kernel_size; ++i) {
		for (int j=0; j<kernel_size; ++j) {
			kernel_matrix[i*kernel_size + j] = 1/kernel_size2;
		}
	}
	
	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}


void weighted_kernel(kernel* k, int kernel_size, float kernel_weight) {
	
	int kernel_size2 = kernel_size*kernel_size;
	
	float w = ( 1 - kernel_weight)/(kernel_size2-1);
	
	float* kernel_matrix = (float*)calloc( kernel_size*kernel_size , sizeof(float));
	
	for (int i=0; i<kernel_size; ++i) {
		for (int j=0; j<kernel_size; ++j) {
			kernel_matrix[i*kernel_size + j] = w;
		}
	}
	
	int mid_idx = (kernel_size-1)/2;
	kernel_matrix[mid_idx*(kernel_size+1)] = kernel_weight;
	
	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}


void gaussian_kernel(kernel* k, int kernel_size) {
	
	int mid_idx = (kernel_size-1)/2;
	
	float* kernel_matrix = (float*)calloc( kernel_size*kernel_size , sizeof(float));
	
	float stdev = 2*pow( (kernel_size/2) ,2);
	int x2,y2;
	float newval,norm;
	for (int i=0; i<kernel_size; ++i) {
		for (int j=0; j<kernel_size; ++j) {
			x2 = pow( (i - mid_idx),2);
			y2 = pow( (j - mid_idx),2);
			newval = (float) exp( -(x2 + y2)/stdev  )/ (M_PI*stdev);
			norm+= newval;
			kernel_matrix[i*kernel_size + j] =  newval;
		}
	}
	for (int i=0; i<kernel_size; ++i) {
		for (int j=0; j<kernel_size; ++j) {
			//kernel_matrix[i*kernel_size + j] = kernel_matrix[i*kernel_size + j]/norm;
		}
	}
	
	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}


void delete_kernel( kernel* k) {
	free(k->ker);
	k->ker=NULL;
	k->size=0;
}