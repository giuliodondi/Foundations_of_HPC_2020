#ifndef KERNEL_H
#define KERNEL_H

//kernel struct definition
typedef struct {
	float* ker;
	int size;
} kernel;



void average_kernel(kernel* k, int kernel_size);
void weighted_kernel(kernel* k, int kernel_size, float kernel_weight);
void gaussian_kernel(kernel* k, int kernel_size);

void delete_kernel( kernel* k) ;

#endif