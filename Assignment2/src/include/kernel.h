#ifndef KERNEL_H
#define KERNEL_H

//kernel struct definition
typedef struct {
	float* ker;
	int size;
} kernel_t;



void average_kernel(kernel_t* k, int kernel_size);
void weighted_kernel(kernel_t* k, int kernel_size, float kernel_weight);
//void gaussian_kernel(kernel_t* k, int kernel_size);
void gaussian_kernel_simple(kernel_t* k, int kernel_size);

void delete_kernel( kernel_t* k) ;

#endif