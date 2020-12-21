#include <stdint.h>

#ifndef KERNEL_H
#define KERNEL_H

//kernel struct definition
typedef struct {
	double* ker;
	double* kernorm;
	int size;
} kernel_t;


int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int kernel_size, const float kernel_weight);
void average_kernel(kernel_t* k);
void weighted_kernel(kernel_t* k, const float kernel_weight);
void gaussian_kernel_simple(kernel_t* k);
//void gaussian_kernel(kernel_t* k);

void get_kernel_normalisations(kernel_t* k);

void delete_kernel( kernel_t* k) ;

#endif