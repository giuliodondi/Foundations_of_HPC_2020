
#include <stdint.h>

#ifndef KERNEL_H
#define KERNEL_H

//kernel struct definition
typedef struct {
	double* ker;
	double* kernorm;
	unsigned int size[2];
	unsigned int halfsize[2];
} kernel_t;

//kernel func headers
int8_t alloc_kernel( kernel_t* k, const unsigned int* ker_s );
int8_t copy_kernel(kernel_t* new_ker, const kernel_t *old_ker) ;
void delete_kernel( kernel_t* k) ;
int8_t kernel_init_from_file(kernel_t* k, const  char* kernel_fname );
int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int* kernel_size, const float kernel_weight);

void average_kernel(kernel_t* k);
void weighted_kernel(kernel_t* k, const double kernel_weight);
void gaussian_kernel_simple(kernel_t* k);

void kernel_normalisations(kernel_t* k);
void normalise_kernel(kernel_t* k);


#endif