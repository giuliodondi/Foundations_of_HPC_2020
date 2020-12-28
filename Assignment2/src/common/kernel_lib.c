#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>


int8_t alloc_kernel( kernel_t* k, const unsigned int ker_s ) {
	
	k->size = ker_s;
	int ker_size2 = ker_s*ker_s;
	k->ker = (double*)calloc( ker_size2 , sizeof(double));	
	k->kernorm = (double*)calloc( ker_size2 , sizeof(double));
	if ( k->ker && k->kernorm ) {
		return 0;
	} else {
		return -1;	
	}
	
}


int8_t copy_kernel(kernel_t* new_ker, kernel_t *old_ker) {
	
	if (alloc_kernel(new_ker, old_ker->size)== -1 ) {
		return -1;
	}
	
	int ker_size2 = old_ker->size*old_ker->size;
	memcpy( new_ker->ker , old_ker->ker , ker_size2*sizeof(double) );
	memcpy( new_ker->kernorm , old_ker->kernorm , ker_size2*sizeof(double) );
	return 0;
}



void delete_kernel( kernel_t* k) {
	free(k->ker);
	free(k->kernorm);
	k->ker=NULL;
	k->kernorm=NULL;
	k->size=0;
}

int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int ker_s, const float kernel_weight) {
	
	alloc_kernel( k , ker_s);
	
	//initialise the kernel matrix
	switch  (kernel_type) {
		case 0:
			average_kernel(k);
			break;
		case 1:
			//check if the weight was given as argument
			if (kernel_weight<0) {
				printf("Error: weighted kernel chosen but the weight was not specified.\n");	
				return -1;
			}
			else {
				weighted_kernel(k, kernel_weight);
				break;
				
			}
		case 2 :
			gaussian_kernel_simple(k);
			//gaussian_kernel(k);
			break;
		default:
			printf("Error: unknown kernel type.\n");
			printf("Currently supported: 0 (average), 1 (weighted), 2 (gaussian).\n");
			return -1;
			
	}
	
	
	kernel_normalisations(k);

	
	return 0;
}


void average_kernel(kernel_t* k) {
	
	double* kernel_matrix = k->ker;
	const size_t ker_s2 = (k->size*k->size);
	
	for (size_t i=0; i<ker_s2; ++i) {
		kernel_matrix[i] = (float)1/ker_s2;
	}
	
	
	// = kernel_matrix;
	//k->size = ker_s;
}


void weighted_kernel(kernel_t* k, const float kernel_weight) {
	
	double* kernel_matrix = k->ker;
	const size_t ker_s2 = (k->size*k->size);
	
	float w = ( 1 - kernel_weight)/(ker_s2-1);
	
	
	for (size_t i=0; i<ker_s2; ++i) {
		kernel_matrix[i] = w;
	}
	
	int mid_idx = (k->size-1)/2;
	kernel_matrix[mid_idx*(k->size+1)] = kernel_weight;
	
	
	//k->ker = kernel_matrix;
}


//this uses the binomial coefficients
//uses the tgamma function in math.h for the factorial
void gaussian_kernel_simple(kernel_t* k) {
	
	double* kernel_matrix = k->ker;
	const size_t ker_s = k->size;
	
	float binomial[ker_s];
	float newval, norm=0, num=tgamma(ker_s );

	for (size_t k=0; k<ker_s; ++k) {
		newval = num/ ( tgamma(k + 1)*tgamma(ker_s - k ) );
		binomial[k]	= newval;
		norm += newval;
	}
	norm = norm*norm;


	for (size_t i=0; i<ker_s; ++i) {
		for (size_t j=0; j<ker_s; ++j) {
			kernel_matrix[i*ker_s + j] =  binomial[i]*binomial[j]/norm;
		}
	}

	
	k->ker = kernel_matrix;
	k->size = ker_s;
}




void kernel_normalisations(kernel_t* k) {
	/*
	calculates a matrix of normalisation constants to use vor vignetting removal
	explanation : consider a 3x3 kernel and overlap its central point on the top edge and corners
	of the image -  here's my attempt at a pictorial representation:
	
	o o o      o o o       o o o                              o o o      o o o       o o o
	o + + ---  + + +  ---  + + o              -->             o o o ---  o o o  ---  o o o 
	o + +      + + +       + + o                              o o x      o x o       x o o
	
	
	only the kernel entries marked with '+' should count , in the blurring function we use four indices to calculate
	the four extreme indices and loop only over the right subset
	the same four indices can be used to identify a unique cell in a similar 3x3 matrix, marked with 'x'
	
	here we compute a normalisation matrix of the same size as the kernel, compute all the normalisation constants and store 
	them in their unique place, so that vignetting removal entails just a memory access and a float division
	*/
	
	register const int ker_s = k->size ;	
	register const int ker_hsize = ( k->size - 1)/2 ;
	
	double* kernel = k->ker;
	
	register int offs_l, offs_r, offs_u, offs_d;
	register double normc;
	
	double* kernel_norm = k->kernorm;
	
	for (int i=0; i<ker_s; ++i) {
		for (int j=0; j<ker_s; ++j) {
			
			offs_l = -min( ker_hsize, j );
			offs_u = -min( ker_hsize, i );
			offs_r = min( ker_hsize, ker_s - 1 - j );
			offs_d = min( ker_hsize, ker_s - 1 - i );
						
			normc=0;
			for (int k = offs_u; k<= offs_d; ++k) {
				for (int t =  offs_l; t<= offs_r; ++t) {
					normc += kernel[ ker_s*(ker_hsize + k) +  ker_hsize + t  ] ;
				}

			}
			kernel_norm[ker_s*(ker_hsize + offs_u + offs_d) +  ker_hsize + offs_l + offs_r] = 1/normc;
				
		}
	}

}


