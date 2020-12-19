#include <kernel_t.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>



int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int kernel_size, const float kernel_weight) {
	
	//initialise the kernel
	switch  (kernel_type) {
		case 0:
			average_kernel(k, kernel_size);
			break;
		case 1:
			//check if the weight was given as argument
			if (kernel_weight<0) {
				printf("Error: weighted kernel chosen but the weight was not specified.\n");	
				return -1;
			}
			else {
				weighted_kernel(k, kernel_size, kernel_weight);
				break;
				
			}
		case 2 :
			gaussian_kernel_simple(k, kernel_size);
			//gaussian_kernel(k, kernel_size);
			break;
		default:
			printf("Error: unknown kernel type.\n");
			printf("Currently supported: 0 (average), 1 (weighted), 2 (gaussian).\n");
			return -1;
			
	}
	printf("Kernel initialised.\n");
	
	return 0;
}


void average_kernel(kernel_t* k, const unsigned int kernel_size) {
	
	int kernel_size2 = kernel_size*kernel_size;
	
	double* kernel_matrix = (double*)calloc( kernel_size2 , sizeof(double));
	
	for (unsigned int i=0; i<kernel_size; ++i) {
		for (unsigned int j=0; j<kernel_size; ++j) {
			kernel_matrix[i*kernel_size + j] = (float)1/kernel_size2;
		}
	}
	
	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}


void weighted_kernel(kernel_t* k, const unsigned int kernel_size, const float kernel_weight) {
	
	int kernel_size2 = kernel_size*kernel_size;
	
	float w = ( 1 - kernel_weight)/(kernel_size2-1);
	
	double* kernel_matrix = (double*)calloc( kernel_size*kernel_size , sizeof(double));
	
	for (unsigned int i=0; i<kernel_size; ++i) {
		for (unsigned int j=0; j<kernel_size; ++j) {
			kernel_matrix[i*kernel_size + j] = w;
		}
	}
	
	int mid_idx = (kernel_size-1)/2;
	kernel_matrix[mid_idx*(kernel_size+1)] = kernel_weight;
	
	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}


//this uses the binomial coefficients
//uses the tgamma function in math.h for the factorial
void gaussian_kernel_simple(kernel_t* k, const unsigned  int kernel_size) {
	
	float binomial[kernel_size];
	float newval, norm=0, num=tgamma(kernel_size );

	for (unsigned int k=0; k<kernel_size; ++k) {
		newval = num/ ( tgamma(k + 1)*tgamma(kernel_size - k ) );
		binomial[k]	= newval;
		norm += newval;
	}
	norm = norm*norm;
	
	
	double* kernel_matrix = (double*)calloc( kernel_size*kernel_size , sizeof(double));

	for (unsigned int i=0; i<kernel_size; ++i) {
		for (unsigned int j=0; j<kernel_size; ++j) {
			kernel_matrix[i*kernel_size + j] =  binomial[i]*binomial[j]/norm;
		}
	}

	
	k->ker = kernel_matrix;
	k->size = kernel_size;
}




//void gaussian_kernel(kernel_t* k, const int kernel_size) {
//	
//	int mid_idx = (kernel_size-1)/2;
//	
//	float* kernel_matrix = (float*)calloc( kernel_size*kernel_size , sizeof(float));
//	
//	float stdev = 2*pow( (kernel_size/2) ,2);
//	int x2,y2;
//	float newval,norm;
//	for ( unsigned int i=0; i<kernel_size; ++i) {
//		for (unsigned int j=0; j<kernel_size; ++j) {
//			x2 = pow( (i - mid_idx),2);
//			y2 = pow( (j - mid_idx),2);
//			newval = (float) exp( -(x2 + y2)/stdev  )/ (M_PI*stdev);
//			norm+= newval;
//			kernel_matrix[i*kernel_size + j] =  newval;
//		}
//	}
//	for (int i=0; i<kernel_size; ++i) {
//		for (int j=0; j<kernel_size; ++j) {
//			//kernel_matrix[i*kernel_size + j] = kernel_matrix[i*kernel_size + j]/norm;
//		}
//	}
//	
//	
//	k->ker = kernel_matrix;
//	k->size = kernel_size;
//}


void delete_kernel( kernel_t* k) {
	free(k->ker);
	k->ker=NULL;
	k->size=0;
}