#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>



int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int ker_s, const float kernel_weight) {
	
	//initialise the kernel
	switch  (kernel_type) {
		case 0:
			average_kernel(k, ker_s);
			break;
		case 1:
			//check if the weight was given as argument
			if (kernel_weight<0) {
				printf("Error: weighted kernel chosen but the weight was not specified.\n");	
				return -1;
			}
			else {
				weighted_kernel(k, ker_s, kernel_weight);
				break;
				
			}
		case 2 :
			gaussian_kernel_simple(k, ker_s);
			//gaussian_kernel(k, ker_s);
			break;
		default:
			printf("Error: unknown kernel type.\n");
			printf("Currently supported: 0 (average), 1 (weighted), 2 (gaussian).\n");
			return -1;
			
	}
	
	get_kernel_normalisations(k);
	
	printf("Kernel initialised.\n");
	
	return 0;
}


void average_kernel(kernel_t* k, const unsigned int ker_s) {
	
	int ker_s2 = ker_s*ker_s;
	
	double* kernel_matrix = (double*)calloc( ker_s2 , sizeof(double));
	
	for (unsigned int i=0; i<ker_s; ++i) {
		for (unsigned int j=0; j<ker_s; ++j) {
			kernel_matrix[i*ker_s + j] = (float)1/ker_s2;
		}
	}
	
	
	k->ker = kernel_matrix;
	k->size = ker_s;
}


void weighted_kernel(kernel_t* k, const unsigned int ker_s, const float kernel_weight) {
	
	int ker_s2 = ker_s*ker_s;
	
	float w = ( 1 - kernel_weight)/(ker_s2-1);
	
	double* kernel_matrix = (double*)calloc( ker_s*ker_s , sizeof(double));
	
	for (unsigned int i=0; i<ker_s; ++i) {
		for (unsigned int j=0; j<ker_s; ++j) {
			kernel_matrix[i*ker_s + j] = w;
		}
	}
	
	int mid_idx = (ker_s-1)/2;
	kernel_matrix[mid_idx*(ker_s+1)] = kernel_weight;
	
	
	k->ker = kernel_matrix;
	k->size = ker_s;
}


//this uses the binomial coefficients
//uses the tgamma function in math.h for the factorial
void gaussian_kernel_simple(kernel_t* k, const unsigned  int ker_s) {
	
	float binomial[ker_s];
	float newval, norm=0, num=tgamma(ker_s );

	for (unsigned int k=0; k<ker_s; ++k) {
		newval = num/ ( tgamma(k + 1)*tgamma(ker_s - k ) );
		binomial[k]	= newval;
		norm += newval;
	}
	norm = norm*norm;
	
	
	double* kernel_matrix = (double*)calloc( ker_s*ker_s , sizeof(double));

	for (unsigned int i=0; i<ker_s; ++i) {
		for (unsigned int j=0; j<ker_s; ++j) {
			kernel_matrix[i*ker_s + j] =  binomial[i]*binomial[j]/norm;
		}
	}

	
	k->ker = kernel_matrix;
	k->size = ker_s;
}



//void gaussian_kernel(kernel_t* k, const int ker_s) {
//	
//	int mid_idx = (ker_s-1)/2;
//	
//	float* kernel_matrix = (float*)calloc( ker_s*ker_s , sizeof(float));
//	
//	float stdev = 2*pow( (ker_s/2) ,2);
//	int x2,y2;
//	float newval,norm;
//	for ( unsigned int i=0; i<ker_s; ++i) {
//		for (unsigned int j=0; j<ker_s; ++j) {
//			x2 = pow( (i - mid_idx),2);
//			y2 = pow( (j - mid_idx),2);
//			newval = (float) exp( -(x2 + y2)/stdev  )/ (M_PI*stdev);
//			norm+= newval;
//			kernel_matrix[i*ker_s + j] =  newval;
//		}
//	}
//	for (int i=0; i<ker_s; ++i) {
//		for (int j=0; j<ker_s; ++j) {
//			//kernel_matrix[i*ker_s + j] = kernel_matrix[i*ker_s + j]/norm;
//		}
//	}
//	
//	
//	k->ker = kernel_matrix;
//	k->size = ker_s;
//}


void get_kernel_normalisations(kernel_t* k) {
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
	
	double* kernel_norm = (double*)calloc( ker_s*ker_s , sizeof(double));
	
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
			kernel_norm[ker_s*(ker_hsize + offs_u + offs_d) +  ker_hsize + offs_l + offs_r] = normc;
				
		}
	}
	
	k->kernorm = kernel_norm;
}


void delete_kernel( kernel_t* k) {
	free(k->ker);
	k->ker=NULL;
	k->size=0;
}