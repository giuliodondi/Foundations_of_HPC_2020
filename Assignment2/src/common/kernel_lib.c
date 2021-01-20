#include <kernel_t.h>
#include <common_headers.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

//allocate the kernel main matrix and the normalisation kernel
//iniitalise the size and halfsize arrays
int8_t alloc_kernel( kernel_t* k, const unsigned int* restrict ker_s ) {
	
	k->size[0] = ker_s[0];
	k->size[1] = ker_s[1];
	const int ker_size2 = ker_s[0]*ker_s[1];
	k->halfsize[0] = (k->size[0] - 1)/2;
	k->halfsize[1] = (k->size[1] - 1)/2;
	k->ker = (double*)calloc( ker_size2 , sizeof(double));	
	k->kernorm = (double*)calloc( ker_size2 , sizeof(double));
	if ( k->ker && k->kernorm ) {
		return 0;
	} else {
		return -1;	
	}
	
}

//copy over a kernel
int8_t copy_kernel(kernel_t* new_ker, const kernel_t *old_ker) {
	
	if (alloc_kernel(new_ker, old_ker->size)== -1 ) {
		return -1;
	}
	
	const int ker_size2 = old_ker->size[0]*old_ker->size[1];
	memcpy( new_ker->ker , old_ker->ker , ker_size2*sizeof(double) );
	memcpy( new_ker->kernorm , old_ker->kernorm , ker_size2*sizeof(double) );
	return 0;
}


//free the memory and reset danglign pointers
void delete_kernel( kernel_t* k) {
	free(k->ker);
	free(k->kernorm);
	k->ker=NULL;
	k->kernorm=NULL;
	k->size[0] = 0;
	k->size[1] = 0;
	k->halfsize[0] = 0;
	k->halfsize[1] = 0;
}


//read a kernel from formatted file
int8_t kernel_init_from_file(kernel_t* k, const  char* kernel_fname ) {
	FILE* kernel_file; 
	kernel_file = fopen(kernel_fname, "r"); 
	
	char   *line = NULL;
	char * kerval;
	size_t  a, n=0,i = 0;
	unsigned int kernel_size[2]={0,0};

	// skip all the comments
	a = getline( &line, &n, kernel_file);
	while ( (a > 0) && (line[0]=='#') ) {
    	a = getline( &line, &n, kernel_file);
	}	
	if (a<=0) {
	  	printf("Error while reading the kernel file header.\n");
	  	free( line );
		fclose(kernel_file);
	  	return -1;
	}
	
	// read the kernel sizes
	a = fscanf(kernel_file, "%d%*c %d%*c", &kernel_size[0], &kernel_size[1] );
	
	//check kernel dimensions
	if ( (kernel_size[0]%2)==0 || (kernel_size[0]<3) || (kernel_size[1]%2)==0 || (kernel_size[1]<3) ) {
		printf("Invalid kernel dimensions specified.\n");
		printf("Kernel size(s) must be an odd integer 3 or greater.\n");
		free( line );
		fclose(kernel_file);
	  	return -1;
		
	}
	
	alloc_kernel( k , kernel_size);
	
	//read line by line for the kernel values separated either by comma or space
	while ( fgets( line, sizeof(line), kernel_file ) != NULL ) {
		kerval=strtok(line," ,");
		while (kerval!=NULL) {
			k->ker[i] = atof(kerval);
			++i;
			kerval=strtok(NULL," ,");
			
		}
	  
	}
	
	fclose(kernel_file);
	free( line );
	
	normalise_kernel(k);
	kernel_normalisations(k);
	
	return 0;
}

//initialise kernel of the standard kinds
int8_t kernel_init(kernel_t* k, const unsigned int kernel_type, const unsigned int* restrict ker_s, const float kernel_weight) {
	
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


//initialise the averaging kernel
void average_kernel(kernel_t* k) {
	
	double* kernel = k->ker;
	const size_t ker_s2 = k->size[0]*k->size[1];
	
	for (size_t i=0; i<ker_s2; ++i) {
		kernel[i] = 1/(double)ker_s2;
	}
	

}

//initialise the weighted kernel
void weighted_kernel(kernel_t* k, const double kernel_weight) {
	
	double* kernel = k->ker;
	const size_t ker_s2 = (k->size[0]*k->size[1]);
	
	const double w = ((double)( 1 - kernel_weight))/((double)(ker_s2-1));
	
	
	for (size_t i=0; i<ker_s2; ++i) {
		kernel[i] = w;
	}
	
	kernel[k->size[0]*k->halfsize[1] + k->halfsize[0]] = kernel_weight;

}


//initialise the gaussian kernel with binomial coefficients
//uses the tgamma function in math.h for the factorial
void gaussian_kernel_simple(kernel_t* k) {
	
	double* kernel = k->ker;
	
	double binomial_h[k->size[0]];
	double binomial_v[k->size[1]];
	double newval, norm1=0, norm2=0, num1=tgamma(k->size[0] ), num2=tgamma(k->size[1] );

	for (size_t i=0; i<k->size[0]; ++i) {
		newval = num1/ ( tgamma(i + 1)*tgamma((int)k->size[0] - i ) );
		binomial_h[i]	= newval;
		norm1 += newval;
	}
	for (size_t i=0; i<k->size[1]; ++i) {
		newval = num2/ ( tgamma(i + 1)*tgamma((int)k->size[1] - i ) );
		binomial_v[i]	= newval;
		norm2 += newval;
	}
	norm1 = 1/(norm1*norm2);


	for (size_t i=0; i<k->size[1]; ++i) {
		for (size_t j=0; j<k->size[0]; ++j) {
			kernel[i*k->size[0] + j] =  binomial_v[i]*binomial_h[j]*norm1;
		}
	}

	k->ker = kernel;
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
	
	//const int ker_s = k->size ;	
	//const int ker_hsize = ( k->size - 1)/2 ;
	
	
	double* kernel = k->ker;
	
	register int offs_l, offs_r, offs_u, offs_d;
	register double normc;
	
	double* kernel_norm = k->kernorm;
	
	for (size_t i=0; i<k->size[1]; ++i) {
		for (size_t j=0; j<k->size[0]; ++j) {
			
			offs_l = -min( k->halfsize[0], j );
			offs_u = -min( k->halfsize[1], i );
			offs_r = min( k->halfsize[0], k->size[0] - 1 - j );
			offs_d = min( k->halfsize[1], k->size[1] - 1 - i );

			//x2 loop unrolling 
			normc=0;			
			if (( offs_r - offs_l)%2) {
				for (int u = offs_u; u<= offs_d; ++u) {
					for (int t =  offs_l ; t<= offs_r; t+=2) {
						normc += kernel[ k->size[0]*(k->halfsize[1] + u) +  k->halfsize[0] + t  ] 
							+ kernel[ k->size[0]*(k->halfsize[1] + u) +  k->halfsize[0] + t + 1 ] ;
					}

				}
			}
			else {
				for (int u = offs_u; u<= offs_d; ++u) {
					normc += kernel[ k->size[0]*(k->halfsize[1] + u) +  k->halfsize[0] + offs_l  ] ;
					for (int t =  offs_l + 1; t<= offs_r; t+=2) {
						normc += kernel[ k->size[0]*(k->halfsize[1] + u) +  k->halfsize[0] + t  ] 
							+ kernel[ k->size[0]*(k->halfsize[1] + u) +  k->halfsize[0] + t + 1 ] ;
					}

				}
				
			}
			kernel_norm[k->size[0]*(k->halfsize[1] + offs_u + offs_d) +  k->halfsize[0] + offs_l + offs_r] = 1/normc;
			
		}
	}

}

//normalises a given kernel
//needed for the kernel from file, not neded for the built-in kernels
void normalise_kernel(kernel_t* k) {
	double* kernel = k->ker;
	double norm=0;
	for (size_t i=0; i<k->size[1]; ++i) {
		for (size_t j=0; j<k->size[0]; ++j) {
			norm += kernel[i*k->size[0] + j];
		}
	}
	for (size_t i=0; i<k->size[1]; ++i) {
		for (size_t j=0; j<k->size[0]; ++j) {
			kernel[i*k->size[0] + j] /=  norm;
		}
	}
	
}
