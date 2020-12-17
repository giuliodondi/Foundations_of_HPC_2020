#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <kernel.h>



int read_params_initialise_kernel( const int argc, const char **argv , char* infile, char* outfile , kernel* k ) {
	
	unsigned int kernel_type;
	unsigned int kernel_size;
	float kernel_weight = -1;
	
	//switch doesn't work with strings unfortunately
	for( int arg = 1; arg < argc; arg+=2 )  {
		if (strcmp(argv[arg], "-input")==0 ) {
			strcat(infile,argv[arg+1]);
		}
		else if (strcmp(argv[arg], "-output")==0 ) {
			strcat(outfile,argv[arg+1]);	
		}
		else if (strcmp(argv[arg], "-kernel-type")==0 ) {
			kernel_type = atoi(argv[arg+1]);
		}
		else if (strcmp(argv[arg], "-kernel-size")==0 ) {
			kernel_size = atoi(argv[arg+1]);
			//check kernel size is odd	
			if ( (kernel_size%2)==0) {
				printf("Error: Kernel size must be a positive odd integer.\n");
				return -1;
			}
		}
		else if (strcmp(argv[arg], "-kernel-weight")==0 ) {
			kernel_weight = atof(argv[arg+1]);
			//check the kernel weight is bounded
			//the check must be done here so that if the weight is not given as argument
			//we can tell because the value is initialised negative
			if ((kernel_weight<0) || (kernel_weight>1)) {
				printf("Kernel weight must be a positive float less than 1.\n");
				return -1;
			}
		}
		else {
			printf("Illegal or missing argument.\n");
			printf("Usage: %s -input path/name1.pgm -output path/name2.pgm .\n",argv[0]);
			return -1;
		  }
		 
		 
	}
	

	
	
	
	
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
			gaussian_kernel(k, kernel_size);
			break;
		default:
			printf("Error: unknown kernel type.\n");
			printf("Currently supported: 0 (average), 1 (weighted), 2 (gaussian).\n");
			return -1;
			
	}
	printf("Kernel initialised.\n");
	
	//check if file exists 
	if( access( infile, F_OK ) == 0 ) {
		printf("Processing file \"%s\".\n",infile); 
	} else {
		printf("The input file does not exists.\n");
		return -1;
	} 
	
	return 0;
	
}


