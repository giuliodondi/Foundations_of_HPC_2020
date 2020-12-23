#include <kernel_t.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

int min( const int a, const int b) {
	if (a<b) {return a;}
	else {return b;}
}





void print_usage(char **argv) {
	printf("Usage: %s -input input_img.pgm -kernel_type t -kernel-size s (optional) -output output_img.pgm -kernel-weight w.\n",argv[0]);
}

//handles the program parameters and initialises the kernel struct
int8_t read_params_initialise_kernel( int argc, char **argv , char* infile, char* outfile , kernel_t* k ) {
	
	int kernel_type=-1;
	int kernel_size=-1;
	float kernel_weight = -1;
	

	
	//switch doesn't work with strings unfortunately
	for( int arg = 1; arg < argc; arg+=2 )  {
		if (strcmp(argv[arg], "-input")==0 ) {
			strcpy(infile,argv[arg+1]);
		}
		else if (strcmp(argv[arg], "-output")==0 ) {
			strcpy(outfile,argv[arg+1]);	
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
			printf("Illegal argument.\n");
			print_usage(argv);
			return -1;
		  }
		 
		 
	}
	
	//check if file exists 
	if( !access( infile, F_OK ) == 0 ) {
		printf("The input file does not exists.\n");
		return -1;
	} 
	
	if (kernel_type == -1 ) {
		printf("The kernel type was not specified.\n");	
		print_usage(argv);
		return -1;
	}
	
	if (kernel_size == -1 ) {
		printf("The kernel size was not specified.\n");	
		print_usage(argv);
		return -1;
	}
	
	if (kernel_init( k, kernel_type, kernel_size, kernel_weight) == -1 ) {
		printf("Error during kernel initialisation.\n");	
		return -1;
	}
	
	
	

	
	return 0;
	
}


