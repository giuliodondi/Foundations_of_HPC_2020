#include <kernel_t.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdint.h>

int min( const int a, const int b) {
	if (a<b) {return a;}
	else {return b;}
}

int max( const int a, const int b) {
	if (a>b) {return a;}
	else {return b;}
}


void print_usage(char **argv) {
	printf("Usage: %s -input input_img.pgm -kernel_type t -kernel-size s (optional) -output output_img.pgm -kernel-weight w.\n",argv[0]);
}

char is_number( char arg[] ) {
	int i=0;
	if (arg[0] == '-' ) {
		++i	;
	}
	for (; arg[i] != 0; ++i ) {
		if (arg[i] == '.' ) {
			++i	;
			continue;
		}
		if (! isdigit(arg[i]) ) {
			return 0;	
		}
	}
	return 1;
}

//handles the program parameters and initialises the kernel struct
int8_t read_params_initialise_kernel( int argc, char **argv , char* infile, char* outfile , kernel_t* k ) {
	
	char kernel_file_flag=0;
	char kernel_fname[80] = "";
	int kernel_type=-1;
	int kernel_size[2]={-1,-1};
	int tmp;
	float kernel_weight = -1;
	

	
	//switch doesn't work with strings unfortunately
	for( int arg = 1; arg < argc; arg+=2 )  {
		if (strcmp(argv[arg], "-input")==0 ) {
			strcpy(infile,argv[arg+1]);
		}
		else if (strcmp(argv[arg], "-output")==0 ) {
			strcpy(outfile,argv[arg+1]);	
		}
		else if (strcmp(argv[arg], "-kernel-file")==0 ) {
			kernel_file_flag = 1;
			strcpy(kernel_fname,argv[arg+1]);
		}
		else if (strcmp(argv[arg], "-kernel-type")==0 ) {
			if (is_number(argv[arg + 1]) ) {
				kernel_type = atoi(argv[arg+1]);
			}
			else {
				printf("Kernel type is not a number.\n");	
				return -1;
			}
		}
		else if (strcmp(argv[arg], "-kernel-size")==0 ) {
			tmp=0;
			//get the first kernel dimension
			if (is_number(argv[arg + 1]) ) {
				tmp = atoi(argv[arg+1]);
			}
			else {
				printf("Kernel size is not a number.\n");	
				return -1;
			}
			
			//check kernel size is odd	
			if ( (tmp%2)==0 || (tmp<3)) {
				printf("Error: Kernel size(s) must be an odd integer 3 or greater.\n");
				return -1;
			}
			kernel_size[0]=tmp;
			
			//see if the next argument is a second integer, the second kernel dimension
			if (is_number(argv[arg + 2]) ) {
				tmp = atoi(argv[arg+2]);
				if ( (tmp%2)==0 || (tmp<3)) {
					printf("Error: Kernel size(s) must be an odd integer 3 or greater.\n");
					printf("Defaulting to a square kernel of size %d.\n", kernel_size[0] );
					tmp = kernel_size[0];
				}
				++arg;
			} 			
			
			kernel_size[1]=tmp;
			
			
		}
		else if (strcmp(argv[arg], "-kernel-weight")==0 ) {
			if (is_number(argv[arg + 1]) ) {
				tmp = atoi(argv[arg+1]);
			}
			else {
				printf("Kernel weight is not a number.\n");	
				return -1;
			}
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
	
	if (kernel_file_flag) {
		
		if (kernel_init_from_file( k, kernel_fname) == -1 ) {
			printf("Error during kernel initialisation.\n");	
			return -1;
		}
		
	} else {
		
		if (kernel_type == -1 ) {
			printf("The kernel type was not specified.\n");	
			print_usage(argv);
			return -1;
		}

		if (kernel_size[0] == -1 ) {
			printf("The kernel size was not specified.\n");	
			print_usage(argv);
			return -1;
		}

		if (kernel_init( k, kernel_type, (unsigned int*) kernel_size, kernel_weight) == -1 ) {
			printf("Error during kernel initialisation.\n");	
			return -1;
		}
	}
	
	
	

	
	return 0;
	
}


