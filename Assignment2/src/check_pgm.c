#include <pgm.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>



int main( int argc, char **argv )  { 
		
	
	char file1[80] = "";
	char file2[80] = "";
	char* outfile = "checklog.txt";
	long int header_offs=0;
	
	if (argc<2) {
		printf("Error, not enough arguments. Specify two image names to compare.\n");	
	}
	
	strcpy(file1,argv[1]);
	strcpy(file2,argv[2]);


	pgm image1 = new_pgm();
	pgm image2 = new_pgm();


	//read the blurred image
	if (read_pgm_header( &image1 , file1, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image1);
		return -1;
	}
	if (allocate_pgm_memory( &image1)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image1);
		return -1;
	}
	if (read_pgm_data( &image1 , file1, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image1);
		return -1;
	}
	
	//read the reference image
	if (read_pgm_header( &image2 , file2, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image2);
		return -1;
	}
	if (allocate_pgm_memory( &image2)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image2);
		return -1;
	}
	if (read_pgm_data( &image2 , file2, &header_offs)== -1 ) {
		printf("Aborting.\n");
		clear_pgm( &image2);
		return -1;
	}
	
	endian_swap(&image1);
	endian_swap(&image2);
	
	FILE* f; 
	f = fopen(outfile, "w"); 
	
	fprintf(f,"Image 1 : \"%s\"\n",file1);
	fprintf(f,"Image 2 : \"%s\"\n",file2);
	
	fclose(f); 
	
	compare_pgm( &image1, &image2, outfile) ;

	clear_pgm( &image1);
	clear_pgm( &image2);
    return 0;
} 