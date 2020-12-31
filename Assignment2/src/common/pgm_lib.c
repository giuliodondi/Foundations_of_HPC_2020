#include <pgm.h>

#include <stdlib.h>
#include <stdio.h> 
#include <string.h>



#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif



pgm new_pgm() {
	pgm new;
	new.size[0]=0;
	new.size[1]=0;
	new.maxval=0;
	new.data=NULL;
	
	return new;	
}

char allocate_pgm_memory(pgm* image) {
	//allocate the image memory, the threads will touch it later
	unsigned int img_size =  image->size[0]*image->size[1]*image->pix_bytes;
	
	image->data = (uint8_t*)malloc( img_size*sizeof(uint8_t) );
	if (! image->data) {
		printf("Error during image memory allocation.\n");
		return -1;
	}	
	return 0;
}


void copy_pgm( pgm *image1, pgm* image2) {
	
	image2->size[0] = image1->size[0];
	image2->size[1] = image1->size[1];
	image2->maxval = image1->maxval;
	image2->pix_bytes = image1->pix_bytes;
	
	if (image2->data) {
		free(image2->data);
	}
	
	unsigned int elements = image1->size[0] * image1->size[1];
	unsigned int img_size =  elements*image1->pix_bytes;
	
	image2->data = (uint8_t*)malloc( img_size*sizeof(uint8_t) );
	
	memcpy( image2->data , image1->data, img_size*sizeof(uint8_t)  );

}


void clear_pgm( pgm *image) {
	
	image->size[0]=0;
	image->size[1]=0;
	image->maxval = 0;
	image->pix_bytes = 0;
	
	if (image->data) {
		free(image->data);
		image->data=NULL;
	}

}

char write_pgm( pgm *output_img, const char *image_name) {
	
	unsigned int elements = output_img->size[0] * output_img->size[1];
	if ( ( output_img->maxval > 255 ) && I_M_LITTLE_ENDIAN) {
		for (unsigned int i=0; i<elements; ++i) {
			
			((uint16_t*) output_img->data)[i] = swap(  ((uint16_t*) output_img->data)[i] );
			
		}
		
	}
	
	FILE* image_file; 
	image_file = fopen(image_name, "w"); 
	fprintf(image_file, "P5\n# generated by\n# Giulio Dondi\n%d %d\n%d\n", output_img->size[0], output_img->size[1], output_img->maxval);
	
	unsigned int img_size =  elements*output_img->pix_bytes;
	

	if (! fwrite( output_img->data , 1, img_size, image_file) ) {
		printf("Error writing image.\n");
		clear_pgm(output_img);
		fclose(image_file);
		return -1;
	}


	fclose(image_file); 
	return 0;

}




char read_pgm_header( pgm* input_img, const char *image_name, long int* header_offs) {
	
	FILE* image_file; 
	image_file = fopen(image_name, "r"); 

	char    MagicN[2];
	char   *line = NULL;
	size_t  k, n = 0;
  
	// get the Magic Number
	k = fscanf(image_file, "%2s%*c", MagicN );
	
	if (strcmp(MagicN, "P5")) {
		printf("Error: the image is not PGM.\n");
		return -1;
	}

 	
	// skip all the comments
	k = getline( &line, &n, image_file);
	while ( (k > 0) && (line[0]=='#') ) {
    	k = getline( &line, &n, image_file);
	}	
	if (k<=0) {
	  	printf("Error while reading the image header.\n");
	  	free( line );
		fclose(image_file);
	  	return -1;
	}
	
	k = sscanf(line, "%d%*c%d%*c%d%*c", &input_img->size[0], &input_img->size[1], &input_img->maxval);
	if (k<=0) {
	  	printf("Error while reading the image header.\n");
	  	free( line );
		fclose(image_file);
	  	return -1;
	}
	else if ( k < 3 ) {
		//if k=2 the line contained only widt hand height of the image and the 
		//max brightness is on the next line
		k = fscanf(image_file, "%d%*c", &input_img->maxval);
	} 
	free( line );
	
	input_img->pix_bytes = 1 + (input_img->maxval > 255);
	*header_offs = ftell(image_file);
	fclose(image_file);
	return 0;
	
}

char read_pgm_data( pgm* input_img, const char *image_name, long int* skip_header_offs) {
		
	
	FILE* image_file; 
	image_file = fopen(image_name, "r"); 
	
	fseek(image_file, *skip_header_offs*sizeof(uint8_t), SEEK_SET);
	
	unsigned int img_size =  input_img->size[0] * input_img->size[1]*input_img->pix_bytes;
	
	
	if ( fread( input_img->data , 1, img_size, image_file) != img_size )     {
		printf("Error reading image.\n");
		clear_pgm(input_img);
		fclose(image_file);
		return -1;
    }  

	
	endian_swap(input_img);

	
	fclose(image_file);
	return 0;
}


char write_pgm_header( pgm *output_img, const char *image_name, long int* header_offs) {
	
	
	FILE* image_file; 
	image_file = fopen(image_name, "w"); 
	
	fprintf(image_file, "P5\n# generated by\n# Giulio Dondi\n%d %d\n%d\n", output_img->size[0], output_img->size[1], output_img->maxval);
	
	
	*header_offs = ftell(image_file);

	fclose(image_file); 
	return 0;
}


char write_pgm_data( pgm *output_img, const char *image_name) {
	
	FILE* image_file; 
	image_file = fopen(image_name, "a"); 
	
	unsigned int elements = output_img->size[0] * output_img->size[1];
	unsigned int img_size =  elements*output_img->pix_bytes;
	
	endian_swap(output_img);
	

	if (! fwrite( output_img->data , 1, img_size, image_file) ) {
		printf("Error writing image.\n");
		fclose(image_file);
		return -1;
	}
	
	fclose(image_file); 
	return 0;
}


//take care of endianness
void endian_swap(pgm *image) {
	
	unsigned int elements = image->size[0] * image->size[1];
	
	if ( ( image->maxval > 255 ) && I_M_LITTLE_ENDIAN) {
		for (unsigned int i=0; i<elements; ++i) {
			((uint16_t*) image->data)[i] = swap(  ((uint16_t*) image->data)[i] );
		}
	}
	
}


//converts a row-column pair of indices into a linear array index for the image
int img_idx_convert(pgm* image, unsigned int* idx_arr) {
	return  (image->size[0]*idx_arr[1]+ idx_arr[0])*image->pix_bytes;
}



