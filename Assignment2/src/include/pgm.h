#include <stdint.h>

#ifndef PGM_H
#define PGM_H


//pgm struct definition
typedef struct {
	int width;
	int height;
	int maxval;
	uint8_t* data;
	
} pgm;


char write_pgm( pgm *output_img, const char *image_name) ;
char read_pgm( pgm* input_img, const char *image_name) ;
void clear_pgm( pgm* image) ;
void copy_pgm( pgm *image1, pgm* image2) ;

#endif