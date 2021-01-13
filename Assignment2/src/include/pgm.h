#include <stdint.h>

#ifndef PGM_H
#define PGM_H


//pgm struct definition
typedef struct {
	int maxval;
	int size[2];
	uint8_t pix_bytes;
	uint8_t* data;
} pgm;


pgm new_pgm();
char allocate_pgm_memory( pgm* image) ;
void clear_pgm( pgm* image) ;
void copy_pgm( pgm *image1, pgm* image2) ;

char read_pgm_header( pgm* input_img, const char *image_name, long int* header_offs);
char read_pgm_data( pgm* input_img, const char *image_name, long int* skip_header_offs) ;
char write_pgm_header( pgm* input_img, const char *image_name, long int* header_offs);
char write_pgm_data( pgm* input_img, const char *image_name) ;
void endian_swap(pgm *image) ;
int img_idx_convert(pgm* image, int* idx_arr) ;

void compare_pgm( pgm *image1, pgm* image2, const char *outfile) ;

#endif