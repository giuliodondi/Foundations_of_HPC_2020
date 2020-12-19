#include <pgm.h>
#include <kernel_t.h>


void pgm_blur_copy(  pgm* input_img , kernel_t k);
void pgm_blur_linebuf(  pgm* input_img , kernel_t k);