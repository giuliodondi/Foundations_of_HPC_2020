#include <pgm.h>
#include <kernel_t.h>






void pgm_blur_copy(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx2(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx4(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx8(  pgm* input_img , const  kernel_t* k);

//blurring function headers
void pgm_blur_halo_unrolx2(  pgm* input_img , const kernel_t* k,  const int* halos);
void pgm_blur_halo_unrolx4(  pgm* input_img , const kernel_t* k,  const int* halos);

void blur_func_manager(  pgm* input_img , kernel_t* k);
void blur_halo_func_manager(  pgm* input_img , const  kernel_t* k,  const int* halos);