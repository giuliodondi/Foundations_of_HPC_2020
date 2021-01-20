#include <pgm.h>
#include <kernel_t.h>




//blur functions without haloing
void pgm_blur_copy(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx2(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx4(  pgm* input_img , const  kernel_t* k);
void pgm_blur_linebuf_unrolx8(  pgm* input_img , const  kernel_t* k);


//haloed blur functions
void pgm_blur_halo_unrolx2(  pgm* input_img , const kernel_t* k,  const int* halos);
void pgm_blur_halo_unrolx4(  pgm* input_img , const kernel_t* k,  const int* halos);
void pgm_blur_halo_unrolx8(  pgm* input_img , const kernel_t* k,  const int* halos);

