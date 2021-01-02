#include <kernel_t.h>
#include <stdint.h>


inline int max(int a, int b) { return((a) > (b) ? a : b); }
inline int min(int a, int b) { return((a) > (b) ? b : a); }

void print_usage(char **argv) ;
char is_number( char arg[] ) ;

int8_t read_params_initialise_kernel( int argc, char **argv , char* infile, char* outfile , kernel_t* k );