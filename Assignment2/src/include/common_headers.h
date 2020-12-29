#include <kernel_t.h>
#include <stdint.h>



int min( const int a, const int b);
int max( const int a, const int b);

void print_usage(char **argv) ;
char is_number( char arg[] ) ;

int8_t read_params_initialise_kernel( int argc, char **argv , char* infile, char* outfile , kernel_t* k );

