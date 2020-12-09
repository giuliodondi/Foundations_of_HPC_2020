#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#define SEED 35791246

int main ( int argc , char *argv[ ] )
{
  

	// number of points inside the circle
	long long int N, M=0 ; 
	double pi ;
	
	// times 
	double start_time, end_time;   
	int nthreads=1;
	
	
	if ( argc <=1) {
		fprintf (stderr , " Usage :  %s number_of_iterations \n", argv[0] ) ;
		return -1;
	}
	
	N = atoll(argv[1]);
	

	
	#pragma omp parallel shared(M)
	{
		
		start_time = omp_get_wtime(); 

		// initialize random numbers 
		int myid =  omp_get_thread_num();
		printf( "Hello from thread %d.\n", myid);

		
		unsigned int myseed = (SEED*(myid+1)) ; // seed the number generator
		
		// coordinates
		double x, y ;

		
		#pragma omp for  reduction (+:M)
		for (long long i=0; i<N ; i++) {
			
			// take a point P(x,y) inside the unit square
			x = ((double)rand_r(&myseed))/RAND_MAX ;
			y = ((double)rand_r(&myseed))/RAND_MAX;

			// check if the point P(x,y) is inside the circle
			if ((x*x + y*y)<1) {				
				//#pragma omp atomic
				M++;
			}
		}
		
		
		end_time = omp_get_wtime(); 
	
		printf("Thread %d completed work in %f seconds\n",myid, end_time - start_time);
		
	}

	
	pi = 4.0*M/N ;
	printf ( "\nTotal # of throws = %llu , estimate of pi is %1.9f \n", N*nthreads, pi ) ;

	
	return 0;
}