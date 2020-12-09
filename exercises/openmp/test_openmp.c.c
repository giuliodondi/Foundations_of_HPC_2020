#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(){
	
	
	int nthreads;
	
	long int N = 10000;
	
	#if defined(_OPENMP)
	
	double start; 
	double end; 
	start = omp_get_wtime(); 
	
	
	
	#pragma omp parallel 
	{
	
		#pragma omp master
    	nthreads = omp_get_num_threads();
		
		//#pragma omp master
		printf(" %d threads will be spawned.\n",nthreads );
		
		int my_thread_id = omp_get_thread_num();
	
		int nn=N*nthreads;
		
		printf( "thread num %d is doing %d amount of work.\n", my_thread_id, nn);
		
		for (int i=0;i<nn;) {
				++i;
		}
		
	
	}
	
	end = omp_get_wtime(); 
	
	printf("Work took %f seconds\n", end - start);
	
	#else

	  nthreads = 1;
	  printf( "\tgreetings from thread num 0\n");
	
	#endif

	  printf(" %d thread(s) greeted you from the parallel region\n",nthreads );

  return 0;
	

}