#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>

#define SEED 35791246

int cmp(const void* a, const void* b) {
	
	int x=*(int*)a;
	int y=*(int*)b;
	
	return (  x - y );
	
}

int binary_search(int* data, long long int N, int Key) {
	int register low=0;
	int register high=N-1;
	int register mid;

	
	int count = 1;
	
	while(low <= high) {
		mid = (int)(low + high)/2;
		if (data[mid] < Key) {
			low = mid + 1;
		}
		else if (data[mid] > Key) {
			high = mid - 1;	
		}
		else {
			return mid;
		}
		//printf("itern : %d\n",count);
		//++count;
	}
	return -1;
}




int main ( int argc , char *argv[ ] )
{
	
	if ( argc <=1) {
		fprintf (stderr , " Usage :  %s size_of_array \n", argv[0] ) ;
		return -1;
	}
	
	long long int N = atoll(argv[1]);
	double start_time, end_time;  
	
	start_time = clock(); 
	
	//initialise array
	int* arr = (int *)malloc(sizeof(int)*N);
	srand(SEED);
	for(int i=0; i<N; ++i)
        arr[i] = rand()%N;
		//arr[i] = i;
	qsort(arr, N, sizeof(*arr), cmp);
	
	end_time =  clock(); 
	printf("Array initialised and sorted in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC );
	
	//print the entire array
	//for (int i = 0 ; i < N ; i++)
    //    printf ("%d : %d\n ", i, arr[i]);
	
	 
	
	
	
	//element to search for
	int elem;
	
	//read from stdin
	//printf("Enter the value to serch:");
    //scanf("%d", &elem);
	
	//pick a random element from the array itself
	// Use current time as seed for random generator
    srand(time(0));
	//elem = arr[ rand()%N ];
	elem = rand()%N;
	
	start_time = clock(); 
	
	int index = binary_search( arr, N, elem  );
	
	end_time =  clock(); 
	printf("Serial algorithm completed in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC );
	
	if (index==-1) {
		printf ("Element %d is not in the array.\n ", elem);
	} else {
		printf ("Element %d is at index %d.\n ", elem, index);
	}
	
	#if defined(_OPENMP)
	
		start_time = omp_get_wtime(); 
	
		#pragma omp parallel  shared(N, arr, index)
		{
			
			double local_start = omp_get_wtime(); 

			// initialize random numbers 
			int myid =  omp_get_thread_num();
			printf( "Hello from thread %d.\n", myid);

			int p = omp_get_num_threads();
			
			long long int local_N = floor((double)N/p);
			long long int initial_index;
		
			
			if (myid==0) {
				initial_index = local_N*(p-1);
				local_N = N - initial_index;
			} else {
				initial_index = local_N*(myid-1);
			}
			
			
			int local_i = binary_search( &arr[ initial_index ], local_N, elem  );
			
			if (!(local_i==-1)) {
				index = initial_index + local_i;	
			}
			
			double local_end = omp_get_wtime(); 
			printf("Thread %d completed work in %f seconds\n",myid, local_end - local_start);
		}
	

		end_time =  omp_get_wtime(); 
		printf("Parallel algorithm completed in %f seconds\n", (end_time - start_time) );
	
		if (index==-1) {
			printf ("Element %d is not in the array.\n ", elem);
		} else {
			printf ("Element %d is at index %d.\n ", elem, index);
		}
	
	#endif
	
	
	return 0;
	
}