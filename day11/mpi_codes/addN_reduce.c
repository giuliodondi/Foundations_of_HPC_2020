#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#define USE MPI

int sum_int_array(int* array, int len) {
	int summ = 0;
	for (int i=0;i<len; ++i) {
		summ += array[i];
	}	
	
	return summ;
}

int main(int argc, char *argv[])        
{
	
	int nprocs,n_num,i;
	unsigned int irank;
		
	MPI_Status recv_status;
	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);  
	
	if ( argc <=1) {
		fprintf (stderr , "\nToo few arguments. Usage : mpirun -np n %s numbers_file \n", argv[0] ) ;
		MPI_Finalize() ;
		exit(-1) ;
	  }
	

	 // master process
  	int master = 0;
	int numtag = 111;
  	int sendtag = 123;
	
	
	int* num_array;
	
	
	if (irank == master) {
		
		//read the numbers from file 
		int* init_arr;
		FILE * f =fopen(argv[1],"r");
		if (! f ) { 
			printf("\nUnable to open file.\n"); 
			MPI_Finalize();
			return 1;
		}
		else { 
			fscanf(f, "%d\n", &n_num); 
			init_arr = (int *)malloc(sizeof(int)*n_num);

			i=1;
			int readnum;
			while( (fscanf(f, "%d\n", &readnum) == 1) && (i<=n_num) ){
				init_arr[i-1]=readnum;
				++i;
			}

			n_num=i;
			fclose(f);		
		}
		
		//distribute the numbers

		// we need to calculate how many numbers to give each children process
		// and let them know since it indicates the size of the coming number array
		int numxproc = ceil( (double) n_num/nprocs );
		printf("\nEach children process will receive %d numbers.\n", numxproc);
        MPI_Bcast(&numxproc, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		for (i=1; i<nprocs; ++i) {
			//send to each children process its share of numbers to sum up, starting from index 0
			MPI_Ssend(&init_arr[(i-1)*numxproc] , numxproc, MPI_INT, i, sendtag, MPI_COMM_WORLD); 
			printf("Master sent the numbers to process %d \n", i);

		}
		
		//the last elements belong to the Master processor, populare the num_array with them
		n_num = n_num  - numxproc*(nprocs-1);
		num_array = (int *)malloc(sizeof(int)*n_num);
		memcpy(num_array, &init_arr[ numxproc*(nprocs-1) ], n_num*sizeof(*init_arr));	
		//the complete array is no longer needed, deallocate
		free(init_arr);
		
	}
	else {
		//get the message size
		MPI_Bcast(&n_num, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		//create the array to store the incoming numbers
		num_array = (int *)malloc(sizeof(int)*n_num);
		MPI_Recv(num_array, n_num, MPI_INT, master, sendtag, MPI_COMM_WORLD, &recv_status);
		printf("Process %d received the numbers from Master.\n", irank);
		
	}

	
	//calculate the partial sum
	int mysum = sum_int_array(num_array,n_num);
	printf("Process %d's partial sum is : %d\n",irank, mysum);

	int total_sum = 0;
	//collect and reduce all the prtial sums
	MPI_Reduce(&mysum, &total_sum, 1, MPI_INT, MPI_SUM, master, MPI_COMM_WORLD);
		   
	if (irank == master) {
		printf("FINAL RESULT: %d \n", total_sum);	
	}

	
	free(num_array);

	MPI_Finalize();
	return 0;
}
