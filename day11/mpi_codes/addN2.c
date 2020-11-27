#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

int sum_int_array(int* array, int len) {
	int summ = 0;
	for (int i=0;i<len; ++i) {
		summ += array[i];
	}	
	
	return summ;
}

//initialise the numbers array
		//or get them from file TODO
		int init_arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

int main(int argc, char *argv[])        
{
	
	int nprocs,n_num;
	unsigned int irank;
		
	MPI_Status recv_status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);  
	

	 // master process
  	int master = 0;
	int numtag = 111;
  	int sendtag = 123;
	
	
	int* num_array;
	
	
	if (irank == master) {
		n_num = sizeof(init_arr)/sizeof(*init_arr);
		
	
		//distribute the numbers
		int numxproc = ceil( (double) n_num/nprocs );
		
		// we need to broadcast numxproc from the master to the children since it 
		//indicates the size of the coming number array
        MPI_Bcast(&numxproc, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		for (int i=1; i<nprocs; ++i) {
			//send to each children process its share of numbers to sum up, starting from index 0
			MPI_Ssend(&init_arr[(i-1)*numxproc] , numxproc, MPI_INT, i, sendtag, MPI_COMM_WORLD); 
			printf("Master sent the numbers to process %d \n", i);

		}
		
		//the last elements belong to the Master processor, populare the num_array with them
		n_num = n_num  - numxproc*(nprocs-1);
		//int num_array[n_num];
		num_array = malloc(sizeof(int)*n_num);
		memcpy(num_array, &init_arr[ numxproc*(nprocs-1) ], n_num*sizeof(*init_arr));	
		
		
		//calculate the partial sum
		int summ = sum_int_array(num_array,n_num);
		printf("Process %d's partial sum is : %d\n", irank, summ);
		
		//allocate a sums array to store al the results gathered from the children
		int* sums = malloc(sizeof(int)*nprocs);
		MPI_Gather(&summ, 1, MPI_INT, sums, 1, MPI_INT, 0,MPI_COMM_WORLD);
		printf("Master has received the sums.\n");
		
		//calculate the final result
		summ = sum_int_array(sums,nprocs);
		printf("The final sum is: %d \n", summ);
		free(sums);
		
	}
	else {
		//get the message size
		MPI_Bcast(&n_num, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		//create the array to store the incoming numbers
		num_array = malloc(sizeof(int)*n_num);
		MPI_Recv(num_array, n_num, MPI_INT, master, sendtag, MPI_COMM_WORLD, &recv_status);
		printf("Process %d received the numbers from Master.\n", irank);
		
		//calculate the partial sum
		int summ = sum_int_array(num_array,n_num);
		printf("Process %d's partial sum is : %d\n", irank, summ);
		
		//send back the result
		printf("Process %d is sending back the sum.\n", irank);
		MPI_Gather(&summ, 1, MPI_INT, NULL, 1, MPI_INT, 0,MPI_COMM_WORLD);
		
	}


	
	free(num_array);

	MPI_Finalize();
	return 0;
}
