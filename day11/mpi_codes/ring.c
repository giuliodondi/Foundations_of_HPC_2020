#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])        
{
#define MSGLEN 2048
	int idest, isrc, nprocs, rectag;
	unsigned int irank;
	float rmsg[MSGLEN];
	MPI_Status recv_status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);  
	
	//initialise a sample message
	for (int i = 1; i <= MSGLEN; i++)
	{
	  rmsg[i] = 100;
	}
	
	//the message tag will be the rank of the sending process
	//we need to calculate the correct sending and destination ranks
	

	if (irank==0) {
		
		rectag = nprocs - 1;
		isrc   = nprocs - 1;
		idest = irank + 1;
		
	}
	else if ( irank == (nprocs - 1) ) {
		rectag = irank - 1;
		isrc   = irank - 1;
		idest = 0;
	}
	else {
		rectag = irank - 1;
		isrc   = irank - 1;
		idest = irank + 1;
	}
	
	


	if ( irank == 0 )
	{
		printf("Proc %d has sent the message %d to Proc %d\n", irank, irank, idest);
		MPI_Ssend(&rmsg, MSGLEN, MPI_FLOAT, idest, irank, MPI_COMM_WORLD); 
		MPI_Recv(&rmsg, MSGLEN, MPI_FLOAT, isrc, isrc, MPI_COMM_WORLD, &recv_status);
		printf("Proc %d has received the message %d from Proc %d\n", irank, isrc, isrc);

	}
	else
	{
		MPI_Recv(&rmsg, MSGLEN, MPI_FLOAT, isrc, isrc, MPI_COMM_WORLD, &recv_status);
		printf("Proc %d has received the message %d from Proc %d\n", irank, isrc, isrc);
		printf("Proc %d has sent the message %d to Proc %d\n", irank, irank, idest);
		MPI_Ssend(&rmsg, MSGLEN, MPI_FLOAT, idest, irank, MPI_COMM_WORLD); 
	}


	MPI_Finalize();
}
