/////////////////////////////////////////////////////////////////////
// CENG316 Parallel Programming
// Programing exam 1
//
// Author: Group 2
// Group members:  Mert Kara 16050151001 , Burcu Özşahin 15050111026
// Date: 24/03/2019
// Description: implementation of myMPI_Allreduce
/////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h> 

#define NUMBER_OF_TESTS 1000
#define ARRAY_SIZE 100000

void myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int main(int argc, char **argv )
{
    int          rank,  numproc;
    double       t1, t2;
    int		     d_in[ARRAY_SIZE];
    int 		 d_out_sum[ARRAY_SIZE], d_out_max[ARRAY_SIZE],d_out_min[ARRAY_SIZE];
    int          i, j, k, nloop;
	
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numproc );
    srand(rank+1); 
    
    if (rank == 0) 
		printf( "Test launched with %d processor, ARRAY_SIZE: %d \n",numproc, ARRAY_SIZE );
    
    for(i=0;i<ARRAY_SIZE;i++)
		d_in[i] = rand()% 100;
		
	MPI_Barrier(MPI_COMM_WORLD );			
	t1 = MPI_Wtime();		
    for (k=0; k<NUMBER_OF_TESTS; k++) {		
		MPI_Allreduce( d_in, d_out_sum, ARRAY_SIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( d_in, d_out_max, ARRAY_SIZE, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
		MPI_Allreduce( d_in, d_out_min, ARRAY_SIZE, MPI_INT, MPI_MIN, MPI_COMM_WORLD );		
    }
    t2 = MPI_Wtime() - t1;		
    	
    MPI_Barrier(MPI_COMM_WORLD);
    int avg=0;
	for(i=0;i<ARRAY_SIZE;i++){ avg += d_out_sum[i]; } avg = avg / (numproc*ARRAY_SIZE);
	printf("Proc %2d: Average of all elements: %d \n", rank, avg );
	
    if (rank == 0) {		
		printf( "\nAverage time for Allreduce: %.5f msec\n", t2/NUMBER_OF_TESTS );
		printf( "First 3-elements of sum array: %d %d %d \n", d_out_sum[0], d_out_sum[1],d_out_sum[2]);		
		printf( "First 3-elements of max array: %d %d %d \n", d_out_max[0], d_out_max[1], d_out_max[2] );		
		printf( "First 3-elements of min array: %d %d %d \n", d_out_min[0], d_out_min[1], d_out_min[2] );		
		printf( "Range of First 3-elements: \t%d %d %d \n", d_out_max[0]-d_out_min[0], d_out_max[1]-d_out_min[1], d_out_max[2]-d_out_min[2]  );
		printf( "Average of First 3-elements: \t%d %d %d \n\n",  d_out_sum[0]/numproc, d_out_sum[1]/numproc,d_out_sum[2]/numproc );		
    }
    
    // reinitializing output arrays
    for(i=0;i<ARRAY_SIZE;i++){ d_out_sum[i] = 0; d_out_max[i]=0; d_out_min[i]=0;  } 
    
    ///////////////////////////////////////////////////////////////////
    /* Calling your implementation of AllReduce */
    /* myMPI_Allreduce */
    
    MPI_Barrier(MPI_COMM_WORLD );			
	t1 = MPI_Wtime();
    for (k=0; k<NUMBER_OF_TESTS; k++) {				
		myMPI_Allreduce( d_in, d_out_sum, ARRAY_SIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		myMPI_Allreduce( d_in, d_out_max, ARRAY_SIZE, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
		myMPI_Allreduce( d_in, d_out_min, ARRAY_SIZE, MPI_INT, MPI_MIN, MPI_COMM_WORLD );		
    }
    t2 = MPI_Wtime() - t1;	
    
	MPI_Barrier(MPI_COMM_WORLD);
    avg=0;
	for(i=0;i<ARRAY_SIZE;i++){ avg += d_out_sum[i]; } avg = avg / (numproc*ARRAY_SIZE);
	printf("Proc %2d: Average of all elements: %d \n", rank, avg );	
    		
    if (rank == 0) {
		printf( "\nAverage time for myAllreduce: %.5f msec\n", t2/NUMBER_OF_TESTS );		
		printf( "First 3-elements of sum array: %d %d %d \n", d_out_sum[0], d_out_sum[1],d_out_sum[2]);		
		printf( "First 3-elements of max array: %d %d %d \n", d_out_max[0], d_out_max[1], d_out_max[2] );		
		printf( "First 3-elements of min array: %d %d %d \n", d_out_min[0], d_out_min[1], d_out_min[2] );		
		printf( "Range of First 3-elements: \t%d %d %d \n", d_out_max[0]-d_out_min[0], d_out_max[1]-d_out_min[1], d_out_max[2]-d_out_min[2]  );
		printf( "Average of First 3-elements: \t%d %d %d \n\n",  d_out_sum[0]/numproc, d_out_sum[1]/numproc,d_out_sum[2]/numproc );
    }

    MPI_Finalize( );
    return 0;
}

void myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
// write your code here
	int rank;
	int numproc;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numproc );
    int tempsum[ARRAY_SIZE];
    int tempmax[numproc][ARRAY_SIZE];
    int tempmin[numproc][ARRAY_SIZE];
    
	if(op==MPI_SUM)
	{	
		int i,k;
		if(rank!=0)
		{
			MPI_Send(sendbuf,ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else{
			for(k=0;k<ARRAY_SIZE;k++)
				{
					((int *)recvbuf)[k] = ((int *)sendbuf)[k];
				}
			for(i=1;i<numproc;i++)
			{
				MPI_Recv(tempsum,ARRAY_SIZE,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for(k=0;k<ARRAY_SIZE;k++)
				{
					((int *)recvbuf)[k] = tempsum[k] +  ((int *)recvbuf)[k];
				}
			}
		}
		if(rank==0)
		{
			int z;
			for(z=1;z<numproc;z++)
			{
				MPI_Send(recvbuf,ARRAY_SIZE, MPI_INT, z, 0, MPI_COMM_WORLD);
			}		
		}
		else
		{
			MPI_Recv(recvbuf,ARRAY_SIZE,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	
	
	if(op==MPI_MAX)
	{
		int i,k,x,y;
		if(rank!=0)
		{
			MPI_Send(sendbuf,ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			for(x=0;x<ARRAY_SIZE;x++)
			{
				tempmax[0][x] = ((int *)sendbuf)[x];
			}
			
			for(i=1;i<numproc;i++)
			{
				MPI_Recv(tempmax[i],ARRAY_SIZE,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			
			for(y=0;y<ARRAY_SIZE;y++)
			{
				
				int max=0;
				for(k=0;k<numproc;k++)
				{
					if(max<tempmax[k][y])
					{
						max=tempmax[k][y];
					}
				}
				((int *)recvbuf)[y] = max;
			}
			
		}
	}
	
	
	if(op==MPI_MIN)
	{
		int i,k,x,y;
		if(rank!=0)
		{
			MPI_Send(sendbuf,ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			for(x=0;x<ARRAY_SIZE;x++)
			{
				tempmin[0][x] = ((int *)sendbuf)[x];
			}
			
			for(i=1;i<numproc;i++)
			{
				MPI_Recv(tempmin[i],ARRAY_SIZE,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			
			for(y=0;y<ARRAY_SIZE;y++)
			{
				
				int min=200;
				for(k=0;k<numproc;k++)
				{
					if(min>tempmin[k][y])
					{
						min=tempmin[k][y];
					}
				}
				((int *)recvbuf)[y] = min;
			}
			
		}
	}
	
}
