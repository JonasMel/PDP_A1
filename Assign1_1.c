#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	
	if (argc != 3)
	{
		printf("To run the program you need to input number of processors and number of elements");
		return -1;
	}
	int N, p, myid, ndim, rank, root, sqrt_p;
	int dims[2], coords[2], cyclic[2], reorder;
	MPI_Comm comm2D;
	MPI_Datatype blk
	MPI_Request request;
	MPI_Status status;
	p = atoi(argv[1]);
	N = atoi(argv[2]);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	ndim = 2;
	root = 0;
	cyclic[0] = 0;
	cyclic[1] = 0;
	reorder = 0;
	sqrt_p = sqrt(p);
	blk_size = N/sqrt_p;
	MPI_Barrier(MPI_COMM_WORLD);
	dims[0] = sqrt_p;
	dims[1] = sqrt_p;
	
	MPI_Cart_Create(MPI_COMM_WORLD, ndim, dims, cyclic, reorder &comm2D);
	
	
	if (myid == root)
	{


		int i, j;
		double A, B;
		/* Allocating memory for matrices.*/
		A = malloc(N*sizeof(double));
		B = malloc(N*sizeof(double));
		for (i = 0; i < N; i++)
		{
			A[i] = malloc(N*sizeof(double));
			B[i] = malloc(N*sizeof(double));
		}
		
		/* initialize data */
	    for (i=0; i< N; i++) 
	    {
	        for (j=0; j<N; j++) 
	        {
	            A[i*N+j] = double(rand())/RAND_MAX;
	            B[i*N+j] = double(rand())/RAND_MAX;
	        }
	    }
	  MPI_Type_vector(blk_size,blk_size, N, MPI_DOUBLE, &blk);
	  MPI_Type_commit(&blk);
	  int blksqr; 
	  blksqr = blk_size*blk_size;
	  for (i = 0; i < sqrt_p; i++)
	  {
		  for (j = 0; j < sqrt_p; j++)
		  {
			  MPI_Cart_rank(comm2D, coords, &rank)
			  MPI_Isend(A[i*N + blk_size*j], blksqr, blk, rank, rank, comm2D \
			  , &request);
			  MPI_Isend(B[i*N + blk_size*j], blksqr, blk, rank, rank+N, comm2D \
			  , &request);
		  }
	  }
    }
    
}

