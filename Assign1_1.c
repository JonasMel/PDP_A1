#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define kanske 1
int main(int argc, char **argv)
{
	
	if (argc != 3)
	{
		printf("To run the program you need to input number of processors and number of elements");
		return -1;
	}
	int N, p, myid, ndim, rank, root, sqrt_p;
	int dims[2], coords[2], cyclic[2], reorder;
	int i, j, blksqr, blk_size;
	double *A, *B;
	int row[2], col[2];
	MPI_Comm comm2D, cart_row, cart_col;
	MPI_Datatype blk;
	MPI_Request *requestA, *requestB, request[2];
	MPI_Status status[2];
	p = atoi(argv[1]);
	N = atoi(argv[2]);
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	
	row[0] = 0;
	row[1] = 1;
	col[0] = 1;
	col[1] = 0;
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
	
	blksqr = blk_size*blk_size;
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, cyclic, reorder, &comm2D);
	MPI_Comm_rank(comm2D, &myid);
	MPI_Cart_coords(comm2D, myid, ndim, coords);
	MPI_Comm_split(comm2D, coords[0], coords[1], &cart_row);
	/*printf("Processor %d.\n", myid);*/
	if (myid == root)
	{


		

		/* Allocating memory for matrices.*/
		
		A = malloc(N*N*sizeof(double));
		B = malloc(N*N*sizeof(double));
		requestA = malloc(p*sizeof(MPI_Request));
		requestB = malloc(p*sizeof(MPI_Request));		
		/* initialize data */
		printf("A\n");
		for (i=0; i< N; i++) 
		{
			for (j=0; j<N; j++) 
			{
				
				A[i*N+j] = (double)(rand())/RAND_MAX;
				B[i*N+j] = (double)(rand())/RAND_MAX;
				printf(" %lf",A[i*N +j]);
			}
			printf("\n");
		}
		printf("B\n");
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
					printf(" %lf", B[i*N +j]);
			}
			printf("\n");
		}
	  MPI_Type_vector(blk_size,blk_size, N, MPI_DOUBLE, &blk);
	  MPI_Type_commit(&blk);

	  for (i = 0; i < dims[0]; i++)
	  {
		  for (j = 0; j < dims[1]; j++)
		  {

			  coords[0] = i;
			  coords[1] = j;
			  MPI_Cart_rank(comm2D, coords, &rank);
			  MPI_Isend(&A[i*N*blk_size + j*blk_size], 1, blk, rank, rank, comm2D \
			  , &requestA[i*sqrt_p +j]);
			  MPI_Isend(&B[i*N*blk_size + j*blk_size], 1, blk, rank, rank+N, comm2D \
			  , &requestB[i*sqrt_p +j]);
/*			  printf("Processor %d. Sending to %d\n", root, rank);*/
		  }
	  }
	


	}

	double *buffa, *buffb;
	buffa = malloc(blk_size*blk_size*sizeof(double));
	buffb = malloc(blk_size*blk_size*sizeof(double));
	/*printf("Processor %d. After dist\n", myid);*/
	MPI_Recv(buffa, blksqr, MPI_DOUBLE, root, myid, comm2D, &status[0]);
	/*printf("Processor %d.recieved data to buffa\n", myid);*/
	MPI_Recv(buffb, blksqr, MPI_DOUBLE, root, N+myid, comm2D, &status[1]);
	/*printf("Processor %d.recieved data to buffb\n", myid);*/
	int rowID;
		for (i = 0; i < dims[0]; i++)
		{
			coords[0] = i;
			coords[1] = (i+1)%dims[1];
			MPI_Cart_rank(comm2D, coords, &rank);
			if (myid == rank) MPI_Bcast(buffa, 1, blk, rank, cart_row);
			/*printf("Processor %d. Sending to %d\n", root, rank);*/
		}
/*	printf("Processor = %d Buffa\n", myid);
	for (i = 0; i < blk_size; i++)
	{
		for (j = 0; j < blk_size; j++)
		{
				printf(" %lf", buffa[i*blk_size +j]);
		}
		printf("\n");
	}
	printf("Processor = %d Buffb\n", myid);
	for (i = 0; i < blk_size; i++)
	{
		for (j = 0; j < blk_size; j++)
		{
				printf(" %lf", buffb[i*blk_size +j]);
		}
		printf("\n");
	}*/

	if (myid == root)
	{

//		MPI_Waitall(p, requestA, NULL);
		MPI_Waitall(p, requestB, NULL);
		free(A);
		free(B);
		free(requestA);
		free(requestB);

	}
	
	free(buffa);
	free(buffb);
	MPI_Finalize();
	return 0;
}

