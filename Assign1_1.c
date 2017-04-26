/*
 * Partial and distributed programming 
 * Assignment 1. April 13 2017
 * Authors:
 * Jonas Melander
 * Lina Viklund
 * Aleksandra Obeso Duque
 * 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#define print 0

/* Function calculating matrix multiplication and adding to previous value*/
void blk_multi(double *a,double *b, double *c, int dim)
{
	int i, j, k;
	double temp_c;
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			temp_c = 0.0;
			for (k = 0; k < dim; k++)
			{
				c[i*dim + j] = c[i*dim + j] + a[i*dim + k] * b[k*dim + j];
			}
		}
	}
	return;
}


int main(int argc, char **argv)
{
	
	if (argc != 3)
	{
		printf("To run the program you need to input number of processors and number of elements");
		return -1;
	}
	/* Initializing constants.*/
	int N, p, myid, ndim, rank, root, sqrt_p;
	int dims[2], coords[2],  mycoords[2], cyclic[2], reorder;
	int i, j, blksqr, blk_size;
	double *A, *B, *C;
	int row[2], col[2];
	MPI_Comm comm2D, cart_row, cart_col;
	MPI_Datatype blk;
	MPI_Request *requestA, *requestB, requestC,request[2];
	MPI_Status status[2], statusC;
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
	reorder = 1;
	sqrt_p = sqrt(p);
	blk_size = N/sqrt_p;
	blksqr = blk_size*blk_size;
	dims[0] = sqrt_p;
	dims[1] = sqrt_p;
	
	/* Creating processor grids.*/
	MPI_Type_vector(blk_size,blk_size, N, MPI_DOUBLE, &blk);
	MPI_Type_commit(&blk);
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, cyclic, reorder, &comm2D);
	MPI_Comm_rank(comm2D, &myid);
	MPI_Cart_coords(comm2D, myid, ndim, mycoords);
	MPI_Comm_split(comm2D, mycoords[0], mycoords[1], &cart_row);
	MPI_Comm_split(comm2D, mycoords[1], mycoords[0], &cart_col);

	if (myid == root)
	{
		/* Allocating memory for matrices.*/
		A = malloc(N*N*sizeof(double));
		B = malloc(N*N*sizeof(double));
		requestA = malloc(p*sizeof(MPI_Request));
		requestB = malloc(p*sizeof(MPI_Request));
		
		/* Initialize data */
		srand(time(NULL));
		for (i=0; i< N; i++) 
		{
			for (j=0; j<N; j++) 
			{
				A[i*N+j] = (double)(rand())/RAND_MAX;
				B[i*N+j] = (double)(rand())/RAND_MAX;
			}
		}
		
		#if print
		{
			printf("\n");			
			printf("A\n");
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
						printf(" %lf", A[i*N +j]);
				}
				printf("\n");
			}
			printf("\n");
			printf("B\n");
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
						printf(" %lf", B[i*N +j]);
				}
				printf("\n");
			}
			printf("\n");		
		}
		#endif
		
		/*Distributing data*/
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
			}
		}
	}

	/* Recieving data in local buffers.*/
	double *buffa, *buffb, *buffacpy, *buffbcpy, *tmp_c;
	buffa = malloc(blksqr*sizeof(double));
	buffb = malloc(blksqr*sizeof(double));
	buffacpy = malloc(blksqr*sizeof(double));
	buffbcpy = malloc(blksqr*sizeof(double));
	tmp_c = malloc(blksqr*sizeof(double));
	MPI_Recv(buffa, blksqr, MPI_DOUBLE, root, myid, comm2D, &status[0]);
	MPI_Recv(buffb, blksqr, MPI_DOUBLE, root, N+myid, comm2D, &status[1]);

	
	/* Calculating shift for matrix B.*/
	int source, dest;
	source = (mycoords[0]+1) % sqrt_p;
	dest = (mycoords[0] + sqrt_p - 1) % sqrt_p;
	
	if (myid ==root) 	printf("Test with %d processors and %d x %d matrices\n", p, N, N);
	for (int loop = 0; loop < 10; loop++){
	double begin, end;
	begin = MPI_Wtime();
	
	/* Fox's algorithm.*/
	for (int stage = 0; stage < dims[0]; stage++)
	{
			/* Calculating position of block in A to broadcast.*/
			coords[0] = mycoords[0];
			coords[1] = (mycoords[0]+stage)%dims[1];
			
			/* Copying local A block for broadcast to avoid overwriting original.*/
			memcpy(buffacpy, buffa, blksqr*sizeof(double));
			memcpy(buffbcpy, buffb, blksqr*sizeof(double));
			/* Broadcasting.*/
			MPI_Bcast(buffacpy, blksqr, MPI_DOUBLE, coords[1], cart_row);
			

			/* Shiftig B.*/
			MPI_Isend(buffbcpy, blksqr, MPI_DOUBLE, dest, dest, cart_col,  &request[0]);
			MPI_Irecv(buffb, blksqr, MPI_DOUBLE, source, mycoords[0], \
			cart_col, &request[1]);
			/*MPI_Sendrecv_replace(buffb, blksqr, MPI_DOUBLE, dest, 0, \
									source, 0, cart_col, &status[0]);*/
			/* Multiplying local copies of A and B placing them in a local C*/
			blk_multi(buffacpy, buffbcpy, tmp_c, blk_size);
			MPI_Wait(&request[0], &status[0]);
			MPI_Wait(&request[1], &status[1]);
	}
	if (myid == root){
	end = MPI_Wtime();
	printf("Time(%d) = %5.16lf\n", loop, end - begin);}
	}
	/* Sending each processors result to root processor.*/
	MPI_Isend(tmp_c, blksqr, MPI_DOUBLE, root, myid, comm2D, &requestC);
	
	
	#if print
	{
		printf("Processor = %d Buffacpy\n", myid);
		for (i = 0; i < blk_size; i++)
		{
			for (j = 0; j < blk_size; j++)
			{
					printf(" %lf", buffacpy[i*blk_size +j]);
			}
			printf("\n");
		}
		printf("\n");

		printf("Processor = %d tmp_c\n", myid);
		for (i = 0; i < blk_size; i++)
		{
			for (j = 0; j < blk_size; j++)
			{
					printf(" %lf", tmp_c[i*blk_size +j]);
			}
			printf("\n");
		}
		printf("\n");

		printf("Processor = %d Buffb\n", myid);
		for (i = 0; i < blk_size; i++)
		{
			for (j = 0; j < blk_size; j++)
			{
					printf(" %lf", buffb[i*blk_size +j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	#endif

	

	
	if (myid == root)
	{	
		MPI_Status *sA, *sB;
		sA = (MPI_Status*)malloc(p*sizeof(MPI_Status));
		sB = (MPI_Status*)malloc(p*sizeof(MPI_Status));
		/* Waiting for nonblocking send to complete.*/
		MPI_Waitall(p, requestB, sB);
		MPI_Waitall(p, requestA, sA);
		
		/* Allocating space for result matrix, C, in root.*/
		C = malloc(N*N*sizeof(double));
		
		/* Probing for results and recieving them*/
		for (i = 0; i < p; i++)
		{
			MPI_Probe(i, i, comm2D, &statusC);
			MPI_Cart_coords(comm2D, i, ndim, coords);
			MPI_Recv(&C[coords[0]*blk_size*N + blk_size*coords[1]], 1, blk, \
						i, i, comm2D, &statusC);
			
		}

	}
	/* Waiting for nonblocking send to complete.*/
	MPI_Wait(&requestC, &statusC);

	/* Freeing memory local to rooot.*/
	if (myid == root)
	{
		free(A);
		free(B);
		free(requestA);
		free(requestB);
		#if print
		{
			printf("\n");
			printf("C\n");
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
						printf(" %lf", C[i*N +j]);
				}
				printf("\n");
			}
		
		}
		#endif
		free(C);


	}
	
	/* Freeing memory local to all processors.*/
	free(buffacpy);
	free(buffbcpy);
	free(tmp_c);
	free(buffa);
	free(buffb);


	MPI_Finalize();
	return 0;
}

