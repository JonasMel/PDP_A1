#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


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
			//c[i*dim + j] = temp_c;
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
	int N, p, myid, ndim, rank, root, sqrt_p;
	int dims[2], coords[2],  mycoords[2], cyclic[2], reorder;
	int i, j, blksqr, blk_size;
	double *A, *B;
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
	reorder = 0;
	sqrt_p = sqrt(p);
	blk_size = N/sqrt_p;
	MPI_Barrier(MPI_COMM_WORLD);
	dims[0] = sqrt_p;
	dims[1] = sqrt_p;
	MPI_Type_vector(blk_size,blk_size, N, MPI_DOUBLE, &blk);
	MPI_Type_commit(&blk);
	
	blksqr = blk_size*blk_size;
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, cyclic, reorder, &comm2D);
	MPI_Comm_rank(comm2D, &myid);
	MPI_Cart_coords(comm2D, myid, ndim, mycoords);
	MPI_Comm_split(comm2D, mycoords[0], mycoords[1], &cart_row);
	MPI_Comm_split(comm2D, mycoords[1], mycoords[0], &cart_col);
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

	double *buffa, *buffb, *buffacpy, *tmp_c;
	buffa = malloc(blksqr*sizeof(double));
	buffb = malloc(blksqr*sizeof(double));
	buffacpy = malloc(blksqr*sizeof(double));
	tmp_c = malloc(blksqr*sizeof(double));
	/*printf("Processor %d. After dist\n", myid);*/
	MPI_Recv(buffa, blksqr, MPI_DOUBLE, root, myid, comm2D, &status[0]);
	/*printf("Processor %d.recieved data to buffa\n", myid);*/
	MPI_Recv(buffb, blksqr, MPI_DOUBLE, root, N+myid, comm2D, &status[1]);
	/*printf("Processor %d.recieved data to buffb\n", myid);*/
	int source, dest;
	source = (mycoords[0]+1) % sqrt_p;
	dest = (mycoords[0] + sqrt_p - 1) % sqrt_p;
	for (int stage = 0; stage < dims[0]; stage++)
	{
			
			coords[0] = mycoords[0];
			coords[1] = (mycoords[0]+stage)%dims[1];
			memcpy(buffacpy, buffa, blksqr*sizeof(double));
			MPI_Cart_rank(comm2D, coords, &rank);
			//printf("rank = %d, (%d, %d), i = %d\n", coords[1], mycoords[0], coords[1], stage);
			MPI_Bcast(buffacpy, blksqr, MPI_DOUBLE, coords[1], cart_row);
			blk_multi(buffacpy, buffb, tmp_c, blk_size);
			/*printf("Processor %d. Sending to %d\n", root, rank);*/
			MPI_Sendrecv_replace(buffb, blksqr, MPI_DOUBLE, dest, 0, \
									source, 0, cart_col, &status[0]);
	}
	MPI_Isend(tmp_c, blksqr, MPI_DOUBLE, root, myid, comm2D, &requestC);
/*	printf("Processor = %d Buffacpy\n", myid);
	for (i = 0; i < blk_size; i++)
	{
		for (j = 0; j < blk_size; j++)
		{
				printf(" %lf", buffacpy[i*blk_size +j]);
		}
		printf("\n");
	}
	printf("\n");*/

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
/*
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
*/
	
	double *C;

	if (myid == root)
	{
		MPI_Waitall(p, requestB, NULL);
		MPI_Waitall(p, requestA, NULL);
		
		C = malloc(N*N*sizeof(double));
		for (i = 0; i < p; i++)
		{
			MPI_Probe(i, i, comm2D, &statusC);
			MPI_Cart_coords(comm2D, i, ndim, coords);
			MPI_Recv(&C[coords[0]*blk_size*N + blk_size*coords[1]], 1, blk, \
						i, i, comm2D, &statusC);
			//printf("C[%d,%d] = %lf\n", coords[0], coords[1], \
					C[coords[0]*blk_size*N + coords[1]*blk_size]);
			
		}

	}
	MPI_Wait(&requestC, &statusC);

	
	if (myid == root)
	{
		free(A);
		free(B);
		free(requestA);
		free(requestB);
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
		free(C);
	}

	free(buffacpy);
	free(tmp_c);
	free(buffa);
	free(buffb);

	MPI_Finalize();


	return 0;
}

