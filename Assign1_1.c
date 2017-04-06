#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf("To run the program you need to input number of processors and number of elements");
		return -1;
	}
	int N, p;
	p = atoi(argv[1]);
	N = atoi(argv[2]);
	
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
}

