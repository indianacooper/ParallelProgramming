// TODO: Make serial gauss() run
#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
/* Program Parameters */

#define MAXN 2000
/* Max value of N */
int N;
/* Matrix size */
int procs;
/* Number of processors to use */

/* Matrices and vectors */
volatile
float A[MAXN][MAXN], B[MAXN],
 X[MAXN];
/* A * X = B, solve for X */
/* Prototype */
void gauss();
/* The function you will provide.
* It is this routine that is timed.
* It is called only on the parent.
*/
int main3()
{
    /* static matrix A B
    1 -2 9 | 8
    3 1 -1 | 3
    2 -8 1 | -5
    */

    A[0][0]= 1;
    A[0][1]= -2;
    A[0][2]= 9.0;
    A[1][0]= 3.0;
    A[1][1]= 1.0;
    A[1][2]= -1.0;
    A[2][0]= 2.0;
    A[2][1]= -8.0;
    A[2][2]= 1.0;
    B[0]= 8.0;
    B[1]= 3.0;
    B[2]= -5.0;
    N = 3;
    /* Gaussian Elimination */
    gauss();
    /* Display output */
    int i;
    for ( i = 0; i < 3;i++)
    {
        printf("%f\n",X[i]);
    }
    printf("End of main\n");
    return 0;
}

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
* defined in the beginning of this code. X[] is initialized to zeros.
*/
void gauss()
{
    int norm, row, col;
    /* Normalization row, and zeroing
    * element row and col */
    float multiplier;
    printf("Computing Serially.\n");
    /* Gaussian elimination */
    for (norm = 0; norm < N - 1; norm++)
    {
        for (row = norm + 1; row < N; row++)
        {
            multiplier = A[row][norm] / A[norm][norm];
            for (col = norm; col < N; col++)
            {

                A[row][col] -= A[norm][col] * multiplier;

            }
            B[row] -= B[norm] * multiplier;
        }
    }
/* (Diagonal elements are not normalized to 1. This is treated in back
* substitution.)
*/

/* Back substitution */
    for (row = N - 1; row >= 0; row--)
    {
        X[row] = B[row];
        for (col = N-1; col > row; col--)
        {
            X[row] -= A[row][col] * X[col];
        }
    X[row] /= A[row][row];
    }
}
