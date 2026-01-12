#include "mpi.h"
#include <stdio.h>
using namespace std; 
double compute_inner_reduction(long N, double* x, double* y) {
  
  double z;

  #pragma omp parallel for reduction(+:z)
  for (long i = 0; i < N; i++){ 
	  z +=(x[i] * y[i]);}
  return z;
}

int main( int argc, char *argv[]
)
{ 
int rank, size;
MPI_Init( &argc, &argv );
MPI_Comm_rank( MPI_COMM_WORLD, &rank );
MPI_Comm_size( MPI_COMM_WORLD, &size );

long N = 100000000;
len_per_proc = N/size;
len_rest_proc= N%size;
double dotProduct;
if(rank==0){
double* x = (double*) malloc(N * sizeof(double));
double* y = (double*) malloc(N * sizeof(double));
for (long i = 0; i < N; i++) {
        x[i] = i+1;
        y[i] = 2.0 / (i+1);}
}



if(rank == size-1){
len_per_proc ==len_rest_proc;

}

double sum=0;
double *x_i = (double*) malloc(len_per_proc * sizeof(double));
double *y_i = (double*) malloc(len_per_proc * sizeof(double));

    
MPI_Scatter(x, len_per_proc, MPI_DOUBLE, x_i, len_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Scatter(y, len_per_proc, MPI_DOUBLE, y_i, len_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//CODE THAT REALIZES VECTOR RESIDUE ALLOCATION TO THE LAST PROCESSOR.

sum  = compute_inner_reduction(N, x_i, y_i);
MPI_Reduce(&sum, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);


cout<<dotProduct<<endl;

MPI_Finalize();



return 0;
}


