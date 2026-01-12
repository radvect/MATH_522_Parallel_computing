#include "utils.h"
#include <random>
#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;

double* fill_random_array(int n){
	double* x = (double*) malloc(n * sizeof(double));	
	srand(time(0));
	for(int i =0; i<n; i++){
		x[i] = rand()%100000000000; 
	}
 	return x;
}

double compute_inner_nonparallel(long N, double * x, double *y){
  double z = 0;

  for(int i =0; i<N; i++){
     z += sqrt(x[i]*y[i]) ;
  }
 // cout<<"Nonparallel algorithm result "<<z<<endl;	
  return z;
}
double inner1(long n, double* v, double* w) {
  int p = omp_get_max_threads();
  double* partial_sum = (double*) malloc(p * sizeof(double));

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    double local_prod = 0;
    #pragma omp for
    for (long i = 0; i < n; i++) {
      local_prod += sqrt(v[i] * w[i]);
    }
    partial_sum[tid] = local_prod;
  }

  double sum = 0;
  for (int i = 0; i < p; i++) sum += partial_sum[i];
  free(partial_sum);

  return sum;
}



double compute_inner_reduction(long N, double* x, double* y) {
  
  double z;

  #pragma omp parallel for reduction(+:z)
  for (long i = 0; i < N; i++){ 
	  z +=sqrt(x[i] * y[i]);}
  return z;
}




double compute_inner_parallel(long N, double * x, double *y){
  int num_threads = omp_get_max_threads();  
  double* sum_per_thread =  (double*) malloc(num_threads* sizeof(double));
  
  #pragma omp parallel
  {
  int i_thread = omp_get_thread_num();	  
  double i_sum = 0; 
  #pragma omp for
  for(long i =0; i<N; i++){
	  i_sum += sqrt(x[i]*y[i]);
  }

  sum_per_thread[i_thread] = i_sum; 
  }  
  double z =0;
  for(int i = 0; i< num_threads; i++)z += sum_per_thread[i];
  
  free(sum_per_thread);
  return z;
}

int main(int argc, char** argv) {
  
  long N = 100000000; 

//double* x = fill_random_array(N);
//double* y = fill_random_array(N);
	double* x = (double*) malloc(N * sizeof(double));
	double* y = (double*) malloc(N * sizeof(double));
for (long i = 0; i < N; i++) {
	x[i] = i+1;
	y[i] = 2.0 / (i+1);
}
  double sum;
  Timer t;
  t.tic();  
  


  sum = compute_inner_nonparallel(N, x,y);

  
  printf("res nonparallel - malloc: %f time elapsed = %f\n",sum, t.toc());

//  t.tic();

//  sum = inner1(N, x,y);

// printf("res parallel - inner1: %f time elapsed = %f\n",sum, t.toc());



  t.tic();
  
  sum = compute_inner_parallel(N, x,y);
  
  printf("res parallel - malloc: %f time elapsed = %f\n",sum, t.toc());
  t.tic();
  sum = compute_inner_reduction(N, x, y); 

  printf("res reduction: %f time elapsed = %f\n",sum, t.toc());


  free(x);
  free(y);
  return 0;
}
 
