#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#define N 1024 // Grid size (50x50)
#define TOL 1e-6
#define MAX_ITER 100000
#include <unistd.h>
// Boundary Conditions
#define LEFT_TEMP 1.2
#define RIGHT_TEMP 1.3
#define TOP_TEMP 1.5
#define BOTTOM_TEMP 1.1

#include <sys/resource.h>

void debug_print_grid_i(double **grid_i, int local_size, int processor_id, int num_of_proc) {
    MPI_Barrier(MPI_COMM_WORLD);  

    for (int i = 0; i < num_of_proc; i++) {
        MPI_Barrier(MPI_COMM_WORLD);  

        if (processor_id == i) {
            printf("====== Processor %d printing grid_i ======\n", processor_id);
            for (int row = 0; row < local_size; row++) {
                printf("Proc %d | ", processor_id);
                for (int col = 0; col < N; col++) {
                    printf("%6.2f ", grid_i[row][col]);
                }
                printf("\n");
            }
            printf("========================================\n");
            fflush(stdout);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);  
}

// void initialize_grid( double** grid ) {
//     // Set boundary conditions
//     for (int i = 0; i < N; i++) {
//         grid[i][0] = top_TEMP;
//         grid[i][N-1] = bottom_TEMP;
//         grid[0][i] = TOP_TEMP;
//         grid[N-1][i] = BOTTOM_TEMP;
//     }
//     // Interior points initial guess
//     for (int i = 1; i < N-1; i++) {
//         for (int j = 1; j < N-1; j++) {
//             grid[i][j] = 0.0;
//         }
//     }
// }
void print_grid( double** grid ) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%6.2f ", grid[i][j]);
        }
        printf("\n");
    }
}

double update_red_black(double **grid_i, int local_size, int initial_position) {
    int processor_id;
    int num_of_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
    double* send_temp_grid_top  = (double*)malloc(N * sizeof(double));
    double* send_temp_grid_bottom = (double*)malloc(N * sizeof(double));
    double* recv_temp_grid_top  = (double*)malloc(N * sizeof(double));
    double* recv_temp_grid_bottom = (double*)malloc(N * sizeof(double));






    if (num_of_proc!=1){

    if (processor_id == 0) { 
        for (int j = 0; j < N; j++) {
            send_temp_grid_bottom[j] = grid_i[local_size - 1][j];
            recv_temp_grid_top[j] = TOP_TEMP;
        }
       MPI_Send(send_temp_grid_bottom, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
       MPI_Recv(recv_temp_grid_bottom, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } 
    else if (processor_id == num_of_proc - 1) {
        for (int j = 0; j < N; j++) {
            send_temp_grid_top[j] = grid_i[0][j]; 
            recv_temp_grid_bottom[j] = BOTTOM_TEMP;

        }
        MPI_Recv(recv_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(send_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD);
    } 
    else {
        for (int j = 0; j < N; j++) {
            send_temp_grid_top[j] = grid_i[0][j];             
            send_temp_grid_bottom[j] = grid_i[local_size - 1][j]; 
        }
        MPI_Recv(recv_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(send_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD);
        MPI_Send(send_temp_grid_bottom, N, MPI_DOUBLE, processor_id + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_temp_grid_bottom, N, MPI_DOUBLE, processor_id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }    
    }
    else{
        for (int j = 0; j < N; j++) {
            recv_temp_grid_bottom[j] = BOTTOM_TEMP;
            recv_temp_grid_top[j] = TOP_TEMP;
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double local_max_diff = 0.0;
    double max_diff = 0.0;

    for (int i = 0; i < local_size; i++) {
        for (int j = 1; j < N-1; j++) {
            if ((i + j+initial_position) % 2 == 0) { 
                
                double old_val = grid_i[i][j];
                

                if(local_size==1){
                    grid_i[i][j] = 0.25 * (
                        grid_i[i][j+1] + 
                        grid_i[i][j-1] 
                        + recv_temp_grid_bottom[j] +
                        recv_temp_grid_top[j]); 
                }
                else{
                    if(i==0){
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            + grid_i[i+1][j] +
                            recv_temp_grid_top[j]);

                    }
                    else if(i ==local_size-1){
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            + recv_temp_grid_bottom[j]+
                            grid_i[i- 1][j]
                        );}
                    else{
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            +grid_i[i+1][j] 
                            +grid_i[i-1][j] 
                        );
                        
                    }
                }


                double diff = fabs(grid_i[i][j] - old_val);
                if (diff > local_max_diff) local_max_diff = diff;

            }
            
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (num_of_proc!=1){

        if (processor_id == 0) { 
            for (int j = 0; j < N; j++) {
                send_temp_grid_bottom[j] = grid_i[local_size - 1][j];
                recv_temp_grid_top[j] = TOP_TEMP;
            }
           MPI_Send(send_temp_grid_bottom, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
           MPI_Recv(recv_temp_grid_bottom, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } 
        else if (processor_id == num_of_proc - 1) {
            for (int j = 0; j < N; j++) {
                send_temp_grid_top[j] = grid_i[0][j]; 
                recv_temp_grid_bottom[j] = BOTTOM_TEMP;
    
            }
            MPI_Recv(recv_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(send_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD);
        } 
        else {
            for (int j = 0; j < N; j++) {
                send_temp_grid_top[j] = grid_i[0][j];             
                send_temp_grid_bottom[j] = grid_i[local_size - 1][j]; 
            }
            MPI_Recv(recv_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(send_temp_grid_top, N, MPI_DOUBLE, processor_id - 1, 0, MPI_COMM_WORLD);
            MPI_Send(send_temp_grid_bottom, N, MPI_DOUBLE, processor_id + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(recv_temp_grid_bottom, N, MPI_DOUBLE, processor_id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
        }    
        }
        else{
            for (int j = 0; j < N; j++) {
                recv_temp_grid_bottom[j] = BOTTOM_TEMP;
                recv_temp_grid_top[j] = TOP_TEMP;
            }
    
        }
    MPI_Barrier(MPI_COMM_WORLD);


    for (int i = 0; i < local_size; i++) {
        for (int j = 1; j < N-1; j++) {
            if ((i + j+initial_position) % 2 != 0) { 
                double old_val = grid_i[i][j];
                if(local_size==1){
                    grid_i[i][j] = 0.25 * (
                        grid_i[i][j+1] + 
                        grid_i[i][j-1] 
                        + recv_temp_grid_bottom[j] +
                        recv_temp_grid_top[j]); 
                }
                else{
                    if(i==0){
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            + grid_i[i+1][j] +
                            recv_temp_grid_top[j]);
                    }
                    else if(i ==local_size-1){
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            + recv_temp_grid_bottom[j]+
                            grid_i[i- 1][j]
                        );}
                    else{
                        grid_i[i][j] = 0.25 * (
                            grid_i[i][j+1] + 
                            grid_i[i][j-1] 
                            +grid_i[i+1][j] 
                            +grid_i[i-1][j] 
                        );
                        
                    }
                }
                
                double diff = fabs(grid_i[i][j] - old_val);

                if (diff > local_max_diff) local_max_diff = diff;

            }
        }

    }
    free(recv_temp_grid_top);
    free(recv_temp_grid_bottom);
    free(send_temp_grid_top);
    free(send_temp_grid_bottom);

    MPI_Allreduce(&local_max_diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max_diff;
}



int main(int argc, char* argv[]) {
    MPI_Init( &argc, &argv );
    MPI_Status status;
    int iter = 0;
    double max_diff;
    int processor_id;
    int num_of_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);

    int initial_position;
    int final_position;


    int rows_per_proc = (N - 2) / num_of_proc;  
    int remainder = (N - 2) % num_of_proc;      

    if (processor_id < remainder) {
        initial_position = 1 + processor_id * (rows_per_proc + 1);
        final_position = initial_position + (rows_per_proc + 1);
    } else {
        initial_position = 1 + remainder * (rows_per_proc + 1) + (processor_id - remainder) * rows_per_proc;
        final_position = initial_position + rows_per_proc;
    }

    int local_size = final_position - initial_position;
    double** grid_i = (double**)malloc(local_size * sizeof(double*));
    
    for (int i = 0; i < local_size; i++) {
        grid_i[i] = (double*)malloc(N * sizeof(double));
    }
    

    for (int i = initial_position; i < final_position; i++) {
        
        grid_i[i-initial_position][0] = LEFT_TEMP; 
        grid_i[i-initial_position][N-1] = RIGHT_TEMP; 
        
        for (int j = 1; j < N-1; j++) {
            grid_i[i-initial_position][j] = 0.0; 
        }
    }
    clock_t start_time = clock();
    do {
        max_diff = update_red_black(grid_i, local_size, initial_position);
        iter++;
    } while (max_diff > TOL && iter < MAX_ITER);

    clock_t end_time = clock();
    
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    if(processor_id==0){
    printf("Total computation time: %.6f seconds, processor_num = %d \n", elapsed_time, num_of_proc);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (processor_id == 0) {
        FILE* f = fopen("grid_parallel.txt", "w");

        {
            for (int j = 0; j < N; j++) {
                fprintf(f, "%.6f ", TOP_TEMP);
            }
            fprintf(f, "\n");
        }


        for (int i = 0; i < local_size; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(f, "%.6f ", grid_i[i][j]);
            }
            fprintf(f, "\n");
        }


        for (int src_rank = 1; src_rank < num_of_proc; src_rank++) {

            int local_size;
            MPI_Recv(&local_size, 1, MPI_INT, src_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            
            double* buf = (double*)malloc(local_size * N * sizeof(double));
            for (int i = 0; i < local_size; i++) {
                MPI_Recv(buf + i*N, N, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);
            }
            
            for (int row = 0; row < local_size; row++) {
                for (int col = 0; col < N; col++) {
                    fprintf(f, "%.6f ", buf[row * N + col]);

                }
                fprintf(f, "\n");
                
            }
            fflush(f);
            free(buf);
        }

        fprintf(f, "%.6f ", LEFT_TEMP); //forcing the bug with the noncontinuous corners
        for (int j = 1; j < N; j++) {
            fprintf(f, "%.6f ", BOTTOM_TEMP);
        }
        fprintf(f, "\n");
        

        fclose(f);
    } else {
        MPI_Send(&local_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        for (int i = 0; i < local_size; i++) {
            MPI_Send(grid_i[i], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    for (int i = 0; i < local_size; i++) {
        free(grid_i[i]);
    }
    free(grid_i);


    MPI_Finalize();
    return 0;
}