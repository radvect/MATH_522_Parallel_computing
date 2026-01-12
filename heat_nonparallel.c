#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 512// Grid size (50x50)
#define TOL 1e-6
#define MAX_ITER 100000

// Boundary Conditions
#define LEFT_TEMP 1.2
#define RIGHT_TEMP 1.3
#define TOP_TEMP 1.5
#define BOTTOM_TEMP 1.1

void initialize_grid(double grid[N][N]) {
    // Set boundary conditions
    for (int i = 0; i < N; i++) {
        grid[i][0] = LEFT_TEMP;
        grid[i][N-1] = RIGHT_TEMP;
        grid[0][i] = TOP_TEMP;
        grid[N-1][i] = BOTTOM_TEMP;
    }
    // Interior points initial guess
    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            grid[i][j] = 0.0;
        }
    }
}

double update_red_black(double grid[N][N]) {
    double max_diff = 0.0;

    // Red Update (even sum indices)
    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            if ((i + j) % 2 == 0) { 
                double old_val = grid[i][j];
                grid[i][j] = 0.25 * (grid[i+1][j] + grid[i-1][j] +
                                     grid[i][j+1] + grid[i][j-1]);
                double diff = fabs(grid[i][j] - old_val);
                if (diff > max_diff) max_diff = diff;
            }
        }
    }

    // Black Update (odd sum indices)
    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            if ((i + j) % 2 != 0) { 
                double old_val = grid[i][j];
                grid[i][j] = 0.25 * (grid[i+1][j] + grid[i-1][j] +
                                     grid[i][j+1] + grid[i][j-1]);
                double diff = fabs(grid[i][j] - old_val);
                if (diff > max_diff) max_diff = diff;
            }
        }
    }

    return max_diff;
}

void print_grid(double grid[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%6.2f ", grid[i][j]);
        }
        printf("\n");
    }
}

int main() {
    double grid[N][N];
    initialize_grid(grid);

    int iter = 0;
    double max_diff;
    clock_t start_time = clock();
    do {
        max_diff = update_red_black(grid);
        //printf("max diff: %f\n", max_diff);
        iter++;
    } while (max_diff > TOL && iter < MAX_ITER);

    //printf("Converged in %d iterations with max diff = %e\n", iter, max_diff);
    //print_grid(grid);
    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Total computation time: %.6f seconds\n", elapsed_time);


    
    FILE *file = fopen("grid_nonparallel.txt", "w");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%.6f ", grid[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);


    return 0;
}
