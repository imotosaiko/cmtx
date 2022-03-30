#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cmtx.h"

// #define DEBUG
#define REPEAT_TIMES 10

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("error!\n"
               "please pass the order(512 to 2048) of matrixs.\n"
               "usage: ./test num1 num2 num3\n");
        return -1;
    }

    // init two random matrixs and a matrix for result
    order_t m = atoi(argv[1]);
    order_t n = atoi(argv[2]);
    order_t k = atoi(argv[3]);
    matrix_t matrixs = new_matrix((m*n + n*k + m*k));
    matrix_t matrix1 = matrixs;
    matrix_t matrix2 = matrixs + m*n;
    matrix_t matrix_r = matrix2 + n*k;
    rand_matrix(matrix1, m * n);
    rand_matrix(matrix2, n * k);

    // Calculation and Timing
    double average_time;
    for (int i = 0; i < REPEAT_TIMES; i++) {
        clock_t start = clock();

        // matrix_mul_general(matrix_r, matrix1, matrix2, m, n, k);
        matrix_mul_cw(matrix_r, matrix1, matrix2, m, n, k);
        
        clock_t end = clock();
        double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        average_time += total_time / REPEAT_TIMES;
    }    

    #ifdef DEBUG    
    disp_matrix("matrix1", matrix1, m, n);
    disp_matrix("matrix2", matrix2, n, k);
    disp_matrix("matrix_r", matrix_r, m, k);
    #endif
    
    // output statistics
    printf("average time used: %fs\n", average_time);

    // free space
    free(matrixs);
    
    return 0;
}
