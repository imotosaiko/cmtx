#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cmtx.h"

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
    element_t* mtx1 = malloc(sizeof(element_t)*m*n);
    element_t* mtx2 = malloc(sizeof(element_t)*n*k);
    element_t* mtx_r = malloc(sizeof(element_t)*m*k);
    rand_mtx(mtx1, m, n);
    rand_mtx(mtx2, n, k);

    // Calculation and Timing
    double average_time;
    for (int i = 0; i < REPEAT_TIMES; i++) {
        clock_t start = clock();
        mtx_mul_general(mtx_r, mtx1, mtx2, m, n, k);
        clock_t end = clock();
        double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        average_time += total_time / REPEAT_TIMES;
    }
    
    // output statistics
    printf("average time used: %fs\n", average_time);

    // free space
    free(mtx1);
    free(mtx2);
    free(mtx_r);

    return 0;
}
