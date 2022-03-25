#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cmtx.h"

// #define DEBUG
#define REPEAT_TIMES 1

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
    element_t* mtx1 = new_mtx(m, n);
    element_t* mtx2 = new_mtx(n, k);
    element_t* mtx_r = new_mtx(m, k);
    rand_mtx(mtx1, m, n);
    rand_mtx(mtx2, n, k);

    // Calculation and Timing
    double average_time;
    for (int i = 0; i < REPEAT_TIMES; i++) {
        clock_t start = clock();

        // mtx_mul_general(mtx_r, mtx1, mtx2, m, n, k);
        mtx_mul_cw(mtx_r, mtx1, mtx2, m, n, k);
        
        clock_t end = clock();
        double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        average_time += total_time / REPEAT_TIMES;
    }    

    #ifdef DEBUG
    disp_mtx("mtx1", mtx1, m, n);
    disp_mtx("mtx2", mtx2, n, k);
    disp_mtx("mtx_r", mtx_r, m, k);
    #endif
    
    // output statistics
    printf("average time used: %fs\n", average_time);

    // free space
    free(mtx1);
    free(mtx2);
    free(mtx_r);

    return 0;
}
