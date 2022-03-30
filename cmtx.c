#include <stdio.h>
#include <stdlib.h>

#include "cmtx.h"

void cw(matrix_t matrix_r, matrix_t matrix1, matrix_t matrix2, order_t n, matrix_t matrixs) {
    if (n == 1) {
        *matrix_r = *matrix1 * *matrix2;
        return;
    }
    
    matrix_t spaces = malloc(sizeof(element_t)*n/2*n/2*13);
    matrix_t a11 = &spaces[0];
    matrix_t a12 = &spaces[n/2*n/2];
    matrix_t a21 = matrixs;
    matrix_t a22 = &spaces[n/2*n/2*2];
    matrix_t b11 = &spaces[n/2*n/2*3];
    matrix_t b12 = matrixs + n/2 * n/2;
    matrix_t b21 = &spaces[n/2*n/2*4];
    matrix_t b22 = &spaces[n/2*n/2*5];
    
    part_matrix(matrix1, a11, a12, a21, a22, n);
    part_matrix(matrix2, b11, b12, b21, b22, n);

    matrix_t s1 = &matrix1[0];
    matrix_t s2 = &matrix1[n/2*n/2];
    matrix_t s3 = &matrix1[n/2*n/2*2];
    matrix_t s4 = &matrix1[n/2*n/2*3];
    matrix_t t1 = &matrix2[0];
    matrix_t t2 = &matrix2[n/2*n/2];
    matrix_t t3 = &matrix2[n/2*n/2*2];
    matrix_t t4 = &matrix2[n/2*n/2*3];
    
    matrix_add(s1, a21, a22, n/2, n/2);
    matrix_sub(s2, s1, a11, n/2, n/2);
    matrix_sub(s3, a11, a21, n/2, n/2);
    matrix_sub(s4, a12, s2, n/2, n/2);
    matrix_sub(t1, b12, b11, n/2, n/2);
    matrix_sub(t2, b22, t1, n/2, n/2);
    matrix_sub(t3, b22, b12, n/2, n/2);
    matrix_sub(t4, t2, b21, n/2, n/2);

    matrix_t m1 = &spaces[n/2*n/2*6];
    matrix_t m2 = &spaces[n/2*n/2*7];
    matrix_t m3 = &spaces[n/2*n/2*8];
    matrix_t m4 = &spaces[n/2*n/2*9];
    matrix_t m5 = &spaces[n/2*n/2*10];
    matrix_t m6 = &spaces[n/2*n/2*11];
    matrix_t m7 = &spaces[n/2*n/2*12];

    cw(m1, a11, b11, n/2, matrixs);
    cw(m2, a12, b21, n/2, matrixs);
    cw(m3, s4, b22, n/2, matrixs);
    cw(m4, a22, t4, n/2, matrixs);
    cw(m5, s1, t1, n/2, matrixs);
    cw(m6, s2, t2, n/2, matrixs);
    cw(m7, s3, t3, n/2, matrixs);    

    matrix_t u1 = matrixs + n/2 * n/2 * 2;
    matrix_t u2 = matrixs + n/2 * n/2 * 3;
    matrix_t u3 = matrixs + n/2 * n/2 * 4;
    matrix_t u4 = matrixs + n/2 * n/2 * 5;
    matrix_t u5 = matrixs + n/2 * n/2 * 6;
    matrix_t u6 = matrixs + n/2 * n/2 * 7;
    matrix_t u7 = matrixs + n/2 * n/2 * 8;
   
    matrix_add(u1, m1, m2, n/2, n/2);
    matrix_add(u2, m1, m6, n/2, n/2);
    matrix_add(u3, u2, m7, n/2, n/2);
    matrix_add(u4, u2, m5, n/2, n/2);
    matrix_add(u5, u4, m3, n/2, n/2);
    matrix_sub(u6, u3, m4, n/2, n/2);
    matrix_add(u7, u3, m5, n/2, n/2);
    
    matrix_t c11 = u1;
    matrix_t c12 = u5;
    matrix_t c21 = u6;
    matrix_t c22 = u7;

    comb_matrix(matrix_r, c11, c12, c21, c22, n);

    free(spaces);
}
