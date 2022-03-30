#include <stdio.h>
#include <stdlib.h>

#include "cmtx.h"

void cw(matrix_t matrix_r, matrix_t matrix1, matrix_t matrix2, order_t n, matrix_t matrixs) {
    if (n == 1) {
        *matrix_r = *matrix1 * *matrix2;
        return;
    }

    matrix_t a11 = matrixs;
    matrix_t a12 = a11 + n/2 * n/2;
    matrix_t a21 = a12 + n/2 * n/2;
    matrix_t a22 = a21 + n/2 * n/2;
    part_matrix(matrix1, a11, a12, a21, a22, n);
    matrix_t b11 = matrix1;
    matrix_t b12 = b11 + n/2 * n/2;
    matrix_t b21 = b12 + n/2 * n/2;
    matrix_t b22 = b21 + n/2 * n/2;
    part_matrix(matrix2, b11, b12, b21, b22, n);

    matrix_t s1 = matrix2;
    matrix_t s2 = s1 + n/2 * n/2;
    matrix_t s3 = s2 + n/2 * n/2;
    matrix_t s4 = s3 + n/2 * n/2;
    matrix_t t1 = a22 + n/2 * n/2;
    matrix_sub(t1, b12, b11, n/2, n/2);
    matrix_t t2 = t1 + n/2 * n/2;
    matrix_t t3 = t2 + n/2 * n/2;
    matrixs = t3 + n/2 * n/2;
    matrix_sub(t3, b22, b12, n/2, n/2);
    matrix_t t4 = b12;
    
    matrix_add(s1, a21, a22, n/2, n/2);
    matrix_sub(s2, s1, a11, n/2, n/2);
    matrix_sub(s3, a11, a21, n/2, n/2);
    matrix_sub(s4, a12, s2, n/2, n/2);

    matrix_sub(t2, b22, t1, n/2, n/2);

    matrix_sub(t4, t2, b21, n/2, n/2);

    matrix_t m1 = a21;
    cw(m1, a11, b11, n/2, matrixs);
    matrix_t m2 = a11;
    cw(m2, a12, b21, n/2, matrixs);
    matrix_t m3 = a12;
    cw(m3, s4, b22, n/2, matrixs);
    matrix_t m4 = s4;
    cw(m4, a22, t4, n/2, matrixs);
    matrix_t m5 = a22;
    cw(m5, s1, t1, n/2, matrixs);
    matrix_t m6 = s1;
    cw(m6, s2, t2, n/2, matrixs);
    matrix_t m7 = s2;
    cw(m7, s3, t3, n/2, matrixs);

    matrix_t u1 = b11;
    matrix_t u2 = b12;
    matrix_t u3 = b21;
    matrix_t u4 = b22;
    matrix_t u5 = s3;
    matrix_t u6 = s4;
    matrix_t u7 = t1;
   
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
}
