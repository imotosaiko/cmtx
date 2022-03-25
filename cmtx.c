#include <stdio.h>
#include <stdlib.h>

#include "cmtx.h"

void cw(mtx_t mtx_r, mtx_t mtx1, mtx_t mtx2, order_t n, mtx_t mtxs[]) {
    if (n == 1) {
        *mtx_r = *mtx1 * *mtx2;
        return;
    }
    
    mtx_t spaces = malloc(sizeof(element_t)*n/2*n/2*13);
    mtx_t a11 = &spaces[0];
    mtx_t a12 = &spaces[n/2*n/2];
    mtx_t a21 = mtxs[0];
    mtx_t a22 = &spaces[n/2*n/2*2];
    mtx_t b11 = &spaces[n/2*n/2*3];
    mtx_t b12 = mtxs[1];
    mtx_t b21 = &spaces[n/2*n/2*4];
    mtx_t b22 = &spaces[n/2*n/2*5];
    
    part_mtx(mtx1, a11, a12, a21, a22, n);
    part_mtx(mtx2, b11, b12, b21, b22, n);

    mtx_t s1 = &mtx1[0];
    mtx_t s2 = &mtx1[n/2*n/2];
    mtx_t s3 = &mtx1[n/2*n/2*2];
    mtx_t s4 = &mtx1[n/2*n/2*3];
    mtx_t t1 = &mtx2[0];
    mtx_t t2 = &mtx2[n/2*n/2];
    mtx_t t3 = &mtx2[n/2*n/2*2];
    mtx_t t4 = &mtx2[n/2*n/2*3];
    
    mtx_add(s1, a21, a22, n/2, n/2);
    mtx_sub(s2, s1, a11, n/2, n/2);
    mtx_sub(s3, a11, a21, n/2, n/2);
    mtx_sub(s4, a12, s2, n/2, n/2);
    mtx_sub(t1, b12, b11, n/2, n/2);
    mtx_sub(t2, b22, t1, n/2, n/2);
    mtx_sub(t3, b22, b12, n/2, n/2);
    mtx_sub(t4, t2, b21, n/2, n/2);

    mtx_t m1 = &spaces[n/2*n/2*6];
    mtx_t m2 = &spaces[n/2*n/2*7];
    mtx_t m3 = &spaces[n/2*n/2*8];
    mtx_t m4 = &spaces[n/2*n/2*9];
    mtx_t m5 = &spaces[n/2*n/2*10];
    mtx_t m6 = &spaces[n/2*n/2*11];
    mtx_t m7 = &spaces[n/2*n/2*12];

    cw(m1, a11, b11, n/2, mtxs);
    cw(m2, a12, b21, n/2, mtxs);
    cw(m3, s4, b22, n/2, mtxs);
    cw(m4, a22, t4, n/2, mtxs);
    cw(m5, s1, t1, n/2, mtxs);
    cw(m6, s2, t2, n/2, mtxs);
    cw(m7, s3, t3, n/2, mtxs);    

    mtx_t u1 = mtxs[2];
    mtx_t u2 = mtxs[3];
    mtx_t u3 = mtxs[4];
    mtx_t u4 = mtxs[5];
    mtx_t u5 = mtxs[6];
    mtx_t u6 = mtxs[7];
    mtx_t u7 = mtxs[8];
   
    mtx_add(u1, m1, m2, n/2, n/2);
    mtx_add(u2, m1, m6, n/2, n/2);
    mtx_add(u3, u2, m7, n/2, n/2);
    mtx_add(u4, u2, m5, n/2, n/2);
    mtx_add(u5, u4, m3, n/2, n/2);
    mtx_sub(u6, u3, m4, n/2, n/2);
    mtx_add(u7, u3, m5, n/2, n/2);
    
    mtx_t c11 = u1;
    mtx_t c12 = u5;
    mtx_t c21 = u6;
    mtx_t c22 = u7;

    comb_mtx(mtx_r, c11, c12, c21, c22, n);

    free(spaces);
}
