#include <stdio.h>
#include <stdlib.h>

#include "cmtx.h"

void cw(mtx_t mtx_r, mtx_t mtx1, mtx_t mtx2, order_t n, mtx_t mtxs[]) {
    if (n == 1) {
        *mtx_r = *mtx1 * *mtx2;
        return;
    }

    mtx_t a11 = new_mtx(n/2, n/2);
    mtx_t a12 = new_mtx(n/2, n/2);
    mtx_t a21 = mtxs[0];
    mtx_t a22 = new_mtx(n/2, n/2);
    mtx_t b11 = new_mtx(n/2, n/2);
    mtx_t b12 = mtxs[1];
    mtx_t b21 = new_mtx(n/2, n/2);
    mtx_t b22 = new_mtx(n/2, n/2);
    
    part_mtx(mtx1, a11, a12, a21, a22, n);
    part_mtx(mtx2, b11, b12, b21, b22, n);

    mtx_t s1 = new_mtx(n/2, n/2);
    mtx_t s2 = new_mtx(n/2, n/2);
    mtx_t s3 = new_mtx(n/2, n/2);
    mtx_t s4 = new_mtx(n/2, n/2);
    mtx_t t1 = new_mtx(n/2, n/2);
    mtx_t t2 = new_mtx(n/2, n/2);
    mtx_t t3 = new_mtx(n/2, n/2);
    mtx_t t4 = new_mtx(n/2, n/2);
    
    mtx_add(s1, a21, a22, n/2, n/2);
    mtx_sub(s2, s1, a11, n/2, n/2);
    mtx_sub(s3, a11, a21, n/2, n/2);
    mtx_sub(s4, a12, s2, n/2, n/2);
    mtx_sub(t1, b12, b11, n/2, n/2);
    mtx_sub(t2, b22, t1, n/2, n/2);
    mtx_sub(t3, b22, b12, n/2, n/2);
    mtx_sub(t4, t2, b21, n/2, n/2);
    
    mtx_t m1 = new_mtx(n/2, n/2);
    mtx_t m2 = new_mtx(n/2, n/2);
    mtx_t m3 = new_mtx(n/2, n/2);
    mtx_t m4 = new_mtx(n/2, n/2);
    mtx_t m5 = new_mtx(n/2, n/2);
    mtx_t m6 = new_mtx(n/2, n/2);
    mtx_t m7 = new_mtx(n/2, n/2);

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

    free(a11);
    free(a12);
    free(a22);
    free(b11);
    free(b21);
    free(b22);

    free(s1);
    free(s2);
    free(s3);
    free(s4);
    free(t1);
    free(t2);
    free(t3);
    free(t4);

    free(m1);
    free(m2);
    free(m3);
    free(m4);
    free(m5);
    free(m6);
    free(m7);
}
