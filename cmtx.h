#ifndef CMTX_H
#define CMTX_H

typedef unsigned int order_t;
typedef int element_t;

#define ELEMENT_VALUE_MAX 255

// display a matrix
#define disp_mtx(name, mtx, m, n) {                 \
        printf("the matrix %s: \n", name);          \
        for (order_t i = 0; i < m; i++) {           \
            for (order_t j = 0; j < n; j++) {       \
                printf("%d ", *(mtx + i*n+j));      \
            }                                       \
            printf("\n");                           \
        }                                           \
    }

// get a random matrix
#define rand_mtx(mtx, m, n) {                           \
        for (order_t i = 0; i < m*n; i++) {             \
            *(mtx+i) = rand() % ELEMENT_VALUE_MAX;      \
        }                                               \
    }

// general matrix multiplication
#define mtx_mul_general(mtx_r, mtx1, mtx2, m, n, k) {   \
        for (order_t i = 0; i < m; i++) {               \
            for (order_t j = 0; j < n; j++) {           \
                for(order_t f = 0; f < k; f++) {                   \
                    *(mtx_r + i*k+f) += *(mtx1 + i*n+j) * *(mtx2 + j*k+f); \
                }                                                       \
            }                                                           \
        }                                                               \
    }

#endif
