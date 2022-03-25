#ifndef CMTX_H
#define CMTX_H

typedef unsigned int order_t;
typedef int element_t;
typedef element_t* mtx_t;

#define ELEMENT_VALUE_MAX 255

// create a new matrix
#define new_mtx(m, n) malloc(sizeof(element_t)*m*n)
   

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
#define mtx_mul_general(mtx_r, mtx1, mtx2, m, n, k) {                   \
        for (order_t i = 0; i < m; i++) {                               \
            for (order_t j = 0; j < n; j++) {                           \
                for(order_t f = 0; f < k; f++) {                        \
                    *(mtx_r + i*k+f) += *(mtx1 + i*n+j) * *(mtx2 + j*k+f); \
                }                                                       \
            }                                                           \
        }                                                               \
    }

// Find the nearest integer power of 2 to num
#define nearby_2_power(num, expanded_num) {    \
        order_t tmp = num;                     \
        tmp -= 1;                              \
        tmp |= tmp >> 16;                      \
        tmp |= tmp >> 8;                       \
        tmp |= tmp >> 4;                       \
        tmp |= tmp >> 2;                       \
        tmp |= tmp >> 1;                       \
        tmp += 1;                              \
        expanded_num = tmp;                    \
    }

// Expand the matrix into a square matrix with the extra elements complemented by 0
#define expand_mtx(mtx_r, mtx, m, n, expanded_order) {                    \
        for (order_t i = 0; i < m; i++) {                               \
            for(order_t j = 0; j < n; j++) {                            \
                *(mtx_r + i*expanded_order+j) = *(mtx + i*n+j);                       \
            }                                                           \
        }                                                               \
    }

// Accept any two matrices, multiply them with the Coppersmith–Winograd algorithm
#define mtx_mul_cw(mtx_r, mtx1, mtx2, m, n, k) {                        \
        int max_order = m > n ? m : n;                                  \
        max_order = max_order > k ? max_order : k;                      \
        int expanded_order;                                             \
        nearby_2_power(max_order, expanded_order);                      \
        mtx_t expanded_mtx1 = new_mtx(expanded_order, expanded_order);  \
        mtx_t expanded_mtx2 = new_mtx(expanded_order, expanded_order);  \
        mtx_t expanded_mtx_r = new_mtx(expanded_order, expanded_order); \
        expand_mtx(expanded_mtx1, mtx1, m, n, expanded_order);          \
        expand_mtx(expanded_mtx2, mtx2, n, k, expanded_order);          \
        int mtxs_size = 21;                                       \
        mtx_t mtxs[mtxs_size];                              \
        for (int i = 0; i < mtxs_size; i++) {                    \
            mtxs[i] = new_mtx(expanded_order/2, expanded_order/2); \
        }                                                               \
        cw(expanded_mtx_r, expanded_mtx1, expanded_mtx2, expanded_order, mtxs); \
        for (order_t i = 0; i < m; i++) {                               \
            for (order_t j = 0; j < k; j++) {                           \
                *(mtx_r + i*k+j) = *(expanded_mtx_r + i*expanded_order+j); \
            }                                                           \
        }                                                               \
        for (int i = 0; i < mtxs_size; i++) {                     \
            free(mtxs[i]);                                         \
        }                                                               \
        free(expanded_mtx1);                                            \
        free(expanded_mtx2);                                            \
        free(expanded_mtx_r);                                           \
    }

// Dividing a matrix equally into four pieces
#define part_mtx(mtx, sub_mtx11, sub_mtx12, sub_mtx21, sub_mtx22, n) {  \
        for (order_t i = 0; i < n/2; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(sub_mtx11 + i*n/2+j) = *(mtx + i*n+j);                \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(sub_mtx12 + i*n/2+(j-n/2)) = *(mtx + i*n+j);          \
            }                                                           \
        }                                                               \
        for (order_t i = n/2; i < n; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(sub_mtx21 + (i-n/2)*n/2+j) = *(mtx + i*n+j);         \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(sub_mtx22 + (i-n/2)*n/2+(j-n/2)) = *(mtx + i*n+j);    \
            }                                                           \
        }                                                               \
    }

// Combine four matrices of the same size into one
#define comb_mtx(mtx, sub_mtx11, sub_mtx12, sub_mtx21, sub_mtx22, n) {  \
        for (order_t i = 0; i < n/2; i++) {                               \
            for (order_t j = 0; j < n/2; j++) {                           \
                *(mtx + i*n+j) = *(sub_mtx11 + i*n/2+j);              \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(mtx + i*n+j) = *(sub_mtx12 + i*n/2+(j-n/2));          \
            }                                                           \
        }                                                               \
        for (order_t i = n/2; i < n; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(mtx + i*n+j) = *(sub_mtx21 + (i-n/2)*n/2+j);          \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(mtx + i*n+j) = *(sub_mtx22 + (i-n/2)*n/2+(j-n/2));    \
            }                                                           \
        }                                                               \
    }

#define mtx_add(mtx_r, mtx1, mtx2, m, n) {                              \
        for (order_t i = 0; i < m; i++) {                               \
            for (order_t j = 0; j < n; j++) {                           \
                *(mtx_r + i*n+j) = *(mtx1 + i*n+j) + *(mtx2 + i*n+j);   \
            }                                                           \
        }                                                               \
    }

#define mtx_sub(mtx_r, mtx1, mtx2, m, n) {      \
        for (order_t i = 0; i < m; i++) {                               \
            for (order_t j = 0; j < n; j++) {                           \
                *(mtx_r + i*n+j) = *(mtx1 + i*n+j) - *(mtx2 + i*n+j);   \
            }                                                           \
        }                                                               \
    }

#define mtx_cp(mtx_dist, mtx_src, m, n) {       \
        for (order_t i = 0; i < m; i++) {       \
            for (order_t j = 0; j < n; j++) {   \
                *(mtx_dist + i*n+j) = *(mtx_src + i*n+j);       \
            }                                   \
        }                                       \
    }

// Coppersmith–Winograd algorihm, Accepts two square matrices of order 2 integer powers
void cw(mtx_t mtx_r, mtx_t mtx1, mtx_t mtx2, order_t n, mtx_t mtxs[]);

#endif
