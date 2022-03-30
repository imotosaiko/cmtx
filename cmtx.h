#ifndef CMTX_H
#define CMTX_H

typedef unsigned int order_t;
typedef int element_t;
typedef element_t* matrix_t;

#define ELEMENT_VALUE_MAX 255

// create a new matrix with counts elements
#define new_matrix(counts) malloc(sizeof(element_t)*counts)

// display a matrix
#define disp_matrix(name, matrix, m, n) {       \
        printf("the matrix %s: \n", name);      \
        for (order_t i = 0; i < m; i++) {       \
            for (order_t j = 0; j < n; j++) {   \
                printf("%d ", *(matrix + i*n+j));  \
            }                                   \
            printf("\n");                       \
        }                                       \
    }

// get a random matrix
#define rand_matrix(matrix, counts) {                      \
        for (order_t i = 0; i < counts; i++) {              \
            *(matrix + i) = rand() % ELEMENT_VALUE_MAX;  \
        }                                                \
    }

// general matrix multiplication
#define matrix_mul_general(matrix_r,              \
                           matrix1,               \
                           matrix2,               \
                           m, n, k) {             \
        for (order_t i = 0; i < m; i++) {         \
            for (order_t j = 0; j < n; j++) {     \
                for(order_t f = 0; f < k; f++) {  \
                    *(matrix_r + i*k+f) +=        \
                        *(matrix1 + i*n+j) *      \
                        *(matrix2 + j*k+f);       \
                }                                 \
            }                                     \
        }                                         \
    }

// Find the nearest integer power of 2 to num
#define nearby_2_power(num,             \
                       expanded_num) {  \
        expanded_num = num;             \
        expanded_num -= 1;              \
        expanded_num |= expanded_num >> 16;      \
        expanded_num |= expanded_num >> 8;       \
        expanded_num |= expanded_num >> 4;       \
        expanded_num |= expanded_num >> 2;       \
        expanded_num |= expanded_num >> 1;       \
        expanded_num += 1;              \
    }

// Expand the matrix into a square matrix with the extra elements complemented by 0
#define expand_matrix(expanded_matrix, matrix,                          \
                      m, n, expanded_order) {                           \
        for (order_t i = 0; i < m; i++) {                               \
            for(order_t j = 0; j < n; j++) {                            \
                *(expanded_matrix + i*expanded_order+j) = *(matrix + i*n+j); \
            }                                                           \
        }                                                               \
    }

// Accept any two matrices, multiply them with the Coppersmith–Winograd algorithm
#define matrix_mul_cw(matrix_r,                                         \
                      matrix1,                                          \
                      matrix2,                                          \
                      m, n, k) {                                        \
        int max_order = m > n ? m : n;                                  \
        max_order = max_order > k ? max_order : k;                      \
        int expanded_order;                                             \
        nearby_2_power(max_order, expanded_order);                      \
        matrix_t expanded_matrixs = malloc(sizeof(element_t) *          \
                                           expanded_order * expanded_order * 3); \
        matrix_t expanded_matrix1 = expanded_matrixs;                   \
        matrix_t expanded_matrix2 = expanded_matrix1 + expanded_order * expanded_order; \
        matrix_t expanded_matrix_r = expanded_matrix2 + expanded_order * expanded_order; \
        expand_matrix(expanded_matrix1, matrix1, m, n, expanded_order); \
        expand_matrix(expanded_matrix2, matrix2, n, k, expanded_order); \
        size_t matrixs_size = (32 * expanded_order * expanded_order / 4 - 32) / 3; \
        matrix_t matrixs = malloc(sizeof(element_t) * matrixs_size); \
        cw(expanded_matrix_r,                                           \
           expanded_matrix1, expanded_matrix2,                          \
           expanded_order, matrixs);                                  \
        for (order_t i = 0; i < m; i++) {                               \
            for (order_t j = 0; j < k; j++) {                           \
                *(matrix_r + i*k+j) = *(expanded_matrix_r + i*expanded_order+j); \
            }                                                           \
        }                                                               \
        free(expanded_matrixs);                                         \
        free(matrixs);                                                  \
    }

// Dividing a matrix equally into four pieces
#define part_matrix(matrix,       \
                    sub_matrix11, \
                    sub_matrix12, \
                    sub_matrix21, \
                    sub_matrix22, n) {                                                \
        for (order_t i = 0; i < n/2; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(sub_matrix11 + i*n/2+j) = *(matrix + i*n+j);          \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(sub_matrix12 + i*n/2+(j-n/2)) = *(matrix + i*n+j);    \
            }                                                           \
        }                                                               \
        for (order_t i = n/2; i < n; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(sub_matrix21 + (i-n/2)*n/2+j) = *(matrix + i*n+j);    \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(sub_matrix22 + (i-n/2)*n/2+(j-n/2)) = *(matrix + i*n+j); \
            }                                                           \
        }                                                               \
    }

// Combine four matrices of the same size into one
#define comb_matrix(matrix,                                             \
                    sub_matrix11,                                       \
                    sub_matrix12,                                       \
                    sub_matrix21,                                       \
                    sub_matrix22, n) {                                  \
        for (order_t i = 0; i < n/2; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(matrix + i*n+j) = *(sub_matrix11 + i*n/2+j);          \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(matrix + i*n+j) = *(sub_matrix12 + i*n/2+(j-n/2));    \
            }                                                           \
        }                                                               \
        for (order_t i = n/2; i < n; i++) {                             \
            for (order_t j = 0; j < n/2; j++) {                         \
                *(matrix + i*n+j) = *(sub_matrix21 + (i-n/2)*n/2+j);    \
            }                                                           \
            for (order_t j = n/2; j < n; j++) {                         \
                *(matrix + i*n+j) = *(sub_matrix22 + (i-n/2)*n/2+(j-n/2)); \
            }                                                           \
        }                                                               \
    }

#define matrix_add(matrix_r, matrix1, matrix2, m, n) {  \
        matrix_t matrix_r_p = matrix_r;                    \
        matrix_t matrix_r_end = matrix_r + m*n;            \
        matrix_t matrix1_p = matrix1;                      \
        matrix_t matrix2_p = matrix2;                         \
        while (matrix_r_p < matrix_r_end) {          \
            *matrix_r_p = *matrix1_p + *matrix2_p;   \
            matrix_r_p++;                            \
            matrix1_p++;                             \
            matrix2_p++;                             \
        }                                            \
    }

#define matrix_sub(matrix_r, matrix1, matrix2, m, n) {  \
        matrix_t matrix_r_p = matrix_r;                 \
        matrix_t matrix_r_end = matrix_r + m*n;         \
        matrix_t matrix1_p = matrix1;                   \
        matrix_t matrix2_p = matrix2;                   \
        while (matrix_r_p < matrix_r_end) {             \
            *matrix_r_p = *matrix1_p - *matrix2_p;      \
            matrix_r_p++;                               \
            matrix1_p++;                                \
            matrix2_p++;                                \
        }                                               \
    }

// Coppersmith–Winograd algorihm, Accepts two square matrices of order 2 integer powers
void cw(matrix_t matrix_r, \
        matrix_t matrix1, \
        matrix_t matrix2, \
        order_t n, matrix_t matrixs);

#endif
