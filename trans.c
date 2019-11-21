/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

void print_2Darray(int m, int (*mat)[m]) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            printf("%d\t", mat[i][j]);
        }
        printf("\n");
    }
}
    
int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */

int min(int a, int b) {
    return (a < b) ? a : b;
}

char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    /* Finding the transpose of a matrix with as few cache misses as possible. The algorithm below uses the technique of
     * blocking to improve the spatial locality of the accesses. The idea is to minimize conflict misses by ensuring that 
     * all accesses to the cache block are completed while it is in the cache rather than evicting it and bringing it
     * back again and again. This helps reduce conflict misses but compulsory misses are still needed to bring the data 
     * into the cache. */

    /* The cache under consideration is a 1KB direct mapped cache with 32 byte blocks, hence there are 32 sets. Since 8 ints
     * fit in one block, the elements in each row map to M/8 distinct sets. Since the 2D arrays are laid out in row-major order,
     * this means that the access pattern for the stores in the matrix transpose would involve going down a column in
     * matrix B and the elements would map to set indices separated by M/8. In case of the 32*32 matrix, M/8 = 4 and we
     * may have a set access pattern like 8, 12, 16 ... But since we only have 32 sets, going beyond 8 rows would evict a 
     * block which is still needed resulting in conflict miss when the evicted block is needed again. The same holds for the 64*64 
     * matrix except that in this case, 4 rows is the threshold rather than 8 */

    /* The above argument only holds for matrices where number of columns is a multiple of 8 (32*32 and 64*64) since
     * blocks of 8 ints will not neatly conform to the rows of the matrix. The best block size is found to be 18 by trial 
     * and error for the 61*67 matrix */

    int num_sets = 32;
    int sets_per_row = M/8;
    int block_dim = (M == N) ? num_sets/sets_per_row : 18;
    int idx1, idx2, diag_flag = 0;
    for (int i = 0; i < N; i += block_dim) {
        for (int j = 0; j < M; j += block_dim) {
            
            /* Transposing the smaller block matrices */
            for (int ib = i; ib < min(i+block_dim, N); ib++) {
                for (int jb = j; jb < min(j+block_dim, M); jb++) {
                    /* Conflict misses along the diagonal are particularly tricky since the the element in A being
                     * "load"ed and the index in B where it is being "store"d map to the same set and one operation evicts the
                     * other data which will be accessed in the near future. So, the diagonal elements are dealt with separately 
                     * at the end after all the other accesses are processed and the the blocks are no longer needed */
                    if ((ib != jb) && ((M != 64) || ((ib != jb-4) && (ib != jb+4)))) {
                        B[jb][ib] = A[ib][jb];
                    } else {
                        diag_flag = 1;
                        idx1 = ib;
                        idx2 = jb;
                    }
                }
                if (diag_flag) {
                    B[idx2][idx1] = A[idx1][idx2];
                    diag_flag = 0;
                }
            }
        }
    }
}

  
/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 
/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            B[j][i] = A[i][j];
        }
    }    

}

char transpose_blocking_desc[] = "Simple blocking";
void transpose_blocking(int M, int N, int A[N][M], int B[M][N])
{
    int num_sets = 32;
    int sets_per_row = M/8;
    int block_dim = num_sets/sets_per_row;
    for (int i = 0; i < N; i += block_dim) {
        for (int j = 0; j < M; j += block_dim) {

            for (int ib = i; ib < min(i+block_dim, N); ib++) {
                for (int jb = j; jb < min(j+block_dim, M); jb++) {
                    B[jb][ib] = A[ib][jb];
                }
            }
        }
    }
}

char transpose_square_matrix_desc[] = "Transpose blocking with block size 8 for 32*32 and 4 for 64*64";
void transpose_square_matrix(int M, int N, int A[N][M], int B[M][N])
{
    int num_sets = 32;
    int sets_per_row = M/8;
    int block_dim = num_sets/sets_per_row;
    int idx1, idx2, diag_flag = 0;
    for (int i = 0; i < N; i += block_dim) {
        for (int j = 0; j < M; j += block_dim) {

            for (int ib = i; ib < min(i+block_dim, N); ib++) {
                for (int jb = j; jb < min(j+block_dim, M); jb++) {
                    if ((ib != jb) && ((M == 32) || ((ib != jb-4) && (ib != jb+4)))) {
                        B[jb][ib] = A[ib][jb];
                    } else {
                        diag_flag = 1;
                        idx1 = ib;
                        idx2 = jb;
                    }
                }
                if (diag_flag) {
                    B[idx2][idx1] = A[idx1][idx2];
                    diag_flag = 0;
                }
            }
        }
    }
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 
    registerTransFunction(transpose_blocking, transpose_blocking_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

