#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gauss_fixed.h"
#if VERBOSE
#include "io.h"
#endif

#define COLS 2 * SIZE
#define AT(mat, i, j) (mat[(i) * (COLS) + (j)])

void swap(float_t exmat[TOTALSIZE], idx_t i1, idx_t i2)
{
swap_for:
    for (idx_t j = 0; j < 2 * SIZE; ++j)
    {
//#pragma HLS PIPELINE
#pragma HLS INLINE off
        float_t temp = AT(exmat, i1, j);
        AT(exmat, i1, j) = AT(exmat, i2, j);
        AT(exmat, i2, j) = temp;
    }
}

void add(float_t exmat[TOTALSIZE], idx_t i, idx_t j, float_t cj)
{
add_for:
    for (idx_t k = 0; k < 2 * SIZE; ++k)
//#pragma HLS PIPELINE
#pragma HLS INLINE off
    {
        float_t offset = cj * AT(exmat, j, k);
        AT(exmat, i, k) += offset;
    }
}

void mul(float_t exmat[TOTALSIZE], idx_t i, float_t ci)
{
mul_for:
    for (idx_t k = 0; k < 2 * SIZE; ++k)
    {
//#pragma HLS PIPELINE
#pragma HLS INLINE off
        AT(exmat, i, k) *= ci;
    }
}

void find_max_pivot_col(float_t exmat[TOTALSIZE], idx_t col, idx_t i_start, idx_t *i_best, float_t *val_best)
{
    idx_t i_max = NOT_FOUND;
    float_t val_max = -1;
fmp_for:
    // for (idx_t i = i_start; i < SIZE; ++i)
    for (idx_t i = 0; i < SIZE; ++i)
    {
        if (i < i_start)
        {
            continue;
        }
        float_t val = fabs(AT(exmat, i, col));
        // float_t val = fabs(exmat[i][col]);
        if (val > val_max)
        {
            //new max !!
            val_max = val;
            i_max = i;
        }
    }
    *i_best = i_max;
    *val_best = val_max;
}

void find_next_pivot_col(float_t exmat[TOTALSIZE], idx_t i_piv, idx_t j_piv, idx_t *i_next)
{
    idx_t i_next_piv;
    float_t valabs_next_piv;
    find_max_pivot_col(exmat, j_piv, i_piv, &i_next_piv, &valabs_next_piv);
    if (valabs_next_piv > 0)
    {
        // there is a pivot
        *i_next = i_next_piv;
    }
    else
    {
        // no valid pivot
        *i_next = NOT_FOUND;
    }
}

void gauss(float_t exmat[TOTALSIZE], idx_t *rank, float_t *determinant)
{

#if VERBOSE
    printf("extended initial matrix\n");
    print_exmat(SIZE, exmat);
#endif

    *determinant = 1;
    *rank = 0;
    idx_t i_piv = 0, j_piv = 0;
    idx_t i_piv_list[SIZE];
g_pivlist:
    for (idx_t k_piv_list = 0; k_piv_list < SIZE; ++k_piv_list)
    {
        i_piv_list[k_piv_list] = NOT_FOUND;
    }

g_global:
    // while (i_piv < SIZE && j_piv < SIZE)
    for (j_piv = 0; j_piv < SIZE; ++j_piv)
    {
        idx_t i_next;
        find_next_pivot_col(exmat, i_piv, j_piv, &i_next);
#if VERBOSE
        printf("-----STARTING STEP WITH i_piv=%d, j_piv=%d\n", i_piv, j_piv);
        printf("pivot search gave i=%d    \n", i_next);
#endif
        if (i_next == NOT_FOUND)
        {
            // means that no pivot are available for this column
            // start again at the next
            continue;
        }

        // j_piv = j_next;
        i_piv_list[j_piv] = i_piv;
        // 1 - swap
        swap(exmat, i_piv, i_next);
#if VERBOSE
        printf("after swap\n");
        print_exmat(SIZE, exmat);
#endif
        // 2 - scale
        float_t f_piv = AT(exmat, i_piv, j_piv);
        // float_t f_piv = exmat[i_piv][j_piv];
        mul(exmat, i_piv, 1 / f_piv);
        *determinant *= f_piv;
#if VERBOSE
        printf("after scale\n");
        print_exmat(SIZE, exmat);
#endif

    // 3 - eliminate
    g_eliminate:
        // for (idx_t i_line = i_piv + 1; i_line < SIZE; ++i_line)
        for (idx_t i_line = 0; i_line < SIZE; ++i_line)
        {
            if (i_line > i_piv)
            {

                float_t ci = -AT(exmat, i_line, j_piv);
                // float_t ci = -exmat[i_line][j_piv];
                add(exmat, i_line, i_piv, ci);
                // not strictly necessary but ensures numerical stability
                AT(exmat, i_line, j_piv) = 0;
                // exmat[i_line][j_piv] = 0;
            }
        }
        ++i_piv;
        // ++j_piv;
#if VERBOSE
        printf("after add\n");
        print_exmat(SIZE, exmat);
        printf("------END OF STEP ----------\n");
#endif
    }
#if VERBOSE
    printf("## END OF 1st PHASE\n");
    printf("## Pivot list\n");
    for (idx_t i = 0; i < SIZE; i++)
    {
        printf("%d ", i_piv_list[i]);
    }
    printf("\n");
#endif

// Reduce phase
g_reduce_j:
    for (idx_t j_rpiv = SIZE - 1; j_rpiv >= 0; --j_rpiv)
    {
        if (i_piv_list[j_rpiv] != NOT_FOUND)
        {
            // there was a pivot
            ++(*rank);
            idx_t i_rpiv = i_piv_list[j_rpiv];
#if VERBOSE
            printf("pivot found at %d, %d\n", i_rpiv, j_rpiv);
#endif
        g_reduct_i:
            // for (idx_t i_line = 0; i_line < i_rpiv; ++i_line)
            for (idx_t i_line = 0; i_line < SIZE; ++i_line)
            {
                if (i_line < i_rpiv)
                {

                    float_t ci = -AT(exmat, i_line, j_rpiv);
                    // float_t ci = -exmat[i_line][j_rpiv];
                    add(exmat, i_line, i_rpiv, ci);
                    // not strictly necessary but ensures numerical stability
                    AT(exmat, i_line, j_rpiv) = 0;
                    // exmat[i_line][j_rpiv] = 0;
                }
            }
        }
        else
        {
#if VERBOSE
            printf("no pivot\n");
#endif
        }
#if VERBOSE
        print_exmat(SIZE, exmat);
        printf("END OF REDUCE STEP\n");
#endif
    }

// determinant calculus
g_det:
    for (idx_t k = 0; k < SIZE; ++k)
    {
        *determinant *= AT(exmat, k, k);
        // *determinant *= exmat[k][k];
    }
#if VERBOSE
    printf("determinant is %lf\n", *determinant);

    // rank
    printf("rank is %d\n", *rank);

    printf("## END OF GAUSS\n");
#endif
}
