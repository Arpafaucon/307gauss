// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gauss_fixed.h"
#if VERBOSE
#include "io.h"
#endif

#define COLS 2 * SIZE
#define AT(mat, i, j) (mat[(i) * (COLS) + (j)])


void swap(float64 exmat[TOTALSIZE], int8 i1, int8 i2)
{

swap_for:
    for (int8 j = 0; j < 2 * SIZE; ++j)
    {
        float64 temp = AT(exmat, i1, j);
        AT(exmat, i1, j) = AT(exmat, i2, j);
        AT(exmat, i2, j) = temp;
        // float64 temp = exmat[i1][j];
        // exmat[i1][j] = exmat[i2][j];
        // exmat[i2][j] = temp;
    }
}

void add(float64 exmat[TOTALSIZE], int8 i, int8 j, float64 cj)
{
add_for:
    for (int8 k = 0; k < 2 * SIZE; ++k)
    {
        AT(exmat, i, k) += cj * AT(exmat, j, k);
        // exmat[i][k] += cj * exmat[j][k];
    }
}

void mul(float64 exmat[TOTALSIZE], int8 i, float64 ci)
{
mul_for:
    for (int8 k = 0; k < 2 * SIZE; ++k)
    {
        AT(exmat, i, k) *= ci;
        // exmat[i][k] *= ci;
    }
}

void find_max_pivot_col(float64 exmat[TOTALSIZE], int8 col, int8 i_start, int8 *i_best, float64 *val_best)
{
    int8 i_max = -1;
    float64 val_max = -1;
fmp_for:
    for (int8 i = i_start; i < SIZE; ++i)
    {
        float64 val = fabs(AT(exmat, i, col));
        // float64 val = fabs(exmat[i][col]);
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

void find_next_pivot_col(float64 exmat[TOTALSIZE], int8 i_piv, int8 j_piv, int8 *i_next)
{
    int8 i_next_piv;
    float64 valabs_next_piv;
    find_max_pivot_col(exmat, j_piv, i_piv, &i_next_piv, &valabs_next_piv);
    if (valabs_next_piv > 0)
    {
        // there is a pivot
        *i_next = i_next_piv;
    }
    else
    {
        // no valid pivot
        *i_next = -1;
    }
}

void gauss(float64 exmat[TOTALSIZE], int8 *rank, float64 *determinant)
{

#if VERBOSE
    printf("extended initial matrix\n");
    print_exmat(SIZE, exmat);
#endif

    *determinant = 1;
    *rank = 0;
    int8 i_piv = 0, j_piv = 0;
    int8 i_piv_list[SIZE];
g_pivlist:
    for (int8 k_piv_list = 0; k_piv_list < SIZE; ++k_piv_list)
    {
        i_piv_list[k_piv_list] = -1;
    }

g_global:
    // while (i_piv < SIZE && j_piv < SIZE)
    for (j_piv = 0; j_piv < SIZE; ++j_piv)
    {
        int8 i_next;
        find_next_pivot_col(exmat, i_piv, j_piv, &i_next);
#if VERBOSE
        printf("-----STARTING STEP WITH i_piv=%d, j_piv=%d\n", i_piv, j_piv);
        printf("pivot search gave i=%d    \n", i_next);
#endif
        if (i_next == -1)
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
        float64 f_piv = AT(exmat, i_piv, j_piv);
        // float64 f_piv = exmat[i_piv][j_piv];
        mul(exmat, i_piv, 1 / f_piv);
        *determinant *= f_piv;
#if VERBOSE
        printf("after scale\n");
        print_exmat(SIZE, exmat);
#endif

    // 3 - eliminate
    g_eliminate:
        for (int8 i_line = i_piv + 1; i_line < SIZE; ++i_line)
        {
            float64 ci = -AT(exmat, i_line, j_piv);
            // float64 ci = -exmat[i_line][j_piv];
            add(exmat, i_line, i_piv, ci);
            // not strictly necessary but ensures numerical stability
            AT(exmat, i_line, j_piv) = 0;
            // exmat[i_line][j_piv] = 0;
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
    for (int8 i = 0; i < SIZE; i++)
    {
        printf("%d ", i_piv_list[i]);
    }
    printf("\n");
#endif

// Reduce phase
g_reduce_j:
    for (int8 j_rpiv = SIZE - 1; j_rpiv >= 0; --j_rpiv)
    {
        if (i_piv_list[j_rpiv] >= 0)
        {
            // there was a pivot
            ++(*rank);
            int8 i_rpiv = i_piv_list[j_rpiv];
#if VERBOSE
            printf("pivot found at %d, %d\n", i_rpiv, j_rpiv);
#endif
        g_reduct_i:
            for (int8 i_line = 0; i_line < i_rpiv; ++i_line)
            {
                float64 ci = -AT(exmat, i_line, j_rpiv);
                // float64 ci = -exmat[i_line][j_rpiv];
                add(exmat, i_line, i_rpiv, ci);
                // not strictly necessary but ensures numerical stability
                AT(exmat, i_line, j_rpiv) = 0;
                // exmat[i_line][j_rpiv] = 0;
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
    for (int8 k = 0; k < SIZE; ++k)
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
