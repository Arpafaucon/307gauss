#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "gauss.h"
#if VERBOSE
#include "main.h"
#endif

void swap(double exmat[SIZE][2 * SIZE], int i1, int i2)
{

swap_for:
    for (int j = 0; j < 2 * SIZE; ++j)
    {
        double temp = exmat[i1][j];
        exmat[i1][j] = exmat[i2][j];
        exmat[i2][j] = temp;
    }
}

void add(double exmat[SIZE][2 * SIZE], int i, int j, double cj)
{
add_for:
    for (int k = 0; k < 2 * SIZE; ++k)
    {
        exmat[i][k] += cj * exmat[j][k];
    }
}

void mul(double exmat[SIZE][2 * SIZE], int i, double ci)
{
mul_for:
    for (int k = 0; k < 2 * SIZE; ++k)
    {
        exmat[i][k] *= ci;
    }
}

void find_max_pivot_col(double exmat[SIZE][2 * SIZE], int col, int i_start, int *i_best, double *val_best)
{
    int i_max = -1;
    double val_max = -1;
fmp_for:
    for (int i = i_start; i < SIZE; ++i)
    {
        double val = fabs(exmat[i][col]);
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

void find_next_pivot_col(double exmat[SIZE][2 * SIZE], int i_piv, int j_piv, int *i_next)
{
    int i_next_piv;
    double valabs_next_piv;
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

#ifndef HLS
void find_next_pivot(double exmat[SIZE][2 * SIZE], int i_piv, int j_piv, int *i_next, int *j_next)
{
fnxtpiv_col_while:
    while (j_piv < SIZE)
    {
        int i_next_piv;
        double valabs_next_piv;
        find_max_pivot_col(exmat, j_piv, i_piv, &i_next_piv, &valabs_next_piv);
        if (valabs_next_piv > 0)
        {
            // there is a pivot
            *i_next = i_next_piv;
            *j_next = j_piv;
            return;
        }
        ++j_piv;
    }
    *i_next = -1;
    *j_next = -1;
}
#endif


void gauss(double inmat[SIZE][SIZE], double exmat[SIZE][2 * SIZE], int *rank, double *determinant)
{
#pragma HLS INTERFACE s_axilite port=determinant bundle=METRICS name=exdet
#pragma HLS INTERFACE s_axilite port=rank bundle=METRICS name=exrank
#pragma HLS INTERFACE m_axi depth=32 port=exmat name=exmat
#pragma HLS INTERFACE m_axi depth=32 port=inmat name=inmat
#pragma HLS INTERFACE ap_ctrl_none port=return
// first, extend inmat into exmat
g_ext_i:
    for (int i_ext = 0; i_ext < SIZE; ++i_ext)
    {
    g_ext_j:
        for (int j_ext = 0; j_ext < SIZE; ++j_ext)
        {
            // copy
            exmat[i_ext][j_ext] = inmat[i_ext][j_ext];
            // identity
            exmat[i_ext][SIZE + j_ext] = (i_ext == j_ext ? 1 : 0);
        }
    }
#if VERBOSE
    printf("extended initial matrix\n");
    print_exmat(exmat);
#endif

    *determinant = 1;
    *rank = 0;
    int i_piv = 0, j_piv = 0;
    int i_piv_list[SIZE];
g_pivlist:
    for (int k_piv_list = 0; k_piv_list < SIZE; ++k_piv_list)
    {
        i_piv_list[k_piv_list] = -1;
    }

g_global:
    // while (i_piv < SIZE && j_piv < SIZE)
    for(j_piv = 0; j_piv<SIZE; ++j_piv)
    {
        int i_next;
        find_next_pivot_col(exmat, i_piv, j_piv, &i_next);
#if VERBOSE
        printf("-----STARTING STEP WITH i_piv=%d, j_piv=%d\n", i_piv, j_piv);
        printf("pivot search gave %d    %d\n", i_next, j_next);
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
        print_exmat(exmat);
#endif
        // 2 - scale
        double f_piv = exmat[i_piv][j_piv];
        mul(exmat, i_piv, 1 / f_piv);
        *determinant *= f_piv;
#if VERBOSE
        printf("after scale\n");
        print_exmat(exmat);
#endif

    // 3 - eliminate
    g_eliminate:
        for (int i_line = i_piv + 1; i_line < SIZE; ++i_line)
        {
            double ci = -exmat[i_line][j_piv];
            add(exmat, i_line, i_piv, ci);
            // not strictly necessary but ensures numerical stability
            exmat[i_line][j_piv] = 0;
        }
        ++i_piv;
        // ++j_piv;
#if VERBOSE
        printf("after add\n");
        print_exmat(exmat);
        printf("------END OF STEP ----------\n");
#endif
    }
#if VERBOSE
    printf("## END OF 1st PHASE\n");
    printf("## Pivot list\n");
    for (int i = 0; i < SIZE; i++)
    {
        printf("%d ", i_piv_list[i]);
    }
    printf("\n");
#endif

// Reduce phase
g_reduce_j:
    for (int j_rpiv = SIZE - 1; j_rpiv >= 0; --j_rpiv)
    {
        if (i_piv_list[j_rpiv] >= 0)
        {
            // there was a pivot
            ++(*rank);
            int i_rpiv = i_piv_list[j_rpiv];
#if VERBOSE
            printf("pivot found at %d, %d\n", i_rpiv, j_rpiv);
#endif
        g_reduct_i:
            for (int i_line = 0; i_line < i_rpiv; ++i_line)
            {
                double ci = -exmat[i_line][j_rpiv];
                add(exmat, i_line, i_rpiv, ci);
                // not strictly necessary but ensures numerical stability
                exmat[i_line][j_rpiv] = 0;
            }
        }
        else
        {
#if VERBOSE
            printf("no pivot\n");
#endif
        }
#if VERBOSE
        print_exmat(exmat);
        printf("END OF REDUCE STEP\n");
#endif
    }

// determinant calculus
g_det:
    for (int k = 0; k < SIZE; ++k)
    {
        *determinant *= exmat[k][k];
    }
#if VERBOSE
    printf("determinant is %lf\n", *determinant);

    // rank
    printf("rank is %d\n", *rank);

    printf("## END OF GAUSS\n");
#endif
}
