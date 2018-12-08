// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gauss_var.h"
#if VERBOSE
#include "io.h"
#endif
namespace gvar
{

inline double& AT(double *mat, int i, int j, int size)
{
    return mat[i*size+j];
}

inline double& ATX(int size, double *mat, int i, int j)
{
    return mat[i*2*size+j];
}

void swap(int size, double *exmat, int i1, int i2)
{
swap_for:
    for (int j = 0; j < 2 * size; ++j)
    {
        double temp = ATX(size, exmat, i1, j);
        ATX(size, exmat, i1, j) = ATX(size, exmat, i2, j);
        ATX(size, exmat, i2, j) = temp;
        // double temp = exmat[i1][j];
        // exmat[i1][j] = exmat[i2][j];
        // exmat[i2][j] = temp;
    }
}

void add(int size, double *exmat, int i, int j, double cj)
{
add_for:
    for (int k = 0; k < 2 * size; ++k)
    {
        ATX(size, exmat, i, k) += cj * ATX(size, exmat, j, k);
        // exmat[i][k] += cj * exmat[j][k];
    }
}

void mul(int size, double *exmat, int i, double ci)
{
mul_for:
    for (int k = 0; k < 2 * size; ++k)
    {
        ATX(size, exmat, i, k) *= ci;
        // exmat[i][k] *= ci;
    }
}

void find_max_pivot_col(int size, double *exmat, int col, int i_start, int *i_best, double *val_best)
{
    int i_max = -1;
    double val_max = -1;
fmp_for:
    for (int i = i_start; i < size; ++i)
    {
        double val = fabs(ATX(size, exmat, i, col));
        // double val = fabs(exmat[i][col]);
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

void find_next_pivot_col(int size, double *exmat, int i_piv, int j_piv, int *i_next)
{
    int i_next_piv;
    double valabs_next_piv;
    find_max_pivot_col(size, exmat, j_piv, i_piv, &i_next_piv, &valabs_next_piv);
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

void gauss(int size, double *exmat, int *rank, double *determinant)
{
#if VERBOSE
    printf("extended initial matrix\n");
    print_exmat(size, exmat);
#endif

    *determinant = 1;
    *rank = 0;
    int i_piv = 0, j_piv = 0;

    // int *i_piv_list = (int *)malloc(size * sizeof(int));
    int *i_piv_list = new int[size];
g_pivlist:
    for (int k_piv_list = 0; k_piv_list < size; ++k_piv_list)
    {
        i_piv_list[k_piv_list] = -1;
    }

g_global:
    for (j_piv = 0; j_piv < size; ++j_piv)
    {
        int i_next;
        find_next_pivot_col(size, exmat, i_piv, j_piv, &i_next);
#if VERBOSE
        printf("-----STARTING STEP WITH i_piv=%d, j_piv=%d\n", i_piv, j_piv);
        printf("pivot search gave %d    %d\n", i_next, j_piv);
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
        swap(size, exmat, i_piv, i_next);
#if VERBOSE
        printf("after swap\n");
        print_exmat(size, exmat);
#endif
        // 2 - scale
        double f_piv = ATX(size, exmat, i_piv, j_piv);
        // double f_piv = exmat[i_piv][j_piv];
        mul(size, exmat, i_piv, 1 / f_piv);
        *determinant *= f_piv;
#if VERBOSE
        printf("after scale\n");
        print_exmat(size, exmat);
#endif

    // 3 - eliminate
    g_eliminate:
        for (int i_line = i_piv + 1; i_line < size; ++i_line)
        {
            double ci = -ATX(size, exmat, i_line, j_piv);
            // double ci = -exmat[i_line][j_piv];
            add(size, exmat, i_line, i_piv, ci);
            // not strictly necessary but ensures numerical stability
            ATX(size, exmat, i_line, j_piv) = 0;
            // exmat[i_line][j_piv] = 0;
        }
        ++i_piv;
        // ++j_piv;
#if VERBOSE
        printf("after add\n");
        print_exmat(size, exmat);
        printf("------END OF STEP ----------\n");
#endif
    }
#if VERBOSE
    printf("## END OF 1st PHASE\n");
    printf("## Pivot list\n");
    for (int i = 0; i < size; i++)
    {
        printf("%d ", i_piv_list[i]);
    }
    printf("\n");
#endif

// Reduce phase
g_reduce_j:
    for (int j_rpiv = size - 1; j_rpiv >= 0; --j_rpiv)
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
                double ci = -ATX(size, exmat, i_line, j_rpiv);
                // double ci = -exmat[i_line][j_rpiv];
                add(size, exmat, i_line, i_rpiv, ci);
                // not strictly necessary but ensures numerical stability
                ATX(size, exmat, i_line, j_rpiv) = 0;
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
        print_exmat(size, exmat);
        printf("END OF REDUCE STEP\n");
#endif
    }

// determinant calculus
g_det:
    for (int k = 0; k < size; ++k)
    {
        *determinant *= ATX(size, exmat, k, k);
        // *determinant *= exmat[k][k];
    }
#if VERBOSE
    printf("determinant is %lf\n", *determinant);

    // rank
    printf("rank is %d\n", *rank);

    printf("## END OF GAUSS\n");
#endif
    delete[] i_piv_list;
    // free(i_piv_list);
}

} // namespace gvar