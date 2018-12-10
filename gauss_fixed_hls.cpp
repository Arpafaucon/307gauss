#include <cstdlib>
#include <hls_math.h>
// #include <cmath>

#include "gauss_fixed_hls.h"
#include "gauss.h"
#include <cstdio>
#define VERBOSE 0
// #if VERBOSE
// #include "io.h"
// #endif

#define COLS 2 * SIZE_HLS
#define AT(mat, i, j) (mat[(i) * (COLS) + (j)])

#if VERBOSE
void print_data(coeff_t data)
{
    if (hls::abs(data) < EPS)
    {
        printf("  ..... ");
    }
    else
    {
        printf("%7.3f ", data.to_float());
    }
}

void print_exmat(hidx_t size, coeff_t const *exmat)
{
    for (hidx_t i = 0; i < size; ++i)
    {
        for (hidx_t j = 0; j < size; ++j)
        {
            coeff_t val = AT(exmat, i, j);
            print_data(val);
        }
        printf(" | ");
        for (hidx_t j = 0; j < size; ++j)
        {
            coeff_t val = AT(exmat, i, size + j);
            print_data(val);
        }
        printf("\n");
    }
}
#endif


coeff_t coeff_abs(coeff_t val){
    return hls::abs(val);
    // return val;
}

void swap(coeff_t exmat[TOTALSIZE_HLS], hidx_t i1, hidx_t i2)
{
swap_for:
   for (hidx_t j = 0; j < 2 * SIZE_HLS; ++j)
   {
//#pragma HLS PIPELINE
       coeff_t temp = AT(exmat, i1, j);
       AT(exmat, i1, j) = AT(exmat, i2, j);
       AT(exmat, i2, j) = temp;
   }
}

void add(coeff_t exmat[TOTALSIZE_HLS], hidx_t i, hidx_t j, coeff_t cj)
{
add_for:
   for (hidx_t k = 0; k < 2 * SIZE_HLS; ++k)
//#pragma HLS PIPELINE
   {
       coeff_t offset = cj * AT(exmat, j, k);
       AT(exmat, i, k) += offset;
   }
}

void mul(coeff_t exmat[TOTALSIZE_HLS], hidx_t i, coeff_t ci)
{
//#pragma HLS INLINE off
mul_for:
   for (hidx_t k = 0; k < 2 * SIZE_HLS; ++k)
   {
//#pragma HLS PIPELINE
       AT(exmat, i, k) *= ci;
   }
}

void find_max_pivot_col(coeff_t exmat[TOTALSIZE_HLS], hidx_t col, hidx_t i_start, hidx_t &i_best, coeff_t &val_best)
{
    hidx_t i_max = NOT_FOUND;
    coeff_t val_max = -1;
fmp_for:
    // for (hidx_t i = i_start; i < SIZE_HLS; ++i)
    for (hidx_t i = 0; i < SIZE_HLS; ++i)
    {
        if (i < i_start)
        {
            continue;
        }
        coeff_t val = coeff_abs(AT(exmat, i, col));
        // coeff_t val = coeff_abs(exmat[i][col]);
        if (val > val_max)
        {
            //new max !!
            val_max = val;
            i_max = i;
        }
    }
    i_best = i_max;
    val_best = val_max;
}

void find_next_pivot_col(coeff_t exmat[TOTALSIZE_HLS], hidx_t i_piv, hidx_t j_piv, hidx_t &i_next)
{
    hidx_t i_next_piv;
    coeff_t valabs_next_piv;
    find_max_pivot_col(exmat, j_piv, i_piv, i_next_piv, valabs_next_piv);
    if (valabs_next_piv > 0)
    {
        // there is a pivot
        i_next = i_next_piv;
    }
    else
    {
        // no valid pivot
        i_next = NOT_FOUND;
    }
}

void gauss(coeff_t exmat[TOTALSIZE_HLS], hidx_t &rank, det_t &determinant)
{

#if VERBOSE
    printf("extended initial matrix\n");
    print_exmat(SIZE_HLS, exmat);
#endif

    determinant = 1.;
    rank = 0;
    hidx_t i_piv = 0, j_piv = 0;
    hidx_t i_piv_list[SIZE_HLS];
g_pivlist:
    for (hidx_t k_piv_list = 0; k_piv_list < SIZE_HLS; ++k_piv_list)
    {
        i_piv_list[k_piv_list] = NOT_FOUND;
    }

g_global:
    // while (i_piv < SIZE_HLS && j_piv < SIZE_HLS)
    for (j_piv = 0; j_piv < SIZE_HLS; ++j_piv)
    {
        hidx_t i_next = NOT_FOUND;
        find_next_pivot_col(exmat, i_piv, j_piv, i_next);
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
        print_exmat(SIZE_HLS, exmat);
#endif
        // 2 - scale
        coeff_t f_piv = AT(exmat, i_piv, j_piv);
        // coeff_t f_piv = exmat[i_piv][j_piv];
        const coeff_t one = 1.;
        mul(exmat, i_piv, one / f_piv);
        determinant *= f_piv;
#if VERBOSE
        printf("after scale\n");
        printf("scale %f, %lf\n", f_piv.to_float(), determinant.to_double());         print_exmat(SIZE_HLS, exmat);
#endif

    // 3 - eliminate
    g_eliminate:
        // for (hidx_t i_line = i_piv + 1; i_line < SIZE_HLS; ++i_line)
        for (hidx_t i_line = 0; i_line < SIZE_HLS; ++i_line)
        {
            if (i_line > i_piv)
            {
                coeff_t ci = -AT(exmat, i_line, j_piv);
                // coeff_t ci = -exmat[i_line][j_piv];
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
        print_exmat(SIZE_HLS, exmat);
        printf("------END OF STEP ----------\n");
#endif
    }
#if VERBOSE
    printf("## END OF 1st PHASE\n");
    printf("## Pivot list\n");
    for (hidx_t i = 0; i < SIZE_HLS; i++)
    {
        printf("%d ", i_piv_list[i]);
    }
    printf("\n");
#endif

// Reduce phase
g_reduce_j:
    for (hidx_t j_rpiv = SIZE_HLS ; j_rpiv -- > 0; )
    {
        if (i_piv_list[j_rpiv] != NOT_FOUND)
        {
            // there was a pivot
            ++rank;
            hidx_t i_rpiv = i_piv_list[j_rpiv];
#if VERBOSE
            printf("pivot found at %d, %d\n", i_rpiv, j_rpiv);
#endif
        g_reduct_i:
            // for (hidx_t i_line = 0; i_line < i_rpiv; ++i_line)
            for (hidx_t i_line = 0; i_line < SIZE_HLS; ++i_line)
            {
                if (i_line < i_rpiv)
                {

                    coeff_t ci = -AT(exmat, i_line, j_rpiv);
                    // coeff_t ci = -exmat[i_line][j_rpiv];
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
        print_exmat(SIZE_HLS, exmat);
        printf("END OF REDUCE STEP\n");
#endif
    }

// determinant calculus
g_det:
    for (hidx_t k = 0; k < SIZE_HLS; ++k)
    {
        determinant *= AT(exmat, k, k);
        // *determinant *= exmat[k][k];
    }
#if VERBOSE
    printf("determinant is %lf\n", determinant.to_double());

    // rank
    printf("rank is %d\n", rank);

    printf("## END OF GAUSS\n");
#endif
}
