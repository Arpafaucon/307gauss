#include <cstdio>
#include <cmath>

#include "gauss.h"

inline float_t AT(idx_t size, float_t const *mat, idx_t i, idx_t j)
{
    return mat[i * size + j];
}

inline float_t ATX(idx_t size, float_t const *mat, idx_t i, idx_t j)
{
    return mat[i * 2 * size + j];
}

void print_data(float_t data)
{
    if (fabs(data) < EPS)
    {
        printf("  ..... ");
    }
    else
    {
        printf("%7.3lf ", data);
    }
}

void print_exmat(idx_t size, float_t const *exmat)
{
    for (idx_t i = 0; i < size; ++i)
    {
        for (idx_t j = 0; j < size; ++j)
        {
            float_t val = ATX(size, exmat, i, j);
            print_data(val);
        }
        printf(" | ");
        for (idx_t j = 0; j < size; ++j)
        {
            float_t val = ATX(size, exmat, i, size + j);
            print_data(val);
        }
        printf("\n");
    }
}

void print_inmat(idx_t size, float_t const *inmat)
{
    for (idx_t i = 0; i < size; ++i)
    {
        for (idx_t j = 0; j < size; ++j)
        {
            float_t val = AT(size, inmat, i, j);
            print_data(val);
            // printf("%5.2lf ", inmat[i][j]);
        }
        printf("\n");
    }
}