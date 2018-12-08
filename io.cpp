#include <cstdio>
#include <cmath>

#include "gauss.h"

// #include "bench.h"
// #include "gauss.h"

inline double AT(int size, double const *mat, int i, int j)
{
    return mat[i * size + j];
}

inline double ATX(int size, double const *mat, int i, int j)
{
    return mat[i * 2 * size + j];
}

void print_data(double data)
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

void print_exmat(int size, double const *exmat)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            double val = ATX(size, exmat, i, j);
            print_data(val);
        }
        printf(" | ");
        for (int j = 0; j < size; ++j)
        {
            double val = ATX(size, exmat, i, size + j);
            print_data(val);
        }
        printf("\n");
    }
}

void print_inmat(int size, double const *inmat)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            double val = AT(size, inmat, i, j);
            print_data(val);
            // printf("%5.2lf ", inmat[i][j]);
        }
        printf("\n");
    }
}