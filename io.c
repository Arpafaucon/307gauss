#include <stdio.h>

#include "main.h"


void print_data(double data)
{
    if (data < EPS)
    {
        printf(" ____ ");
    }
    else
    {
        printf("%5.2lf ", data);
    }
}


void print_exmat(double exmat[SIZE][2 * SIZE])
{
    for (int i = 0; i < SIZE; ++i)
    {
        for (int j = 0; j < SIZE; ++j)
        {
            print_data(exmat[i][j]);
        }
        printf(" | ");
        for (int j = 0; j < SIZE; ++j)
        {
            print_data(exmat[i][SIZE + j]);
        }
        printf("\n");
    }
}

void print_inmat(double inmat[SIZE][SIZE])
{
    for (int i = 0; i < SIZE; ++i)
    {
        for (int j = 0; j < SIZE; ++j)
        {
            print_data(inmat[i][j]);
            // printf("%5.2lf ", inmat[i][j]);
        }
        printf("\n");
    }
}