#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "main.h"

void fill_matrix(double inmat[SIZE][SIZE], int const cap)
{

    for (int i = 0; i < SIZE; ++i)
    {
        for (int j = 0; j < SIZE; ++j)
        {
            int num = rand() % cap;
            inmat[i][j] = num;
        }
    }
}

double det_mat33(double inmat[SIZE][SIZE])
{
    assert(SIZE == 3);
    // loi de malus
    double det = 0;
    det += inmat[0][0] * inmat[1][1] * inmat[2][2];
    det += inmat[0][1] * inmat[1][2] * inmat[2][0];
    det += inmat[0][2] * inmat[1][0] * inmat[2][1];
    det -= inmat[0][2] * inmat[1][1] * inmat[2][0];
    det -= inmat[0][1] * inmat[1][0] * inmat[2][2];
    det -= inmat[0][0] * inmat[1][2] * inmat[2][1];
    return det;
}
// void mat_prod(double mat_a[SIZE][2*SIZE])

int test_matrix(double inmat[SIZE][SIZE], double exmat[SIZE][2 * SIZE], double determinant)
{
    int i, j;
#if TEST_VERBOSE
    printf("\n---- TEST -------\n");
    print_inmat(inmat);
    printf("\n");
    print_exmat(exmat);

    if (SIZE == 3)
    {
        double malus_det = det_mat33(inmat);
        printf("det33 : %lf\n", malus_det);
    }
#endif
    if (fabs(determinant) > EPS)
    {
#if TEST_VERBOSE
        printf("Is inversible (det = %lf)\n", determinant);
#endif
        // test if really inversible
        double prod1[SIZE][SIZE];
        double prod2[SIZE][SIZE];
        for (i = 0; i < SIZE; ++i)
        {
            for (j = 0; j < SIZE; ++j)
            {
                double temp1 = 0;
                double temp2 = 0;
                for (int k = 0; k < SIZE; k++)
                {
                    temp1 += exmat[i][SIZE + k] * inmat[k][j];
                    temp2 += inmat[i][k] * exmat[k][SIZE + j];
                }
                prod1[i][j] = temp1;
                prod2[i][j] = temp2;
            }
        }
#if TEST_VERBOSE
        printf("Product\n");
        print_inmat(prod1);
        print_inmat(prod2);
#endif
        double diff_total = 0;
        for (i = 0; i < SIZE; ++i)
        {
            for (j = 0; j < SIZE; ++j)
            {
                double expected = (i == j ? 1 : 0);
                double diff1 = fabs(expected - prod1[i][j]);
                double diff2 = fabs(expected - prod2[i][j]);
                diff_total += diff1;
                diff_total += diff2;
            }
        }
#if TEST_VERBOSE
        printf("Diff total : %lf\n", diff_total);
#endif
        if (diff_total < 1e-5)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    return -1;
}

void test_battery()
{
    clock_t cstart, cend;
    double inmat[SIZE][SIZE];
    double exmat[SIZE][2 * SIZE];
    int rank;
    double determinant;
    double total_time = 0;
    int test_total = 0, test_correct = 0;

    for (int i_test = 0; i_test < TEST_NUM; i_test++)
    {
        fill_matrix(inmat, TEST_CAP);
        cstart = clock();
        gauss(inmat, exmat, &rank, &determinant);
        // print_exmat(exmat);
        cend = clock();
        double time_taken = (double)(cend - cstart) / CLOCKS_PER_SEC;
        total_time += time_taken;
        int res = test_matrix(inmat, exmat, determinant);
        if (res > -1)
        {
            test_correct += res;
            test_total++;
        }
    }
    printf("## TEST COMPLETE ##\n");
    printf("#tests  : %d\n", TEST_NUM);
    printf("time    : %lf us (%d tests in %lf s )\n", 1.e6 * total_time / TEST_NUM, TEST_NUM, total_time);
    printf("correct : %4lf [%d / %d] \n", 100. * test_correct / test_total, test_correct, test_total);
}

void test_special()
{
    double inmat[SIZE][SIZE] = {{2, 6, 1}, {8, 7, 9}, {2, 0, 2}};
    double exmat[SIZE][2 * SIZE];
    double determinant;
    int rank;
    gauss(inmat, exmat, &rank, &determinant);
}

int main(int argc, char **argv)
{
    test_battery();
    // test_special();
}