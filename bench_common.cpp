#include <cmath>
#include <cassert>
#include <string>
#include <cstdio>
#include <cstdlib>
// #include <stdio.h>
// #include <stdlib.h>

#include "gauss_fixed.h"
#include "bench_common.h"
#include "io.h"

namespace bench_common
{

double &AT(int size, double *matptr, int i, int j)
{
    int index = size * i + j;
    return matptr[index];
}
double AT(int size, double const *matptr, int i, int j)
{
    int index = size * i + j;
    return matptr[index];
}

double &ATX(int size, double *matptr, int i, int j)
{
    int index = (2 * size) * i + j;
    return matptr[index];
}
double ATX(int size, double const *matptr, int i, int j)
{
    int index = (2 * size) * i + j;
    return matptr[index];
}

void gr_fgets(char *read_buffer, int n, FILE *input)
{
    char *ptr = fgets(read_buffer, n, input);
#if TEST_DEBUG
    printf("    [R] : '%s'\n", read_buffer);
    if (ptr == nullptr)
    {
        printf("ERROR: null pointer in fgets. Exiting\n");
        exit(2);
    }
#endif
}

void print_matrix(int row, int col, int *matrix)
{
#if TEST_DEBUG
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%3d ", matrix[i * col + j]);
        }
        printf("\n");
    }
#endif
}

void read_header(FILE *f_dts, int &num_tests, int &max_chars, int &max_coeffs)
{
    char line[LINE_MAX];
    char *current_char = line, *dest_char;
    gr_fgets(line, LINE_MAX, f_dts);
    num_tests = strtol(current_char, &dest_char, 10);
    current_char = dest_char + 1;
#if TEST_DEBUG
    printf("%s\n", current_char);
#endif
    max_chars = strtol(current_char, &dest_char, 10);
    current_char = dest_char + 1;
#if TEST_DEBUG
    printf("%s\n", current_char);
#endif
    max_coeffs = strtol(current_char, &dest_char, 10);
    current_char = dest_char + 1;
#if TEST_DEBUG
    printf("%s\n", current_char);
#endif
    printf("num_tests=%d max_ch=%d max_cf=%d\n", num_tests, max_chars, max_coeffs);
}

/**
 * @brief 
 * 
 * @param[in] f_dts 
 * @param[in] buffer_size 
 * @param[in] buffer_data 
 * @param[out] test_id 
 * @param[out] matrix_size 
 * @param[out] data_size 
 * 
 * @return int 
 */
int read_test(FILE *f_dts, int buffer_size, int *buffer_data, int &test_id, int &matrix_size, int &data_size)
{
#if TEST_DEBUG
    printf("\n--reading test\n");
#endif
    char read_bf[DATASET_UNIT];
    char *current_char, *dest_char;
    // read test id
    gr_fgets(read_bf, DATASET_UNIT, f_dts);
    test_id = strtol(read_bf, nullptr, 10);
    // read matrix size
    gr_fgets(read_bf, DATASET_UNIT, f_dts);
    matrix_size = strtol(read_bf, nullptr, 10);

    data_size = matrix_size * matrix_size;
    if (data_size > buffer_size)
    {
        printf("ERROR: Buffer too small\n");
        // exhaust line until new line
        char car;
        do
        {
            car = fgetc(f_dts);
        } while (car != '\n' || car == EOF);

        return 1;
    }

    for (int i = 0; i < data_size; ++i)
    {
        // read one coeff
        gr_fgets(read_bf, DATASET_UNIT, f_dts);
        buffer_data[i] = strtol(read_bf, nullptr, 10);
    }
#if TEST_DEBUG
    printf("read test: test_id=%d, m_size=%d, data_size=%d \n", test_id, matrix_size, data_size);
    print_matrix(matrix_size, matrix_size, buffer_data);
#endif
#if TEST_DEBUG
    char exp_eof = read_bf[DATASET_UNIT - 2];
    printf("EOL check : %d\n", exp_eof == '\n');
#endif
    return 0;
}

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

int test_matrix(int size, double const *inmat, double const *exmat, double determinant)
{
    int i, j;
#if TEST_VERBOSE
    printf("\n---- TEST -------\n");
    print_inmat(size, inmat);
    printf("\n");
    print_exmat(size, exmat);
#endif
    if (fabs(determinant) < EPS)
    {
        return -1;
    }
#if TEST_VERBOSE
    printf("Is inversible (det = %lf)\n", determinant);
#endif
    // test if really inversible
    double *prod1 = nullptr;
    double *prod2 = nullptr;
    prod1 = new double[size * size];
    prod2 = new double[size * size];
    // double *prod1 = (double *)malloc(size * size * sizeof(double));
    // double *prod2 = (double *)malloc(size * size * sizeof(double));

    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            double temp1 = 0;
            double temp2 = 0;
            for (int k = 0; k < size; k++)
            {
                temp1 += ATX(size, exmat, i, size + k) * AT(size, inmat, k, j);
                temp2 += AT(size, inmat, i, k) * ATX(size, exmat, k, size + j);
            }
            AT(size, prod1, i, j) = temp1;
            AT(size, prod2, i, j) = temp2;
            // prod1[i][j] = temp1;
            // prod2[i][j] = temp2;
        }
    }
#if TEST_VERBOSE
    printf("Product\n");
    print_inmat(size, prod1);
    printf("--\n");
    print_inmat(size, prod2);
#endif
    double diff_total = 0;
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            double expected = (i == j ? 1 : 0);
            double diff1 = fabs(expected - AT(size, prod1, i, j));
            double diff2 = fabs(expected - AT(size, prod2, i, j));
            // double diff2 = fabs(expected - prod2[i][j]);
            diff_total += diff1;
            diff_total += diff2;
        }
    }
#if TEST_VERBOSE
    printf("Diff total : %lf\n", diff_total);
#endif
    delete[] prod1;
    delete[] prod2;
    if (diff_total < 1e-5)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

} // namespace bench_common