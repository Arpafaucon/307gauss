#if !defined(BENCH_COMMON_H)
#define BENCH_COMMON_H

#include <stdio.h>
#include <string>

#include "gauss_fixed.h"

// #define BC_LINE_MAX 512
#define BUFFER_SIZE 90000 //>200*200*2
#define DATASET_UNIT 5

#define TEST_VERBOSE 0
#define TEST_DEBUG 0

namespace bench_common
{

static constexpr int BC_LINE_MAX = 512;

template <typename T>
T &AT(idx_t size, T *matptr, idx_t i, idx_t j)
{
    idx_t index = size * i + j;
    return matptr[index];
}
template <typename T>
T AT(idx_t size, T const *matptr, idx_t i, idx_t j)
{
    idx_t index = size * i + j;
    return matptr[index];
}
template <typename T>
T &ATX(idx_t size, T *matptr, idx_t i, idx_t j)
{
    idx_t index = 2 * size * i + j;
    return matptr[index];
}
template <typename T>
T ATX(idx_t size, T const *matptr, idx_t i, idx_t j)
{
    idx_t index = 2 * size * i + j;
    return matptr[index];
}



std::string const dataset_fname = "/home/arpad/dev/ensta/307/gauss/dataset/100/matrix_dataset.txt";
std::string const golden_fname = "~/dev/ensta/307/gauss/dataset/golden.txt";

float_t &AT(idx_t size, float_t *matptr, idx_t i, idx_t j);
float_t AT(idx_t size, float_t const *matptr, idx_t i, idx_t j);
float_t &ATX(idx_t size, float_t *matptr, idx_t i, idx_t j);
float_t ATX(idx_t size, float_t const *matptr, idx_t i, idx_t j);

void gr_fgets(char *read_buffer, int n, FILE *input);

void print_matrix(idx_t row, idx_t col, idx_t *matrix);
void read_header(FILE *f_dts, int &num_tests, int &max_chars, int &max_coeffs);
int read_test(FILE *f_dts, int buffer_size, int *buffer_data, int &test_id, idx_t &matrix_size, int &data_size);

void fill_matrix(float_t inmat[SIZE][SIZE], idx_t const cap);
int test_matrix(idx_t size, float_t const *inmat, float_t const *exmat, float_t determinant);

} // namespace bench_common

#endif // BENCH_COMMON_H
