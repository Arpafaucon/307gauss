#if !defined(BENCH_COMMON_H)
#define BENCH_COMMON_H

#include <stdio.h>
#include <string>

#include "gauss_fixed.h"

#define LINE_MAX 512
#define BUFFER_SIZE 12000
#define DATASET_UNIT 5

#define TEST_VERBOSE 0
#define TEST_DEBUG 0

namespace bench_common
{

std::string const dataset_fname = "/home/arpad/dev/ensta/307/gauss/dataset/matrix_dataset.txt";
std::string const golden_fname = "~/dev/ensta/307/gauss/dataset/golden.txt";

double &AT(int size, double *matptr, int i, int j);
double AT(int size, double const *matptr, int i, int j);
double &ATX(int size, double *matptr, int i, int j);
double ATX(int size, double const *matptr, int i, int j);

void gr_fgets(char *read_buffer, int n, FILE *input);

void print_matrix(int row, int col, int *matrix);
void read_header(FILE *f_dts, int &num_tests, int &max_chars, int &max_coeffs);
int read_test(FILE *f_dts, int buffer_size, int *buffer_data, int &test_id, int &matrix_size, int &data_size);

void fill_matrix(double inmat[SIZE][SIZE], int const cap);
int test_matrix(int size, double const *inmat, double const *exmat, double determinant);

} // namespace bench_common

#endif // BENCH_COMMON_H
