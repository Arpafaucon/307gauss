#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <string>

#include "bench_common.h"
#include "io.h"
#include "gauss_fixed.h"
#include "gauss_var.h"

#include "bench_pc.h"

namespace bench_pc
{
using namespace bench_common;

void fill_matrices_fixed(int const* buffer, idx_t matrix_dim, float_t *exmat, float_t* inmat)
{
        // Populate matrix for fixed size
    for (int i_row = 0; i_row < SIZE; i_row++)
    {
        for (int i_col = 0; i_col < SIZE; i_col++)
        {
            float_t coeff;
            if (i_row < matrix_dim && i_col < matrix_dim)
            {
                int buf_coeff = buffer[i_row * matrix_dim + i_col];
                coeff = static_cast<float_t>(buf_coeff);
            }
            else
            {
                coeff = (i_row == i_col);
            }
            ATX(SIZE, exmat, i_row, i_col) = coeff;
            AT(SIZE, inmat, i_row, i_col) = coeff;
            ATX(SIZE, exmat, i_row, i_col + SIZE) = (i_row == i_col);

        }
    }
}

void compute_fixed_size_pc(int buffer_size, int const *buffer, idx_t matrix_dim, double &core_time, double &total_time, int &is_correct)
{
    float_t exmat[TOTALSIZE] = {0};
    float_t inmat[SIZE * SIZE] = {0};
    clock_t total_start = clock(), total_end, core_start, core_end;
    float_t determinant;
    idx_t rank;

#if TEST_DEBUG
    printf("test : fixed size\n");
#endif
    fill_matrices_fixed(buffer, matrix_dim, exmat, inmat);

    // 2 - COMPUTE
    core_start = clock();
    gauss(exmat, &rank, &determinant);
    core_end = clock();
    total_end = clock();

    // test

    int test_ok = test_matrix(SIZE, inmat, exmat, determinant);
    is_correct = (test_ok == 1);

    // chrono
    core_time = (double)(core_end - core_start) / CLOCKS_PER_SEC;
    total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    // is_correct = 1;
}

void compute_variable_size_pc(int buffer_size, int const *buffer, idx_t matrix_dim, double &core_time, double &total_time, int &is_correct)
{
#if TEST_DEBUG
    printf("test : variable size\n");
#endif
float_t *inmat = new float_t[matrix_dim*matrix_dim];
float_t *exmat = new float_t[2*matrix_dim*matrix_dim];
if (inmat == nullptr || exmat == nullptr)
{
    printf("Memory allocation failed\n");
    is_correct = 0;
    core_time = total_time = 0;
    return;
}

    // double *exmat = (double *)malloc(2 * matrix_dim * matrix_dim * sizeof(double));
    // double *inmat = (double *)malloc(matrix_dim * matrix_dim * sizeof(double));
    clock_t total_start = clock(), total_end, core_start, core_end;
    float_t determinant;
    idx_t rank;
    // populate matrix for variable size

    for (size_t i_row = 0; i_row < matrix_dim; i_row++)
    {
        for (size_t i_col = 0; i_col < matrix_dim; i_col++)
        {
            float_t coeff = (float_t)buffer[i_row * matrix_dim + i_col];
            ATX(matrix_dim, exmat, i_row, i_col) = coeff;
            AT(matrix_dim, inmat, i_row, i_col) = coeff;
            ATX(matrix_dim, exmat, i_row, i_col+matrix_dim) = (i_col == i_row);
        }
    }

    // 2 - COMPUTE
    core_start = clock();
    gvar::gauss(matrix_dim, exmat, &rank, &determinant);
    core_end = clock();
    total_end = clock();

    // test
    int test_ok = test_matrix(matrix_dim, inmat, exmat, determinant);
    is_correct = (test_ok == 1);

    core_time = (double)(core_end - core_start) / CLOCKS_PER_SEC;
    total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    // is_correct = 1;
    // free(exmat);
    // free(inmat);
    delete[] exmat;
    delete[] inmat;
}

// void mat_prod(double mat_a[SIZE][2*SIZE])
void dataset_battery(std::string dataset, int policy)
{
    int num_tests, num_tests_corrects = 0, max_chars, max_coeffs;
    double total_time = 0, core_time = 0;
    int buffer[BUFFER_SIZE];
    clock_t global_start = clock(), global_end;

    // 0 - SETUP FILES
    FILE *f_dataset;
    // char line[BC_LINE_MAX];
    f_dataset = fopen(dataset.c_str(), "r");
    if (f_dataset == NULL)
    {
        printf("File %s not found\n", dataset.c_str());
        exit(1);
    }
    read_header(f_dataset, num_tests, max_chars, max_coeffs);
    // sscanf(line, "%d,", &num_tests);
    // printf("found %d tests\n", num_tests);


    for (int i = 0; i < num_tests; i++)
    // for (int i = 0; i < 12; i++)
    {
        int test_id, data_size;
        idx_t matrix_dim;
        double core_time_taken, total_time_taken;
        int test_correct;

        int error = read_test(f_dataset, BUFFER_SIZE, buffer, test_id, matrix_dim, data_size);
        if(error){
            printf("skipping test %d\n", test_id);
            continue;
        }

        // if ( matrix_dim <= SIZE && (policy==1 ||  (policy == 2 && matrix_dim > SIZE/2)))
        // {
        //     compute_fixed_size_pc(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        // }
        // else
        // {
            // if(test_id == 11){
            // }

            compute_variable_size_pc(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        // }
        num_tests_corrects += test_correct;
        total_time += total_time_taken;
        core_time += core_time_taken;
        // break;
    }

    global_end = clock();
    double global_time = (double)(global_end - global_start) / CLOCKS_PER_SEC;

    printf("## TEST COMPLETE ##\n");
    // printf("#size  : %d\n", SIZE);
    printf("#tests  : %d\n", num_tests);
    printf("total time   : %lf ms\t[file IO, test setup & comp.]\n", global_time * 1000);
    printf("test time    : %lf ms\t[test setup & comp.] \n", total_time * 1000);
    printf("core time    : %lf ms\t[test comp.] \n", core_time * 1000);
    printf("correct : %4lf [%d / %d] \n", 100. * num_tests_corrects / num_tests, num_tests_corrects, num_tests);
}

} // namespace bench_pc

