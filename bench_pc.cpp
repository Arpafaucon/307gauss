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

void compute_fixed_size_pc(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct)
{
    double exmat[TOTALSIZE] = {0};
    double inmat[SIZE * SIZE] = {0};
    clock_t total_start = clock(), total_end, core_start, core_end;
    double determinant;
    int rank;

#if TEST_DEBUG
    printf("test : fixed size\n");
#endif
    // Populate matrix for fixed size
    for (int i_row = 0; i_row < SIZE; i_row++)
    {
        for (int i_col = 0; i_col < SIZE; i_col++)
        {
            double coeff;
            if (i_row < matrix_dim && i_col < matrix_dim)
            {
                int buf_coeff = buffer[i_row * matrix_dim + i_col];
                coeff = static_cast<double>(buf_coeff);
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
    // print_exmat(SIZE, exmat);

    // for(size_t i = 0; i < TOTALSIZE; i++)
    // {
    //     printf("%lf ", exmat[i]);
    // }

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

void compute_variable_size_pc(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct)
{
#if TEST_DEBUG
    printf("test : variable size\n");
#endif
double *inmat = new double[matrix_dim*matrix_dim];
double *exmat = new double[2*matrix_dim*matrix_dim];

    // double *exmat = (double *)malloc(2 * matrix_dim * matrix_dim * sizeof(double));
    // double *inmat = (double *)malloc(matrix_dim * matrix_dim * sizeof(double));
    clock_t total_start = clock(), total_end, core_start, core_end;
    double determinant;
    int rank;
    // populate matrix for variable size

    for (size_t i_row = 0; i_row < matrix_dim; i_row++)
    {
        for (size_t i_col = 0; i_col < matrix_dim; i_col++)
        {
            double coeff = (double)buffer[i_row * matrix_dim + i_col];
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
void dataset_battery()
{
    int num_tests, num_tests_corrects = 0, max_chars, max_coeffs;
    double total_time = 0, core_time = 0;
    int buffer[BUFFER_SIZE];
    clock_t global_start = clock(), global_end;

    // 0 - SETUP FILES
    FILE *f_dataset;
    char line[LINE_MAX];
    f_dataset = fopen(dataset_fname.c_str(), "r");
    if (f_dataset == NULL)
    {
        printf("File %s not found\n", dataset_fname.c_str());
        exit(1);
    }
    read_header(f_dataset, num_tests, max_chars, max_coeffs);
    // sscanf(line, "%d,", &num_tests);
    printf("found %d tests\n", num_tests);


    for (int i = 0; i < num_tests; i++)
    // for (int i = 0; i < 12; i++)
    {
        int test_id, matrix_dim, data_size;
        double core_time_taken, total_time_taken;
        int test_correct;

        read_test(f_dataset, BUFFER_SIZE, buffer, test_id, matrix_dim, data_size);

        if (matrix_dim <= SIZE)
        {
            compute_fixed_size_pc(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        }
        else
        {
            // if(test_id == 11){
            // }

            compute_variable_size_pc(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        }
        num_tests_corrects += test_correct;
        total_time += total_time_taken;
        core_time += core_time_taken;
        // break;
    }

    global_end = clock();
    double global_time = (double)(global_end - global_start) / CLOCKS_PER_SEC;

    printf("## TEST COMPLETE ##\n");
    printf("#size  : %d\n", SIZE);
    printf("#tests  : %d\n", num_tests);
    printf("total time   : %lf ms\t[file IO, test setup & comp.]\n", global_time * 1000);
    printf("test time    : %lf ms\t[test setup & comp.] \n", total_time * 1000);
    printf("core time    : %lf ms\t[test comp.] \n", core_time * 1000);
    printf("correct : %4lf [%d / %d] \n", 100. * num_tests_corrects / num_tests, num_tests_corrects, num_tests);
}

} // namespace bench_pc

// void test_battery()
// {
//     clock_t cstart, cend;
//     double inmat[SIZE][SIZE];
//     double exmat[SIZE][2 * SIZE];
//     int rank;
//     double determinant;
//     double total_time = 0;
//     int test_total = 0, test_correct = 0;

//     for (int i_test = 0; i_test < TEST_NUM; i_test++)
//     {
//         fill_matrix(inmat, TEST_CAP);

//         // copy phase
//         for (int i_ext = 0; i_ext < SIZE; ++i_ext)
//         {
//             for (int j_ext = 0; j_ext < SIZE; ++j_ext)
//             {
//                 // copy
//                 // AT(exmat, i_ext, j_ext) = AT(inmat, i_ext, j_ext);
//                 exmat[i_ext][j_ext] = inmat[i_ext][j_ext];

//                 // identity
//                 // AT(exmat, i_ext, SIZE + j_ext) = (i_ext == j_ext ? 1 : 0);
//                 exmat[i_ext][SIZE + j_ext] = (i_ext == j_ext ? 1 : 0);
//             }
//         }

//         cstart = clock();
//         gauss((double *)exmat, &rank, &determinant);
//         // print_exmat(exmat);
//         cend = clock();
//         double time_taken = (double)(cend - cstart) / CLOCKS_PER_SEC;
//         total_time += time_taken;
//         int res = test_matrix(inmat, exmat, determinant);
//         if (res > -1)
//         {
//             test_correct += res;
//             test_total++;
//         }
//     }
//     printf("## TEST COMPLETE ##\n");
//     printf("#size  : %d\n", SIZE);
//     printf("#tests  : %d\n", TEST_NUM);
//     printf("time    : %lf us (%d tests in %lf s )\n", 1.e6 * total_time / TEST_NUM, TEST_NUM, total_time);
//     printf("correct : %4lf [%d / %d] \n", 100. * test_correct / test_total, test_correct, test_total);
// }

// for (int i_test = 1; i_test <= num_tests; i_test++)
// {
//     // 1 - READ THE TEST DATA
//     int size, test_id;
//     double exmat[TOTALSIZE];
//     // 1.1 - Read data size
//     // strtol()
//     fscanf(f_dataset, "%d,%d,", &test_id, &size);
//     printf("test %d : size %d\n", test_id, size);

//     if (size > SIZE)
//     {
//         // test cannot be carried : ignore it
//     }

//     // 1.2 - Fill coefficients in SIZE*SIZE matrix
//     for (int i_row = 0; i_row < SIZE; i_row++)
//     {
//         for (int i_col = 0; i_col < SIZE; i_col++)
//         {
//             if (i_row < size && i_col < size)
//             {
//                 double coeff;
//                 fscanf(f_dataset, "%lf,", &coeff);
//                 printf("read %lf\n", coeff);
//                 AT(exmat, i_row, i_col) = coeff;
//             }
//             else
//             {
//                 AT(exmat, i_row, i_col) = (i_row == i_col) ? 1 : 0;
//             }
//             AT(exmat, i_row, i_col + SIZE) = (i_row == i_col) ? 1 : 0;
//         }
//     }
//     // print_exmat(exmat);

// }