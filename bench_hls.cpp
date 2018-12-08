#include <time.h>
#include <string>

#include "bench_pc.h"
#include "gauss_fixed.h"
#include "bench_common.h"
#include "gauss_hls_stream.h"

namespace bench_hls
{
using namespace bench_common;

void compute_fixed_size_ip(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct)
{
#if TEST_DEBUG
    printf("test %d : fixed size\n on ip", test_id);
#endif
    double exmat[TOTALSIZE];
    clock_t total_start = clock(), total_end, core_start, core_end;
    double determinant;
    int rank;

    hls::stream<AXI_VALUE> in_stream;
    hls::stream<AXI_VALUE> out_stream;
    char control_port[GCP_BUS_SIZE] = {0};

    AXI_VALUE aValue;
    aValue.keep = -1;
    aValue.strb = -1;
    aValue.id = 0;
    aValue.user = 0;
    aValue.dest = 0;
    union stream_converter conv;

    // write input
    // matrix is sent as the upper left corner, rest is populated 
    // by identity matrix
    for (int i_row = 0; i_row < SIZE; i_row++)
    {
        for (int i_col = 0; i_col < SIZE; i_col++)
        {
            if (i_row < matrix_dim && i_col < matrix_dim)
            {
                double coeff = static_cast<double>(buffer[i_row * matrix_dim + i_col]);
                aValue.data = coeff;
            }
            else
            {
                aValue.data = (i_row == i_col);
            }
            aValue.last = (i_row == SIZE - 1) && (i_col == SIZE - 1);
            in_stream.write(aValue);
        }
    }

    core_start = clock();
    gauss_stream(in_stream, out_stream, control_port);
    core_end = clock();

    // read outputs
    // read determinant
    out_stream.read(aValue);
    conv.ap_val = aValue.data;
    determinant = conv.d_val;
    // read rank
    out_stream.read(aValue);
    rank = aValue.data;
    // read matrix output
    for (size_t i = 0; i < TOTALSIZE; i++)
    {
        out_stream.read(aValue);
        conv.ap_val = aValue.data;
        exmat[i] = conv.d_val;
    }

    total_end = clock();

    core_time = (double)(core_end - core_start) / CLOCKS_PER_SEC;
    total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    is_correct = 1;
    for (size_t i = 0; i < GCP_BUS_SIZE; i++)
    {
        is_correct *= control_port[i];
    }
    is_correct *= (rank == SIZE);

#if TEST_DEBUG
    print_exmat(exmat);
#endif
}

void dataset_battery()
{
    int num_tests, num_tests_corrects = 0, max_chars, max_coeffs;
    double total_time = 0, core_time = 0;
    int buffer[BUFFER_SIZE];
    int num_tests_ip = 0, num_tests_local = 0;
    clock_t global_start = clock(), global_end;

    // 0 - SETUP FILES
    ::FILE *f_dataset;
    char line[LINE_MAX];
    f_dataset = fopen(bench_common::dataset_fname.c_str(), "r");
    if (f_dataset == NULL)
    {
        printf("File %s not found\n", bench_common::dataset_fname.c_str());
        exit(1);
    }
    bench_common::read_header(f_dataset, num_tests, max_chars, max_coeffs);
    // sscanf(line, "%d,", &num_tests);
    printf("found %d tests\n", num_tests);

    for (int i = 0; i < num_tests; i++)
    {
        int test_id, matrix_dim, data_size;
        double core_time_taken = 0, total_time_taken = 0;
        int test_correct = 0;

        bench_common::read_test(f_dataset, BUFFER_SIZE, buffer, test_id, matrix_dim, data_size);

        if (matrix_dim <= SIZE)
        {
            ++num_tests_ip;
            compute_fixed_size_ip(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        }
        else
        {
            ++num_tests_local;
            bench_pc::compute_variable_size_pc(BUFFER_SIZE, buffer, matrix_dim, core_time_taken, total_time_taken, test_correct);
        }
        num_tests_corrects += test_correct;
        total_time += total_time_taken;
        core_time += core_time_taken;
    }

    global_end = clock();
    double global_time = (double)(global_end - global_start) / CLOCKS_PER_SEC;

    printf("## TEST COMPLETE ##\n");
    printf("#size  : %d\n", SIZE);
    printf("#tests  : %d (%d ip, %d local)\n", num_tests, num_tests_ip, num_tests_local);
    printf("total time   : %lf ms\t[file IO, test setup & comp.]\n", global_time * 1000);
    printf("test time    : %lf ms\t[test setup & comp.] \n", total_time * 1000);
    printf("core time    : %lf ms\t[test comp.] \n", core_time * 1000);
    printf("correct : %4lf [%d / %d] \n", 100. * num_tests_corrects / num_tests, num_tests_corrects, num_tests);
}

} // namespace bench_hls

int main(int argc, char **argv)
{
    printf("check long size : %d (should ==64bits)\n", 8 * sizeof(signed long));
    bench_hls::dataset_battery();
}
