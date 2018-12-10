#include <time.h>
#include <string>

#include "bench_pc.h"
#include "gauss_fixed_hls.h"
#include "bench_common.h"
#include "gauss_hls_stream.h"

namespace bench_hls
{
using namespace bench_common;

void fill_matrices_fixed(int const* buffer, idx_t matrix_dim, coeff_t *exmat, coeff_t* inmat)
{
        // Populate matrix for fixed SIZE_HLS
    for (int i_row = 0; i_row < SIZE_HLS; i_row++)
    {
        for (int i_col = 0; i_col < SIZE_HLS; i_col++)
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
            ATX(SIZE_HLS, exmat, i_row, i_col) = coeff;
            AT(SIZE_HLS, inmat, i_row, i_col) = coeff;
            ATX(SIZE_HLS, exmat, i_row, i_col + SIZE_HLS) = (i_row == i_col);

        }
    }
}

void compute_fixed_size_ip(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct)
{
#if TEST_DEBUG
    printf("test %d : fixed SIZE_HLS\n on ip", test_id);
#endif
    coeff_t exmat[TOTALSIZE_HLS];
    coeff_t inmat[SIZE_HLS * SIZE_HLS];
    clock_t total_start = clock(), total_end, core_start, core_end;
    det_t determinant;
    int rank;

    hls::stream<AXI_VALUE> in_stream;
    hls::stream<AXI_VALUE> out_stream;
    unsigned char control_port[GCP_BUS_SIZE] = {0};

    AXI_VALUE aValue;
    aValue.keep = -1;
    aValue.strb = -1;
    aValue.id = 0;
    aValue.user = 0;
    aValue.dest = 0;
    // union stream_converter conv;

    // write input
    // exmat is sent, extended by the Id matrix
    // by identity matrix
    fill_matrices_fixed(buffer, matrix_dim, exmat, inmat);

    for (size_t i = 0; i < TOTALSIZE_HLS; i++)
    {
        aValue.last = (i == TOTALSIZE_HLS - 1);

        coeff_t cval = exmat[i];
        aValue.data = cval.range();
        // printf("sending %lf\t%x\n", cval.to_float(), aValue.data.to_int());
        in_stream.write(aValue);
    }

    // for (int i_row = 0; i_row < SIZE_HLS; i_row++)
    // {
    //     for (int i_col = 0; i_col < SIZE_HLS; i_col++)
    //     {
    //         if (i_row < matrix_dim && i_col < matrix_dim)
    //         {
    //             double coeff = static_cast<double>(buffer[i_row * matrix_dim + i_col]);
    //             aValue.data = coeff;
    //         }
    //         else
    //         {
    //             aValue.data = (i_row == i_col);
    //         }
    //         aValue.last = (i_row == SIZE_HLS - 1) && (i_col == SIZE_HLS - 1);
    //         in_stream.write(aValue);
    //     }
    // }

    core_start = clock();
    gauss_stream(in_stream, out_stream, control_port);
    core_end = clock();

    // read outputs
    // read determinant
    // out_stream.read(aValue);
    // conv.ap_val = aValue.data;
    // determinant = conv.d_val;
    // // read rank
    // out_stream.read(aValue);
    // rank = aValue.data;
    rank = control_port[GCP_RESULT_RANK];
    ap_int<64>* apint_determinant_ptr = (ap_int<64>*) &control_port[GCP_RESULT_DETERMINANT];
    determinant.range() = *apint_determinant_ptr;
    // conv.ap_val = *apint_determinant_ptr;
    // determinant = conv.f_val;

    // read matrix output
    for (size_t i = 0; i < TOTALSIZE_HLS; i++)
    {
        coeff_t cval;
        out_stream.read(aValue);
        cval.range() = aValue.data;
        exmat[i] = cval;
    }

    total_end = clock();

    core_time = (double)(core_end - core_start) / CLOCKS_PER_SEC;
    total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    is_correct = 1;
    for (size_t i = 0; i < GCP_CHECKS_SIZE; i++)
    {
        is_correct *= (control_port[i]!=0);
    }
    is_correct *= (rank == SIZE_HLS);
    
    
     for(size_t i = 0; i < GCP_BUS_SIZE; i++)
     {
         printf("%3u ", control_port[i]);
     }
     printf("-> det=%lf\n", determinant.to_double());
    
    

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
    char line[BC_LINE_MAX];
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
//     for (int i = 0; i < 15; i++)
    {
        int test_id, data_size;
        double core_time_taken = 0, total_time_taken = 0;
        int test_correct = 0;
        idx_t matrix_dim;
        bench_common::read_test(f_dataset, BUFFER_SIZE, buffer, test_id, matrix_dim, data_size);

        if (matrix_dim <= SIZE_HLS)
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
    printf("#SIZE_HLS  : %d\n", SIZE_HLS);
    printf("#tests  : %d (%d ip, %d local)\n", num_tests, num_tests_ip, num_tests_local);
    printf("total time   : %lf ms\t[file IO, test setup & comp.]\n", global_time * 1000);
    printf("test time    : %lf ms\t[test setup & comp.] \n", total_time * 1000);
    printf("core time    : %lf ms\t[test comp.] \n", core_time * 1000);
    printf("correct : %4lf [%d / %d] \n", 100. * num_tests_corrects / num_tests, num_tests_corrects, num_tests);
}

} // namespace bench_hls

void test_ap_convert()
{
    double d = 12345.6789;
    coeff_t c_in = d;
    ap_int<32> ap_c = c_in.range();
    coeff_t c_out;
    c_out.range() = ap_c;
    double d_out = c_out;
    printf("check ap : %lf=%lf\n", d, d_out);
}

int main(int argc, char **argv)
{
    test_ap_convert();
    printf("check long SIZE_HLS : %d (should == %dbits)\n", 8 * sizeof(float_t), STREAM_WIDTH);
    bench_hls::dataset_battery();
}
