// #include <cmath>

#include <ap_int.h>
#include <ap_axi_sdata.h>
#include <hls_stream.h>

// using namespace std;
#include "gauss_hls_stream.h"
#include "gauss_fixed_hls.h"

#define COLS 2 * SIZE_HLS
#define AT(mat, i, j) (mat[(i) * (COLS) + (j)])

void gauss_stream(hls::stream<AXI_VALUE> &in_stream, hls::stream<AXI_VALUE> &out_stream, unsigned char control_port[GCP_BUS_SIZE])
{
#pragma HLS INTERFACE ap_ctrl_hs port=return name=ipcontrol
#pragma HLS INTERFACE s_axilite port=return bundle=ctrl_bus name=fun_control

#pragma HLS INTERFACE axis register both port=out_stream name=stream_out
#pragma HLS INTERFACE axis register both port=in_stream name=stream_in
#pragma HLS INTERFACE s_axilite port=control_port bundle=ctrl_bus name=status_port
    coeff_t exmat[TOTALSIZE_HLS];
    hidx_t rank;
    det_t determinant;
    // union stream_converter conv;
    // union coeff_stream_converter coeff_conv;
    coeff_t c_val;
    AXI_VALUE aValue;
    hidx_t i;
    // int i, j;
    // has started
    control_port[GCP_STARTED] = 1;
//
// 1 - Read input stream
//
stream_read:
    for (i = 0; i < TOTALSIZE_HLS; ++i)
    {
        // read double value
        in_stream.read(aValue);
        
        c_val.range() = aValue.data;
        // printf("read %lf \t%x\n", c_val.to_float(), aValue.data.to_int());
        exmat[i] = c_val;
        control_port[GCP_STARTED] = i/2/SIZE_HLS;
    }
    control_port[GCP_RECEIVED] = 1;
    //
    // 2 - Compute
    //
    gauss(exmat, rank, determinant);
    // printf("rank = %d\n", rank);
    control_port[GCP_COMPUTED] = 1;
    //
    // 3 - Write Metadata
    //
    aValue.last = 0;
    aValue.strb = -1;
    aValue.keep = -1; //e.strb;
    aValue.user = 0;
    aValue.id = 0;
    aValue.dest = 0;

    // write determinant
    // conv.f_val = determinant;

    ap_int<64>* gcp_det_ptr = (ap_int<64>*)&control_port[GCP_RESULT_DETERMINANT];
    *gcp_det_ptr = determinant.range();
    // conv.d_val = determinant;
    // aValue.data = conv.ap_val;
    // out_stream.write(aValue);
    // write rank - no conversion needed
    control_port[GCP_RESULT_RANK] = rank;
    // aValue.data = rank;
    // out_stream.write(aValue);
    control_port[GCP_TRANSMITTED_META] = 1;
    //
    // 4 - Write output matrix
    //
    // write matrix
stream_write:
    for (i = 0; i < TOTALSIZE_HLS; ++i)
    {
        c_val = exmat[i];
        // conv.f_val = exmat[i];
        aValue.data = c_val.range();
        aValue.last = (i == TOTALSIZE_HLS - 1);
        control_port[GCP_TRANSMITTED_MATRIX] = i/2/SIZE_HLS;
        out_stream.write(aValue);

    }
//    control_port[GCP_TRANSMITTED_MATRIX] = 1;
}
