// #include <cmath>

#include <ap_int.h>
#include <ap_axi_sdata.h>
#include <hls_stream.h>

// using namespace std;
#include "gauss_hls_stream.h"
#include "gauss_fixed.h"

#define COLS 2 * SIZE
#define AT(mat, i, j) (mat[(i) * (COLS) + (j)])

void gauss_stream(hls::stream<AXI_VALUE> &in_stream, hls::stream<AXI_VALUE> &out_stream, char control_port[GCP_BUS_SIZE])
{
#pragma HLS INTERFACE s_axilite port = return bundle = CONTROL_BUS
#pragma HLS INTERFACE axis port = out_stream name = OUTPUT_STREAM
#pragma HLS INTERFACE axis port = in_stream name = INPUT_STREAM
    double exmat[TOTALSIZE];
    int rank;
    double determinant;

    AXI_VALUE aValue;
    int i, j;
    // has started 
    control_port[GCP_STARTED] =1;
    // 
    // 1 - Read input stream
    // 
stream_read_A:
    for (i = 0; i < SIZE; ++i)
    {
    stream_read_B:
        for (j = 0; j < SIZE; ++j)
        {
            // read int64 value
            in_stream.read(aValue);
            double val = (double)aValue.data;
            AT(exmat, i, j) = val;
            AT(exmat, i, j + SIZE) = (i == j);
        }
    }
    control_port[GCP_RECEIVED]=1;
    // 
    // 2 - Compute
    // 
    gauss(exmat, &rank, &determinant);
    control_port[GCP_COMPUTED] = 1;
    // 
    // 3 - Write Metadata
    //     

    union stream_converter conv;
    // write determinant
    conv.d_val = determinant;
    aValue.data = conv.ap_val;
    aValue.last = 0;
    aValue.strb = -1;
    aValue.keep = -1; //e.strb;
    aValue.user = 0;
    aValue.id = 0;
    aValue.dest = 0;
    out_stream.write(aValue);
    // write rank - no conversion needed
    aValue.data = rank;
    aValue.last = 0;
    aValue.strb = -1;
    aValue.keep = -1; //e.strb;
    aValue.user = 0;
    aValue.id = 0;
    aValue.dest = 0;
    out_stream.write(aValue);
    control_port[GCP_TRANSMITTED_META] = 1;
    // 
    // 4 - Write output matrix
    //     
    // write matrix
stream_write_A:
    for (i = 0; i < TOTALSIZE; ++i)
    {
        conv.d_val = exmat[i];
        aValue.data = conv.ap_val;
        aValue.last = (i == TOTALSIZE - 1);
        aValue.strb = -1;
        aValue.keep = -1; //e.strb;
        aValue.user = 0;
        aValue.id = 0;
        aValue.dest = 0;
        //		aValue = OUT[i];
        out_stream.write(aValue);
    }
    control_port[GCP_TRANSMITTED_MATRIX] =1;
}
