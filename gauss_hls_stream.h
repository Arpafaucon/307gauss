#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>

#define GCP_BUS_SIZE 5
// control part
#define GCP_STARTED 0
#define GCP_RECEIVED 1
#define GCP_COMPUTED 2
#define GCP_TRANSMITTED_META 3
#define GCP_TRANSMITTED_MATRIX 4

typedef ap_axis<64, 4, 5, 5> AXI_VALUE;
// that corresponds to int64
// id est signed long

union stream_converter {
    signed long ap_val;
    double d_val;
};

/**
 * @brief 
 * 
 * @param[in] in_stream input stream : SIZE*SIZE int64 in
 * @param[in] out_stream output stream : 2+TOTALSIZE 64bits out:
 *  - determinant (double)
 *  - rank (int64)
 *  - TOTALSIZE coeffs (double)
 */
void gauss_stream(hls::stream<AXI_VALUE> &in_stream, hls::stream<AXI_VALUE> &out_stream, char control_port[GCP_BUS_SIZE]);
