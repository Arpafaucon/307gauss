
#if !defined(GAUSS_HLS_STREAM_H)
#define GAUSS_HLS_STREAM_H





#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>

#define GCP_BUS_SIZE 16
#define GCP_CHECKS_SIZE 5
// control part
#define GCP_STARTED 0
#define GCP_RECEIVED 1
#define GCP_COMPUTED 2
#define GCP_TRANSMITTED_META 3
#define GCP_TRANSMITTED_MATRIX 4
#define GCP_RESULT_RANK 7
#define GCP_RESULT_DETERMINANT 8 // 8 bytes 8->15

// should be sizeof(float_t)
#define STREAM_WIDTH 32

typedef ap_axis<STREAM_WIDTH, 4, 5, 5> AXI_VALUE;
// that corresponds to int64
// id est signed long

// union stream_converter {
//     signed int ap_val;
//     float_t f_val;
// };

// union coeff_stream_converter {
//     ap_int<STREAM_WIDTH> ap_val;
//     coeff_t c_val;
// };

// union det_stream_converter {
//     ap_int<64> ap_val;
//     det_t det_val;
// };



/**
 * @brief 
 * 
 * @param[in] in_stream input stream : TOTALSIZE float in
 * @param[in] out_stream output stream : 2+TOTALSIZE 64bits out:
 *  - determinant (double)
 *  - rank (int64)
 *  - TOTALSIZE coeffs (double)
 */
void gauss_stream(hls::stream<AXI_VALUE> &in_stream, hls::stream<AXI_VALUE> &out_stream, unsigned char control_port[GCP_BUS_SIZE]);

#endif // GAUSS_HLS_STREAM_H