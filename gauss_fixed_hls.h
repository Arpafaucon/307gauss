#if !defined(GAUSS_FIXED_HLS_H)
#define GAUSS_FIXED_HLS_H

#include <ap_fixed.h>
//#include "gauss.h"s


#define SIZE_HLS 100
#define TOTALSIZE_HLS 2 * SIZE_HLS *SIZE_HLS
#define IDX_WIDTH 18 //should have 2^IDX_WIDTH > TOTALSIZE

typedef ap_fixed<32, 16> coeff_t;
typedef ap_fixed<64, 32> det_t;
typedef ap_uint<IDX_WIDTH> hidx_t;
// typedef short int hidx_t;

/**
 * @brief Swap two lines of the ext matrix
 */
void swap(coeff_t exmat[TOTALSIZE_HLS],  hidx_t i1, hidx_t i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(coeff_t exmat[TOTALSIZE_HLS], hidx_t i, hidx_t j, coeff_t cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(coeff_t exmat[TOTALSIZE_HLS], hidx_t i, coeff_t ci);

// /**
//  * @brief Search column by column the next suitable pivot
//  * if the column is non-zero, returns the position of the max in absolute value
//  * otherwise, go to next column
//  * @param[out] i_next row of next pivot, -1 if none was found
//  * @param[out] j_next
//  */
// void find_next_pivot(coeff_t exmat[TOTALSIZE_HLS], int i_piv, int j_piv, int* i_next, int* j_next);

/**
 * @brief Find the best pivot in a column
 * 
 * @param[in] exmat 
 * @param[in] col 
 * @param[in] i_start 
 * @param[in] i_best 
 * @param[in] val_best 
 * 
 */
void find_max_pivot_col(coeff_t exmat[TOTALSIZE_HLS], hidx_t col, hidx_t i_start, hidx_t &i_best, coeff_t &val_best);
/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(coeff_t exmat[TOTALSIZE_HLS], hidx_t &rank, det_t &determinant);



#endif // GAUSS_FIXED_HLS_H
