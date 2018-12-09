#ifndef GAUSS_VAR_H
#define GAUSS_VAR_H

#include "gauss.h"

// #define AT(mat, i, j, size) (mat[(i)*(size)+(j)])
namespace gvar{

// for all below functions, exmat is assumed to be of size 2*size*size
/**
 * @brief Swap two lines of the ext matrix
 */
void swap(idx_t size, float_t *exmat, idx_t i, idx_t i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(idx_t size, float_t *exmat, idx_t i, idx_t j, float_t cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(idx_t size, float_t *exmat, idx_t i, float_t ci);

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
void find_max_pivot_col(idx_t size, float_t *exmat, idx_t col, idx_t i_start, idx_t *i_best, float_t *val_best);
/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(idx_t size, float_t *exmat, idx_t* rank, float_t *determinant);

}
#endif