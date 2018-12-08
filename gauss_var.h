#ifndef GAUSS_VAR_H
#define GAUSS_VAR_H

#include "gauss.h"

// #define AT(mat, i, j, size) (mat[(i)*(size)+(j)])
namespace gvar{

// for all below functions, exmat is assumed to be of size 2*size*size
/**
 * @brief Swap two lines of the ext matrix
 */
void swap(int size, double *exmat, int i, int i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(int size, double *exmat, int i, int j, double cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(int size, double *exmat, int i, double ci);

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
void find_max_pivot_col(int size, double *exmat, int col, int i_start, int *i_best, double *val_best);
/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(int size, double *exmat, int* rank, double *determinant);

}
#endif