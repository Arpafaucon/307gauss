#ifndef GAUSS_FIXED_H
#define GAUSS_FIXED_H
#include "gauss.h"

#define SIZE 50
#if SIZE > 120
#error SIZE should be containable in a char
#endif

#define TOTALSIZE 2 * SIZE *SIZE

/**
 * @brief Swap two lines of the ext matrix
 */
void swap(float_t exmat[TOTALSIZE],  idx_t i1, idx_t i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(float_t exmat[TOTALSIZE], idx_t i, idx_t j, float_t cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(float_t exmat[TOTALSIZE], idx_t i, float_t ci);

// /**
//  * @brief Search column by column the next suitable pivot
//  * if the column is non-zero, returns the position of the max in absolute value
//  * otherwise, go to next column
//  * @param[out] i_next row of next pivot, -1 if none was found
//  * @param[out] j_next
//  */
// void find_next_pivot(float_t exmat[TOTALSIZE], int i_piv, int j_piv, int* i_next, int* j_next);

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
void find_max_pivot_col(float_t exmat[TOTALSIZE], idx_t col, idx_t i_start, idx_t *i_best, float_t *val_best);
/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(float_t exmat[TOTALSIZE], idx_t *rank, float_t *determinant);

#endif
