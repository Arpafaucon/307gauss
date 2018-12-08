#ifndef GAUSS_FIXED_H
#define GAUSS_FIXED_H
#include "gauss.h"


#define SIZE 20
#if SIZE > 120
#error SIZE should be containable in a char
#endif

#define TOTALSIZE 2 * SIZE *SIZE
typedef u_int8_t int8;
typedef double float64;
/**
 * @brief Swap two lines of the ext matrix
 */
void swap(float64 exmat[TOTALSIZE],  int8 i1, int8 i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(float64 exmat[TOTALSIZE], int8 i, int8 j, float64 cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(float64 exmat[TOTALSIZE], int8 i, float64 ci);

// /**
//  * @brief Search column by column the next suitable pivot
//  * if the column is non-zero, returns the position of the max in absolute value
//  * otherwise, go to next column
//  * @param[out] i_next row of next pivot, -1 if none was found
//  * @param[out] j_next
//  */
// void find_next_pivot(float64 exmat[TOTALSIZE], int i_piv, int j_piv, int* i_next, int* j_next);

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
void find_max_pivot_col(float64 exmat[TOTALSIZE], int8 col, int8 i_start, int8 *i_best, float64 *val_best);
/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(float64 exmat[TOTALSIZE], int8 *rank, float64 *determinant);

#endif
