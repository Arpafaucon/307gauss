#ifndef GAUSS_H
#define GAUSS_H

#define SIZE 50
#define VERBOSE 0

#define EPS 1e-10
/**
 * @brief Swap two lines of the ext matrix
 */
void swap(double exmat[SIZE][2 * SIZE], int i, int i2);
/**
 * @brief Add one line to another, with a multiplicative coefficient
 */
void add(double exmat[SIZE][2 * SIZE], int i, int j, double cj);
/**
 * @brief multiply a row by a given coefficient
 * 
 * @param[in] exmat 
 * @param[in] i 
 * @param[in] ci 
 */
void mul(double exmat[SIZE][2 * SIZE], int i, double ci);
/**
 * @brief Search column by column the next suitable pivot
 * if the column is non-zero, returns the position of the max in absolute value
 * otherwise, go to next column
 * @param[out] i_next row of next pivot, -1 if none was found
 * @param[out] j_next 
 */
void find_next_pivot(double exmat[SIZE][2 * SIZE], int i_piv, int j_piv, int* i_next, int* j_next);

/**
 * @brief Total procedure, fill exmat with the extended reduced-eliminated version of inmat. Fill also rank and determinant
 */
void gauss(double inmat[SIZE][SIZE], double exmat[SIZE][2*SIZE], int* rank, double *determinant);

#endif