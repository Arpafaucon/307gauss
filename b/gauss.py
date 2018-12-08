#!/usr/bin/python3
#coding: utf8

import numpy as np


def swap(exmat, i: int, j: int):
    temp_row = exmat[i, :].copy()
    exmat[i, :] = exmat[j, :]
    exmat[j, :] = temp_row


def add(exmat, i: int, j: int, cj: float):
    """ Li <= Li + cj*Lj"""
    exmat[i, :] += cj * exmat[j, :]


def mul(exmat, i: int, ci: float):
    exmat[i, :] *= ci


def find_next_pivot(exmat, i_piv: int, j_piv: int, i_max: int, j_max: int):
    """
    find best pivot row in the submatrix starting from row i_piv and column j_piv
    returns i_best, j_best : best pivot for col j 
    if j_best>j_piv, then all intermediate columns have null coeffs
    if i_piv = -1, then no more pivot exist
    """
    while j_piv < j_max:
        piv_column = abs(exmat[i_piv:, j_piv])
        print('pivot in', piv_column)
        i_amax = np.argmax(piv_column)
        if piv_column[i_amax] > 0:
            # there is a pivot
            i_best = i_piv + i_amax
            j_best = j_piv
            return i_best, j_best
        j_piv += 1
    return -1, -1


def gauss(matrix: 'square matrix', size: int):
    ext_matrix = np.zeros((size, 2*size))
    ext_matrix[:, :size] = matrix
    ext_matrix[:, size:] = np.eye(size)
    print(ext_matrix)
    i_piv = 0
    j_piv = 0
    i_piv_list = [-1]*size
    while i_piv < size and j_piv < size:
        i_next, j_next = find_next_pivot(ext_matrix, i_piv, j_piv, size, size)
        if i_next == -1:
            break
        j_piv = j_next
        # j_piv_list.append(j)
        i_piv_list[j_piv] = i_piv
        print(i_next, j_piv)
        swap(ext_matrix, i_next, i_piv)
        f = ext_matrix[i_piv, j_piv]
        print('swap ', i_next, ' with ', i_piv)
        print(ext_matrix)
        print('pivot coeff is ', f)
        mul(ext_matrix, i_piv, 1/f)
        for i_line in range(i_piv+1, size):
            ci = - ext_matrix[i_line, j_piv]
            add(ext_matrix, i_line, i_piv, ci)
            ext_matrix[i_line, j_piv] = 0
        i_piv += 1
        j_piv += 1

        print(ext_matrix)

    print(i_piv_list)
    for nj, i_piv in enumerate(i_piv_list[::-1]):
        j_piv = size - 1 - nj
        print('inv. pivot is', (i_piv, j_piv), '  ', ext_matrix[i_piv, j_piv])
        for i_line in range(0, i_piv):
            ci = -ext_matrix[i_line, j_piv]
            add(ext_matrix, i_line, i_piv, ci)
            ext_matrix[i_line, j_piv] = 0
        print(ext_matrix)
    return ext_matrix


def test():
    pass


if __name__ == "__main__":
    np.set_printoptions(precision=2)
    t1 = np.array([[1, 2, 3, 4, 5], [2, 1, 3, 4, 5], [
                  2, 4, 7, 4, 4], [4, 1, -4, 4, 4], [1, -1, 0, 0, 0]])
    # t2 = np.zeros()
    # # res = find_next_pivot(t1, 0, 0, 3, 3)
    # print(t1)
    # print('-----------------')
    # res = gauss(t1, 3)
    # print('----------------')
    # print(res)
    exmat = gauss(t1, 5)
    red = exmat[:, :5]
    ext = exmat[:, 5:]
    a1 = np.matmul(t1, ext)
    print(a1)
    a2 = np.matmul(ext, t1)
    print(a2)

