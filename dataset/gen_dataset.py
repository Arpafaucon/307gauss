#!/usr/bin/python3
# coding:utf8
import numpy as np
import random
import time

COEFF_MAX = 20
# DATASET = [(1, 5), (2, 2)]
DATASET = [(1, 5), (5, 5), (10, 50), (50, 100), (100, 5)]


def gen_matrix(size):
    mat = np.random.randint(-COEFF_MAX, COEFF_MAX, (size, size))
    # mat = np.random.uniform(-COEFF_MAX, COEFF_MAX, (size, size))
    det = int(np.linalg.det(mat))
    rank = np.linalg.matrix_rank(mat)
    return (size, (mat.reshape(-1)), (det, rank))

def write_gd_data(gold_file, test_id, det, rank):
    gold_file.write("{:<4},{:<10},{:<4}".format(test_id, det, rank))
    gold_file.write("\n")
    
def write_dts_data(dts_file, test_id, size, data):
    dts_file.write("{:<3},{:<3}".format(test_id, size))
    for data_val in data:
        dts_file.write(",{:<3}".format(data_val))
    dts_file.write("\n")


def header_data():
    zipped_dts = list(zip(*DATASET))
    max_size = max(zipped_dts[0])
    max_chars = 7 + 4*max_size
    sum_counts = sum(zipped_dts[1])

    time_tup = time.localtime(time.time())[1:5]
    header = "{:<3},{:<3},{:<3},{:<2}/{:<2},{:<2}:{:<2}\n".format(sum_counts, max_chars, max_size, *time_tup)
    print(header)
    return header

def gen_dataset(dts_filename, golden_filename):
    test_id = 0
    with open(dts_filename, 'w') as dts_f:
        with open(golden_filename, 'w') as gold_f:
            
            header = header_data()
            dts_f.write(header)
            gold_f.write(header)

            for size, count in DATASET:
                for i in range(count):
                    size, mat, golden = gen_matrix(size)
                    det, rank = golden
                    test_id += 1
                    write_dts_data(dts_f, test_id, size,  mat)
                    write_gd_data(gold_f, test_id, det, rank)
    print("dataset has {} entries".format(test_id))


line = list(gen_matrix(2))
print(line)
gen_dataset('matrix_dataset.txt', 'matrix_golden.txt')
