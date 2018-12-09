

void test_special()
{
    double inmat[SIZE][SIZE] = {{2, 6, 1}, {8, 7, 9}, {2, 0, 2}};
    double exmat[SIZE][2 * SIZE];
    double determinant;
    int rank;
    // gauss(inmat, exmat, &rank, &determinant);
}

// void test_battery()
// {
//     clock_t cstart, cend;
//     double inmat[SIZE][SIZE];
//     double exmat[SIZE][2 * SIZE];
//     int rank;
//     double determinant;
//     double total_time = 0;
//     int test_total = 0, test_correct = 0;

//     for (int i_test = 0; i_test < TEST_NUM; i_test++)
//     {
//         fill_matrix(inmat, TEST_CAP);

//         // copy phase
//         for (int i_ext = 0; i_ext < SIZE; ++i_ext)
//         {
//             for (int j_ext = 0; j_ext < SIZE; ++j_ext)
//             {
//                 // copy
//                 // AT(exmat, i_ext, j_ext) = AT(inmat, i_ext, j_ext);
//                 exmat[i_ext][j_ext] = inmat[i_ext][j_ext];

//                 // identity
//                 // AT(exmat, i_ext, SIZE + j_ext) = (i_ext == j_ext ? 1 : 0);
//                 exmat[i_ext][SIZE + j_ext] = (i_ext == j_ext ? 1 : 0);
//             }
//         }

//         cstart = clock();
//         gauss((double *)exmat, &rank, &determinant);
//         // print_exmat(exmat);
//         cend = clock();
//         double time_taken = (double)(cend - cstart) / CLOCKS_PER_SEC;
//         total_time += time_taken;
//         int res = test_matrix(inmat, exmat, determinant);
//         if (res > -1)
//         {
//             test_correct += res;
//             test_total++;
//         }
//     }
//     printf("## TEST COMPLETE ##\n");
//     printf("#size  : %d\n", SIZE);
//     printf("#tests  : %d\n", TEST_NUM);
//     printf("time    : %lf us (%d tests in %lf s )\n", 1.e6 * total_time / TEST_NUM, TEST_NUM, total_time);
//     printf("correct : %4lf [%d / %d] \n", 100. * test_correct / test_total, test_correct, test_total);
// }

// for (int i_test = 1; i_test <= num_tests; i_test++)
// {
//     // 1 - READ THE TEST DATA
//     int size, test_id;
//     double exmat[TOTALSIZE];
//     // 1.1 - Read data size
//     // strtol()
//     fscanf(f_dataset, "%d,%d,", &test_id, &size);
//     printf("test %d : size %d\n", test_id, size);

//     if (size > SIZE)
//     {
//         // test cannot be carried : ignore it
//     }

//     // 1.2 - Fill coefficients in SIZE*SIZE matrix
//     for (int i_row = 0; i_row < SIZE; i_row++)
//     {
//         for (int i_col = 0; i_col < SIZE; i_col++)
//         {
//             if (i_row < size && i_col < size)
//             {
//                 double coeff;
//                 fscanf(f_dataset, "%lf,", &coeff);
//                 printf("read %lf\n", coeff);
//                 AT(exmat, i_row, i_col) = coeff;
//             }
//             else
//             {
//                 AT(exmat, i_row, i_col) = (i_row == i_col) ? 1 : 0;
//             }
//             AT(exmat, i_row, i_col + SIZE) = (i_row == i_col) ? 1 : 0;
//         }
//     }
//     // print_exmat(exmat);

// }