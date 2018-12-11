#include "bench_pc.h"
#include <cstdio>

int main(int argc, char **argv)
{
    // test_battery();
    printf("%d %u \n", NOT_FOUND, NOT_FOUND);

    for (size_t i = 0; i < 1; i++)
    {
        printf("policy %d\n", i);
        bench_pc::dataset_battery("/home/arpad/dev/ensta/307/gauss/dataset/matrix0.txt", i);
        bench_pc::dataset_battery("/home/arpad/dev/ensta/307/gauss/dataset/matrix1.txt", i);
        bench_pc::dataset_battery("/home/arpad/dev/ensta/307/gauss/dataset/matrix2.txt", i);
        bench_pc::dataset_battery("/home/arpad/dev/ensta/307/gauss/dataset/matrix3.txt", i);
    }
    // test_special();
}