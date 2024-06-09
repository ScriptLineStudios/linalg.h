#include <stdio.h>
#include <stdlib.h>

#define LINALG_IMPLEMENTATION
#include "linalg.h"

int main(void) {
    double data[25] = {
        2.0, 2.0, 3.0, 4.0, 5.0,
        5.0, 7.0, 7.0, 8.0, 5.0,
        9.0, 10.0, 12.0, 12.0, 5.0,
        13.0, 14.0, 19.0, 17.0, 5.0,
        13.0, 14.0, 15.0, 19.0, 5.0,
    };

    Matrix mat1 = matrix_new(5, 5, data);
    printf("determinant %f\n", matrix_determinant(&mat1));

    // printf("%f\n", matrix_determinant(&mat1));
    // Matrix cofactor = matrix_new(3, 3, NULL);
    // matrix_cofactor(&mat1, &cofactor, 2);
    // matrix_print(&cofactor);

    // // float data2[2] = {
    //     // 4.0, 
    //     // 5.0
    // // };
    // // matrix_print(&mat1);
    // // Matrix mat2 = matrix_subset(1, 3, 3, data+5);
    // Matrix mat2 =  matrix_subset(&mat1, 1, 3, 3, 1);
    // matrix_print(&mat2);
    // // printf("stride: %d index: %d\n", mat2.stride, MATRIX_INDEX(&mat2, 0, 3));


    // // Matrix mat2 = matrix_new(2, 1, data2);

    // // float *res = malloc(sizeof(float) * 2);
    // // Matrix result = matrix_new(2, 1, res);

    // // matrix_dot(&mat1, &mat2, &result);
    // // matrix_print(&result);

    return 0;
}