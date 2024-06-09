#ifndef LINALG_H
#define LINALG_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

typedef struct {
    bool allocated;
    double *data;
    int rows, cols;
    int stride;
} Matrix;

Matrix matrix_new(int rows, int cols, double *data);
Matrix matrix_with_stride_new(int stride, int rows, int cols, double *data);
Matrix matrix_subset(Matrix *mat, int start_row, int start_col, int rows, int cols); //take a subset of a matrix, ideally for reading only!

double *matrix_view(Matrix *mat, int row, int col);
double matrix_determinant(Matrix *mat);

void matrix_crout(Matrix *a, Matrix *l, Matrix *u);
void matrix_print(Matrix *mat);
void matrix_dot(Matrix *a, Matrix *b, Matrix *dst);
void matrix_cofactor(Matrix *mat, Matrix *dst, int cofactor);
void matrix_free(Matrix *mat);

#define MATRIX_AT(mat, row, col) *matrix_view((mat), (row), (col))
#define MATRIX_INDEX(mat, row, col) ((row) * (mat)->cols + (col)) + ((mat)->stride * (row)) // maybe not the traditional notion of a stride but whatever

#endif

#ifdef LINALG_IMPLEMENTATION

// thank you https://en.wikipedia.org/wiki/Crout_matrix_decomposition
void matrix_crout(Matrix *a, Matrix *l, Matrix *u) {
    double sum = 0.0f;
    // printf("%d %d %d %d %d %d\n", a->rows, l->rows, u->rows, a->cols, l->cols, u->cols);
    // assert((((((a->rows==l->rows)==u->rows)==a->cols)==l->cols)==u->cols));
    int n = a->rows;
    int i, j, k;
    
    for (i = 0; i < n; i++) {
        MATRIX_AT(u, i, j) = 1;
    }

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0.0f;
			for (k = 0; k < j; k++) {
				sum = sum + MATRIX_AT(l, i, k) * MATRIX_AT(u, k, j);	
			}
			MATRIX_AT(l, i, j) = MATRIX_AT(a, i, j) - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + MATRIX_AT(l, j, k) * MATRIX_AT(u, k, i);
			}
			// if (MATRIX_AT(l, j, j) == 0) {
				// printf("det(L) close to 0!\n Can't divide by 0...\n");
                // exit(1);
			// }
			MATRIX_AT(u, j, i) = (MATRIX_AT(a, j, i) - sum) / (MATRIX_AT(l, j, j)+0.0000000000001);
		}
	}
}

void matrix_cofactor(Matrix *mat, Matrix *dst, int cofactor) {
    assert(dst->rows == mat->rows-1);
    assert(dst->cols == mat->cols-1);
    //we also cross out the top row
    int index = 0;
    for (int row = 1; row < mat->rows; row++) {
        for (int col = 0; col < mat->cols; col++) {
            if (col != cofactor) {
                dst->data[index] = MATRIX_AT(mat, row, col);
                index++;
            }
        }
    }
}

double matrix_determinant(Matrix *mat) {
    Matrix L = matrix_new(mat->rows, mat->cols, NULL);
    Matrix U = matrix_new(mat->rows, mat->cols, NULL);
    matrix_crout(mat, &L, &U);

    matrix_print(&L);
    puts("");
    matrix_print(&U);
    puts("");

    double det = 1.0f; 
    for (int i = 0; i < L.rows; i++) {
        det *= MATRIX_AT(&L, i, i);
    }
    return det;
}

double matrix_determinant_naive(Matrix *mat) {
    assert(mat->rows == mat->cols);

    if (mat->rows == 2) {
        return (MATRIX_AT(mat, 0, 0) * MATRIX_AT(mat, 1, 1)) - (MATRIX_AT(mat, 0, 1) * MATRIX_AT(mat, 1, 0));
    } 
    
    double det = 0.0f;

    int row = 0; // we chose this randomly...
    for (int col = 0; col < mat->cols; col++) {
        Matrix cofactor = matrix_new(mat->rows-1, mat->cols-1, NULL);
        matrix_cofactor(mat, &cofactor, col);
        det += pow(-1, row+col) * MATRIX_AT(mat, row, col) * matrix_determinant(&cofactor);
        matrix_free(&cofactor);
    }
    return det;
}

Matrix matrix_new(int rows, int cols, double *data) {
    if (data) {
        return (Matrix){.allocated=false, .stride=0, .rows=rows, .cols=cols, .data=data};
    }
    double *doubles = (double *)malloc(sizeof(double) * rows * cols);
    return (Matrix){.allocated=true, .stride=0, .rows=rows, .cols=cols, .data=doubles};
}

void matrix_free(Matrix *mat) {
    if (mat->allocated) {
        free(mat->data);
    }
}

Matrix matrix_with_stride_new(int stride, int rows, int cols, double *data) {
    return (Matrix){.allocated=false, .stride=stride, .rows=rows, .cols=cols, .data=data};
}

Matrix matrix_subset(Matrix *mat, int start_row, int start_col, int rows, int cols) {
    return matrix_with_stride_new(mat->cols - cols, rows, cols, mat->data + MATRIX_INDEX(mat, start_row, start_col)); //hell yeah we can take subsets of subsets with this.
}

double *matrix_view(Matrix *mat, int row, int col) {
    return &(mat->data[MATRIX_INDEX(mat, row, col)]);
}

void matrix_print(Matrix *mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            printf("%f, ", MATRIX_AT(mat, i, j));
        }
        puts("");
    }
}

void matrix_dot(Matrix *a, Matrix *b, Matrix *dst) {
    assert(a->cols == b->rows);
    
    for (int row = 0; row < dst->rows; row++) {
        for (int col = 0; col < dst->cols; col++) {
            double accum = 0.0f;
            for (int i = 0; i < a->cols; i++) {
                // printf("a(%d %d): %f b(%d %d)%f\n", row, i, MATRIX_AT(a, row, i), i, col, MATRIX_AT(b, i, col));
                accum += MATRIX_AT(a, row, i) * MATRIX_AT(b, i, col);
            } 
            MATRIX_AT(dst, row, col) = accum;
        }
    }
}

#endif