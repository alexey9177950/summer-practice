#ifndef SOLVE_ITER_INCLUDED
#define SOLVE_ITER_INCLUDED
#include "matrix.cpp"

// метод простой итерации
Matrix iter_method(Matrix A, Matrix B, int iter_num = 100, Matrix X = rand_matrix(0, 0)) {
    if (X.n == 0) {
        X = rand_matrix(A.n, 1);
    }
    for (int i = 0; i < A.n; ++i) {
        for (int j = 0; j < A.m; ++j) {
            A(i, j) *= -1.0;
        }
    }
    for (int i = 0; i < A.n; ++i) {
        A(i, i) += 1;
    }

    for (int iter = 0; iter < iter_num; ++iter) {
        Matrix new_X = A * X;
        for (int i = 0; i < X.n; ++i) {
            new_X(i, 0) += B(i, 0);
        }
        X = new_X;
    }
    return X;
}

// алгоритм Зейделя
Matrix zeidel(Matrix A, Matrix B, int iter_num = 100, Matrix p_X = rand_matrix(0, 0)) {
    if (p_X.n == 0) {
        p_X = rand_matrix(A.n, 1);
    }
    Matrix X(A.n, 1);
    for (int iter = 0; iter < iter_num; ++iter) {
        for (int i = 0; i < A.n; ++i) {
            X(i, 0) = B(i, 0);
            for (int j = 0; j < i; ++j) {
                X(i, 0) -= A(i, j) * X(j, 0);
            }
            for (int j = i + 1; j < A.n; ++j) {
                X(i, 0) -= A(i, j) * p_X(j, 0);
            }
            X(i, 0) /= A(i, i);
        }

        for (int i = 0; i < A.n; ++i) {
            p_X(i, 0) = X(i, 0);
        }
    }
    return X;
}

#endif
