#ifndef SOLVE_SYS_IS_ALREADY_INCLUDED
#define SOLVE_SYS_IS_ALREADY_INCLUDED
#include "matrix.cpp"
#include <math.h>

static double eps = 1e-10;

static void check_args(const Matrix A, Matrix B) {
    if (A.n != A.m || A.n != B.n || B.m != 1) {
        throw std::invalid_argument("Wrong size of matrices");
    }
}

static void swap_rows(Matrix& A, int i_1, int i_2) {
    for (int j = 0; j < A.m; ++j) {
        std::swap(A(i_1, j), A(i_2, j));
    }
}

static void swap_columns(Matrix& A, int j_1, int j_2) {
    for (int i = 0; i < A.n; ++i) {
        std::swap(A(i, j_1), A(i, j_2));
    }
}

// решение системы для верхнетреугольной матрицы
static Matrix get_ans(const Matrix& A, const Matrix& B) {
    Matrix X(B.n, 1);
    for (int i = X.n - 1; i >= 0; --i) {
        double x = B(i, 0);
        for (int j = X.n - 1; j > i; --j) {
            x -= A(i, j) * X(j, 0);
        }
        if (fabs(A(i, i)) > eps) {
            x /= A(i, i);
        }
        X(i, 0) = x;
    }
    return X;
}

// метод Гаусса
Matrix gauss(Matrix A, Matrix B) {
    check_args(A, B);

    for (int iter = 0; iter < A.m; ++iter) {
        // выбираем какую строку менять с текущей
        int to_swap = iter;
        while (to_swap < A.n && fabs(A(to_swap, iter)) < eps) {
            ++to_swap;
        }
        if (to_swap == A.n) {
            continue;
        }
        // меняем
        swap_rows(A, iter, to_swap);
        std::swap(B(iter, 0), B(to_swap, 0));
        // вычитаем из остальных строк с нужным коэффициентом
        int i_1 = iter;
        for (int i_2 = iter + 1; i_2 < A.n; ++i_2) {
            double k = A(i_2, iter) / A(i_1, iter);
            A(i_2, iter) = 0;
            for (int j = iter + 1; j < A.m; ++j) {
                A(i_2, j) -= k * A(i_1, j);
            }
            B(i_2, 0) -= k * B(i_1, 0);
        }
    }
    // A приведена к верхнетреугольному виду
    return get_ans(A, B); 
}

// выбор главного элемента
Matrix main_el(Matrix A, Matrix B) {
    check_args(A, B);

    for (int iter = 0; iter < A.m; ++iter) {
        // выбираем какую строку менять с текущей
        // (здесь единственное отличие от метода Гаусса)
        int to_swap = iter;
        for (int i = 0; i < A.n; ++i) {
            if (fabs(A(i, iter)) > fabs(A(to_swap, iter))) {
                to_swap = iter;
            }
        }
        if (fabs(A(to_swap, iter)) < eps) {
            continue;
        }
        // меняем
        swap_rows(A, iter, to_swap);
        std::swap(B(iter, 0), B(to_swap, 0));
        // вычитаем из остальных строк с нужным коэффициентом
        int i_1 = iter;
        for (int i_2 = iter + 1; i_2 < A.n; ++i_2) {
            double k = A(i_2, iter) / A(i_1, iter);
            A(i_2, iter) = 0;
            for (int j = iter + 1; j < A.m; ++j) {
                A(i_2, j) -= k * A(i_1, j);
            }
            B(i_2, 0) -= k * B(i_1, 0);
        }
    }
    // A приведена к верхнетреугольному виду
    return get_ans(A, B);
}

// выбор главного элемента по подматрице
Matrix main_el_2(Matrix A, Matrix B) {
    check_args(A, B);
    // перестановка столбцов:
    std::vector<int> perm(A.n);
    for (int i = 0; i < A.n; ++i) {
        perm[i] = i;
    }
    for (int iter = 0; iter < A.m; ++iter) {
        // выбираем главный элемент
        int i_sw = iter, j_sw = iter;
        for (int i = iter; i < A.n; ++i) {
            for (int j = iter; j < A.m; ++j) {
                if (fabs(A(i, j)) > fabs(A(i_sw, j_sw))) {
                    i_sw = i;
                    j_sw = j;
                }
            }
        }
        // меняем строки
        swap_rows(A, iter, i_sw);
        std::swap(B(iter, 0), B(i_sw, 0));
        // меняем столбцы
        swap_columns(A, iter, j_sw);
        std::swap(perm[iter], perm[j_sw]);
        // вычитаем из остальных строк с нужным коэффициентом
        int i_1 = iter;
        for (int i_2 = iter + 1; i_2 < A.n; ++i_2) {
            double k = A(i_2, iter) / A(i_1, iter);
            A(i_2, iter) = 0;
            for (int j = iter + 1; j < A.m; ++j) {
                A(i_2, j) -= k * A(i_1, j);
            }
            B(i_2, 0) -= k * B(i_1, 0);
        }
    }
    // A приведена к верхнетреугольному виду
    Matrix X_0(B.n, 1);
    for (int i = X_0.n - 1; i >= 0; --i) {
        double x = B(i, 0);
        for (int j = X_0.n - 1; j > i; --j) {
            x -= A(i, j) * X_0(j, 0);
        }
        if (fabs(A(i, i)) > eps) {
            x /= A(i, i);
        }
        X_0(i, 0) = x;
    }
    Matrix X(B.n, 1);
    for (int i = 0; i < B.n; ++i) {
        X(perm[i], 0) = X_0(i, 0);
    }
    return X;
}

static void rotate(Matrix& A, Matrix& B, int i_1, int i_2) {
    double x = A(i_1, i_1);
    double y = A(i_2, i_1);
    double denum = sqrtl(x * x + y * y);
    double k_x = x / denum;
    double k_y = - y / denum;
    for (int i = 0; i < A.n; ++i) {
        double x_1 = k_x * A(i_1, i) - k_y * A(i_2, i);
        double x_2 = k_y * A(i_1, i) + k_x * A(i_2, i);
        A(i_1, i) = x_1;
        A(i_2, i) = x_2;
    }
    double b_1 = k_x * B(i_1, 0) - k_y * B(i_2, 0);
    double b_2 = k_y * B(i_1, 0) + k_x * B(i_2, 0);
    B(i_1, 0) = b_1;
    B(i_2, 0) = b_2;
}

// метод вращений
Matrix rotation(Matrix A, Matrix B) {
    for (int i = 0; i < A.n; ++i) {
        for (int j = i + 1; j < A.n; ++j) {
            rotate(A, B, i, j);
        }
    }
    return get_ans(A, B);
}

#endif
