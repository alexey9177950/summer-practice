#ifndef MATRIX_IS_ALREADY_INCLUDED
#define MATRIX_IS_ALREADY_INCLUDED
#include <initializer_list>
#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <math.h>

class Matrix {
    std::vector<double> data_;
public:
    int n;
    int m;

    Matrix(int l_y = 0, int l_x = 0)
        : data_(l_x * l_y), n(l_y), m(l_x) {
    }

    Matrix(int l_y, int l_x, std::initializer_list<double> l)
        : data_(l_x * l_y), n(l_y), m(l_x) {
        int i = 0;
        for (double val : l) {
            data_[i++] = val;
        }
    }

    double& operator()(int i, int j) {
        return data_[i * m + j];
    }

    const double& operator()(int i, int j) const {
        return data_[i * m + j];
    }

    Matrix T() {
        Matrix ans(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                ans(j, i) = (*this)(i, j);
            }
        }
        return ans;
    }
};

Matrix operator*(const Matrix& A, const Matrix& B) {
    if (A.m != B.n) {
        throw std::exception();
    }
    Matrix C(A.n, B.m);
    for (int i = 0; i < C.n; ++i) {
        for (int j = 0; j < C.m; ++j) {
            C(i, j) = 0;
            for (int k = 0; k < A.m; ++k) {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }
    return C;
}

Matrix rand_matrix(int n, int m) {
    static std::mt19937 gen;
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    Matrix ans(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            ans(i, j) = dis(gen);
        }
    }
    return ans;
}

double norm_inf(const Matrix& A) {
    double ans = 0;
    for (int i = 0; i < A.n; ++i) {
        double val = 0;
        for (int j = 0; j < A.m; ++j) {
            val += fabs(A(i, j));
        }
        if (val > ans) {
            ans = val;
        }
    }
    return ans;
}

double norm_one(const Matrix& A) {
    double ans = 0;
    for (int j = 0; j < A.m; ++j) {
        double val = 0;
        for (int i = 0; i < A.n; ++i) {
            val += fabs(A(i, j));
        }
        if (val > ans) {
            ans = val;
        }
    }
    return ans;
}

double dist(const Matrix& A, const Matrix& B) {
    if (A.n != B.n || A.m != B.m) {
        throw std::exception();
    }
    int n = A.n, m = B.m;
    Matrix C(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            C(i, j) = A(i, j) - B(i, j);
        }
    }
    return norm_one(C);
}

std::ostream& operator <<(std::ostream& out, const Matrix& M) {
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.m; ++j) {
            out << M(i, j) << ' ';
        }
        std::cout << std::endl;
    }
    return out;
}

#endif
