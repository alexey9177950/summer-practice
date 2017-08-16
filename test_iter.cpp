#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.cpp"
#include "solve_iter.cpp"
#include "solve_sys.cpp"

using iter_m = Matrix(Matrix, Matrix, int, Matrix);
using std::vector;

Matrix rand_m_norm_one(int n, int m, double val = 0.5) {
    Matrix ans = rand_matrix(n, m);
    double k = val / norm_one(ans);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ans(i, j) *= k;
        }
    }
    return ans;
}

void test_convergention(std::string adr, const vector<int>& iter_n, iter_m method) {
	vector<int> dim{5, 10, 50, 200};
    std::ofstream file(adr + "_data.txt");
    for (int i : dim) {
        file << i << ' ';
    }
    file << std::endl;
    for (int i : iter_n) {
        file << i << ' ';
    }
    file << std::endl;
    for (int n : dim) {
        Matrix true_X = rand_matrix(n, 1);
        Matrix A = rand_m_norm_one(n, n, 0.6);
        for (int i = 0; i < A.n; ++i) {
            A(i, i) += 1;
        }
        Matrix B = A * true_X;
        Matrix X = method(A, B, iter_n.front(), B);
        file << dist(X, true_X) << ' ';
        for (size_t i = 1; i < iter_n.size(); ++i) {
            X = iter_method(A, B, iter_n[i] - iter_n[i - 1], X);
            file << dist(X, true_X) << ' ';
        }
        file << std::endl;
    }
}

int iter_num(iter_m method, const Matrix& A, const Matrix X_0) {
    int ans = 1;
    Matrix B = A * X_0;
    Matrix X = method(A, B, 1, rand_matrix(A.n, 1));
    while (dist(X, X_0) / norm_one(X_0) > 1e-10) {
        X = method(A, B, 1, X);
        ++ans;
    }
    return ans;
}

void test_iter_num() {
    int n = 200;
    int att_num = 20;
    std::ofstream out("./iter_num_data.txt");
    for (double k = 0.1; k < 16; k += 0.5) {
        double mean = 0;
        for (int i = 0; i < att_num; ++i) {
            Matrix A = rand_m_norm_one(n, n, k);
            for (int i = 0; i < n; ++i) {
                A(i, i) += 1;
            }
            Matrix X_0 = rand_matrix(A.n, 1);
            mean += iter_num(zeidel, A, X_0);
        }
        out << k << ' ' << mean / double(att_num) << std::endl;
    }
}

void test_not_conv() {
    std::ofstream out("./not_conv.txt");
    Matrix A(2, 2, {1, 2, 3, 4});
    Matrix B(2, 1, {1, 2});
    Matrix true_X = rotation(A, B);
    Matrix X = rand_matrix(2, 1);
    int d = 3;
    for (int i = 0; i < 60; i += d) {
        X = zeidel(A, B, d, X);
        out << i << ' ' << dist(X, true_X) << std::endl;
    }
}

int main() {
    test_not_conv();
    test_iter_num();
    vector<int> iter_n;
    for (int i = 5; i < 100; i += 5) {
        iter_n.push_back(i);
    }
    test_convergention("./data_iter", iter_n, iter_method);
    test_convergention("./data_iter_2", iter_n, zeidel);
    return 0;
}
