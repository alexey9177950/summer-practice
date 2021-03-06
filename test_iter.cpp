#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.cpp"
#include "solve_iter.cpp"
#include "solve_sys.cpp"

using IterM = Matrix(Matrix, Matrix, int, Matrix);
using std::vector;

// случайная матрица с заданной нормой L_1
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

// запись данных о скорости сходимости
void test_convergention(std::string addr, const vector<int>& iter_n, IterM method) {
	vector<int> dim{5, 10, 50, 200};
    std::ofstream file(addr);
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

void test_not_conv(std::string addr = "data_iter/not_conv.txt") {
    std::ofstream out(addr);
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
    double time_0 = clock();
    vector<int> iter_n;
    for (int i = 5; i < 100; i += 5) {
        iter_n.push_back(i);
    }
    test_convergention("data_iter/conv1.txt", iter_n, iter_method);
    test_convergention("data_iter/conv2.txt", iter_n, zeidel);
    test_not_conv();
    uint64_t time = (clock() - time_0) / CLOCKS_PER_SEC;
    std::cout << "TOTAL TIME: " << time / 60 << ":" << time % 60 << std::endl;
    return 0;
}
