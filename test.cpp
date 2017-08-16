#include <fstream>
#include <iostream>
#include <utility>
#include <stdlib.h>
#include "matrix.cpp"
#include "solve_sys.cpp"
#include "solve_iter.cpp"

using std::vector;
using MethodT = Matrix(Matrix, Matrix);

Matrix rand_sq_matrix(int n) {
	return rand_matrix(n, n);
}

// возвращает точность и время
std::pair<double, double> do_test(int t_size, Matrix method(Matrix, Matrix), Matrix gen_m(int) = rand_sq_matrix) {
    Matrix A = gen_m(t_size);
    Matrix true_X = rand_matrix(t_size, 1);
    Matrix B = A * true_X;
    int t_0 = clock();
    Matrix X = method(A, B);
    Matrix dist(t_size, 1);
    for (int i = 0; i < t_size; ++i) {
        dist(i, 0) = X(i, 0) - true_X(i, 0);
    }
    return {log10(norm_one(dist) / norm_one(true_X)), double(clock() - t_0) / 1e6};
}

void get_time_data(std::string adr = "./data_time") {
	std::cout << "Time " << std::endl;
	vector<int> test_size;
	for (int i = 5; i < 50; i += 5) {
		test_size.push_back(i);
	}
	for (int i = 50; i < 100; i += 10) {
		test_size.push_back(i);
	}
	for (int i = 100; i < 200; i += 20){ 
		test_size.push_back(i);
	}
    int m_num = 4; // количество методов
    vector<vector<std::pair<double, double>>> res(m_num); // результаты
    int cnt = 0;
    for (int i : test_size) {
        std::cout << i << ' ' << cnt++ << "/" << test_size.size() << std::endl; 
        res[0].push_back(do_test(i, gauss));
        res[1].push_back(do_test(i, main_el));
        res[2].push_back(do_test(i, main_el_2));
        res[3].push_back(do_test(i, rotation));
    }
    // запись в файлы
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream time_f(adr + "_time.txt");
    for (int i : test_size) {
        dim_f << i << ' ';
    }
    dim_f << std::endl;
    for (int m_i = 0; m_i < m_num; ++m_i) {
        for (std::pair<double, double> p : res[m_i]) {
            time_f << p.second << ' ';
        }
        time_f << std::endl;
    }
}

void get_acc_data(std::string adr, int att_num = 1) {
	std::vector<int> test_size;
	for (int i = 5; i < 50; i += 5) {
		test_size.push_back(i);
	}
	for (int i = 50; i < 200; i += 20) {
		test_size.push_back(i);
	}
	for (int i = 200; i < 400; i += 40) {
		test_size.push_back(i);
	}
	std::cout << "Accuracy" << std::endl;
    int m_num = 4;
    vector<vector<double>> res(m_num);
    int cnt = 0;
    int64_t time_0 = clock();
    for (int i : test_size) {
        vector<double> mean_res(m_num, 0.0);
        for (int i_t = 0; i_t < att_num; ++i_t) {
            mean_res[0] += do_test(i, gauss).first;
            mean_res[1] += do_test(i, main_el).first;
            mean_res[2] += do_test(i, main_el_2).first;
            mean_res[3] += do_test(i, rotation).first;
        }
        for (int i_m = 0; i_m < m_num; ++i_m) {
            res[i_m].push_back(mean_res[i_m] / double(att_num));
        }
        std::cout << i << ' ' << ++cnt << "/" << test_size.size() << std::endl;
    }
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream acc_f(adr + "_acc.txt");
    for (int i : test_size) {
        dim_f << i << ' ';
    }
    dim_f << std::endl;
    for (int m_i = 0; m_i < m_num; ++m_i) {
        for (double mean_val : res[m_i]) {
            acc_f << mean_val << ' ';
        }
        acc_f << std::endl;
    }
    std::cout << double(clock() - time_0) / 1e6 << std::endl;
}

Matrix bad_matrix(int n) {
    double k = 1e-7;
    Matrix row = rand_matrix(1, n);
	Matrix A = rand_matrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = A(i, j) * k + row(0, j);
        }
    }
	return A;
}

Matrix hilb_m(int n) {
    Matrix ans(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ans(i, j) = 1.0 / double(i + j + 1);
        }
    }
    return ans;
}

double do_tests(int i, MethodT method, Matrix gen(int), int t_num = 30) {
    double ans = 0;
    for (int t_i = 0; t_i < t_num; ++t_i) {
        ans += do_test(i, method, gen).first;
    }
    return ans / t_num;
}

void test_acc_bad_m(std::string adr = "./data_bad_acc") {
	vector<int> test_size;
	for (int i = 5; i < 300; i += 20) {
		test_size.push_back(i);
	}
	std::cout << "Bad matrices" << std::endl;
    int m_num = 4;
    vector<vector<double>> res(m_num);
    int cnt = 0;
    int64_t time_0 = clock();
    for (int i : test_size) {
        res[0].push_back( do_tests(i, gauss, bad_matrix) );
        res[1].push_back( do_tests(i, main_el, bad_matrix) );
        res[2].push_back( do_tests(i, main_el_2, bad_matrix) );
        res[3].push_back( do_tests(i, rotation, bad_matrix) );
		std::cout << i << ' ' << ++cnt << "/" << test_size.size() << std::endl;
    }
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream acc_f(adr + "_acc.txt");
    for (int i : test_size) {
        dim_f << i << ' ';
    }
    dim_f << std::endl;
    for (int m_i = 0; m_i < m_num; ++m_i) {
        for (double mean_val : res[m_i]) {
            acc_f << mean_val << ' ';
        }
        acc_f << std::endl;
    }
    std::cout << double(clock() - time_0) / 1e6 << std::endl;
}

void test_hilbert_m() {
    std::ofstream out("./hilb.txt");
    vector<int> t_size;
    vector<double> res[4];
    for (int i = 1; i < 12; i++) {
        t_size.push_back(i);
        out << i << ' ';
    }
    out << std::endl;
    std::cout << "hilbert matrices" << std::endl;
    int cnt = 0;
    for (int i : t_size) {
        int a_n = 10;
        res[0].push_back( do_tests(i, gauss, hilb_m, a_n) );
        res[1].push_back( do_tests(i, main_el, hilb_m, a_n) );
        res[2].push_back( do_tests(i, main_el_2, hilb_m, a_n) );
        res[3].push_back( do_tests(i, rotation, hilb_m, a_n) );
		std::cout << i << ' ' << ++cnt << "/" << t_size.size() << std::endl;
    }
    for (int i = 0; i < 4; ++i) {
        for (double j : res[i]) {
            out << j << ' ';
        }
        out << std::endl;
    }
}

int main() {
    get_time_data(); // время
    get_acc_data("./data_1"); // точность при 1 измерении
    get_acc_data("./data_2", 30); // точность на хороших матрицах
    test_acc_bad_m(); // точность на матрицах, близких к линейно зависимым
    test_hilbert_m(); // матрицы Гильберта
    return 0;
}
