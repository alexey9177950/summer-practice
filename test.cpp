#include <fstream>
#include <functional>
#include <iostream>
#include <utility>
#include <stdlib.h>
#include "matrix.cpp"
#include "solve_sys.cpp"
#include "solve_iter.cpp"

using std::vector;
using MethodT = std::function<Matrix(Matrix, Matrix)>;

// реализованные методы и их количетство
int m_num = 4;
vector<MethodT> methods{gauss, main_el, main_el_2, rotation};

// возвращает квадратную матрицу из равномерного распредения на [0, 1]
Matrix rand_sq_matrix(int n) {
	return rand_matrix(n, n);
}

// возвращает точность и время
std::pair<double, double> test(int t_size, MethodT method, Matrix gen_m(int) = rand_sq_matrix) {
    Matrix A = gen_m(t_size);
    Matrix true_X = rand_matrix(t_size, 1);
    Matrix B = A * true_X;
    int t_0 = clock();
    Matrix X = method(A, B);
    Matrix dist(t_size, 1);
    for (int i = 0; i < t_size; ++i) {
        dist(i, 0) = X(i, 0) - true_X(i, 0);
    }
    return {log10(norm_one(dist) / norm_one(true_X)), double(clock() - t_0) / CLOCKS_PER_SEC};
}

// получает и записывает данные о времени работы
void get_time_data(std::string adr = "data/time") {
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
        for (int m_i = 0; m_i < m_num; ++m_i) {
            res[m_i].push_back(test(i, methods[m_i]));
        }
    }
    // запись в файлы
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream time_f(adr + ".txt");
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

// возвращает среднюю точность после нескольких измерений
double tests(int i, MethodT method, Matrix gen(int), int t_num = 30) {
    double ans = 0;
    for (int t_i = 0; t_i < t_num; ++t_i) {
        ans += test(i, method, gen).first;
    }
    return ans / t_num;
}

// получает и записывает данные о точности
void get_acc_data(std::string adr = "data/acc", int att_num = 1) {
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
	std::cout << "Accuracy, " << att_num << std::endl;
    int m_num = 4;
    vector<vector<double>> res(m_num);
    int cnt = 0;
    int64_t time_0 = clock();
    for (int i : test_size) {
        for (int m_i = 0; m_i < m_num; ++m_i) {
            res[m_i].push_back( tests(i, methods[m_i], rand_sq_matrix, att_num) );
        }
        std::cout << i << ' ' << ++cnt << "/" << test_size.size() << std::endl;
    }
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream acc_f(adr + ".txt");
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
    std::cout << double(clock() - time_0) / CLOCKS_PER_SEC << std::endl;
}

// генерирует квадратную матрицу, близкую к линейно зависимой
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

// возвращает матрицу Гильберта размера n на n
Matrix hilb_m(int n) {
    Matrix ans(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ans(i, j) = 1.0 / double(i + j + 1);
        }
    }
    return ans;
}

// получает и записывает данные о точности на плохих матрицах (близких к линейно зависимым)
void test_acc_bad_m(std::string adr = "data/bad_acc") {
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
        for (int m_i = 0; m_i < m_num; ++m_i) {
            res[m_i].push_back( tests(i, methods[m_i], bad_matrix) );
		}
        std::cout << i << ' ' << ++cnt << "/" << test_size.size() << std::endl;
    }
    std::ofstream dim_f(adr + "_dim.txt");
    std::ofstream acc_f(adr + ".txt");
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
    std::cout << double(clock() - time_0) / CLOCKS_PER_SEC << std::endl;
}

// получает и записывает данные о точности на матрицах Гильберта
void test_hilbert_m(std::string addr = "data/hilb") {
    std::ofstream out(addr + ".txt");
    vector<int> t_size;
    vector<double> res[4];
    for (int i = 1; i < 12; i++) {
        t_size.push_back(i);
        out << i << ' ';
    }
    out << std::endl;
    std::cout << "hilbert matrices" << std::endl;
    int cnt = 0;
    int a_n = 10;
    for (int i : t_size) {
        for (int m_i = 0; m_i < m_num; ++m_i) {
            res[m_i].push_back( tests(i, methods[m_i], hilb_m, a_n) );
		}
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
    uint64_t time_0 = clock();
    get_time_data(); // время
    get_acc_data("data/acc_1"); // точность при 1 измерении
    get_acc_data("data/acc_2", 30); // точность на хороших матрицах
    test_acc_bad_m(); // точность на матрицах, близких к линейно зависимым
    test_hilbert_m(); // матрицы Гильберта
    uint64_t time = (clock() - time_0) / CLOCKS_PER_SEC;
    std::cout << "TOTAL TIME: " << time / 60 << ":" << time % 60 << std::endl;
    return 0;
}
