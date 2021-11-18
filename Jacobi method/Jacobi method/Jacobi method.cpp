#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

bool diagDomination(const std::vector<std::vector<double>> M) {
    for (int i = 0; i < 5; ++i) {
        double sum_i = 0;
        for (int j = 0; j < 5; ++j) {
            if (i != j) {
                sum_i += abs(M[i][j]);
            }
        }
        if (abs(M[i][i]) <= sum_i) {
            return false;
        }
    }

    return true;
}

std::vector<std::vector<double>> transpositionMatrix(std::vector<std::vector<double>> M) {
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::swap(M[i][j], M[j][i]);
        }
    }

    return M;
}

std::vector<std::vector<double>> multMatrix(const std::vector<std::vector<double>> M1, const std::vector<std::vector<double>> M2) {
    std::vector<std::vector<double>> M3 = {
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0} };

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 5; ++k) {
                M3[i][j] += M1[i][k] * M2[k][j];
            }
        }
    }
    return M3;
}

std::vector<double> multVec(const std::vector<std::vector<double>> M, const std::vector<double> f) {
    std::vector<double> b = { 0., 0., 0., 0., 0. };

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            b[i] += M[i][j] * f[j];
        }
    }

    return b;
}

std::vector<std::vector<double>> getB(const std::vector<std::vector<double>> M) {
    std::vector<std::vector<double>> B = {
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0} };

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (i == j) {
                continue;
            }
            else {
                B[i][j] = -M[i][j] / M[i][i];
            }
        }
    }

    return B;
}

std::vector<double> getb(const std::vector<std::vector<double>> M, const std::vector<double> f) {
    std::vector<double> result = { 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int i = 0; i < 5; ++i) {
        result[i] = f[i] / M[i][i];
    }

    return result;
}

double cubeMatNorm(const std::vector<std::vector<double>> M) {
    double max = 0.;
    for (int i = 0; i < 5; ++i) {
        double strSum = 0.;
        for (int j = 0; j < 5; ++j) {
            strSum += abs(M[i][j]);
        }
        if (strSum > max) {
            max = strSum;
        }
    }
    return max;
}

double cubeNorm(std::vector<double> v1, std::vector<double> v2) {
    double max = 0.;

    for (int i = 0; i < 5; ++i) {
        if (abs(v1[i] - v2[i]) > max) {
            max = abs(v1[i] - v2[i]);
        }
    }

    return max;
}

double normVec(std::vector<double> b) {
    int max = 0;
    for (auto i : b) {
        if (max < abs(i)) {
            max = abs(i);
        }
    }
    return max;
}

int main() {
    const double EPS = 10E-5;
    int count = 0;

    std::vector<std::vector<double>> A =
    { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423},
      {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323},
      {0.0162, -0.1454, 0.9367, 0.0178, 0.0565},
      {0.0485, 0.0000, -0.1179, 0.9367, 0.0000},
      {0.0323, -0.0485, 0.2342, -0.0194, 0.6783} };

    std::vector<double> f = { 1.5569, 2.0656, -2.9054, -8.0282, 3.4819 };

    //Протранспонируем А, дабы потом умножить на неё же и получить симметрическую матрицу, не забывая умножить вектор f
   // std::vector<std::vector<double>> A_trans = transpositionMatrix(A);
   // std::vector<std::vector<double>> result_A = multMatrix(A_trans, A);
   // std::vector<double> result_f = multVec(A_trans, f);

    //Получим матрицу B
    std::vector<std::vector<double>> B = getB(A);
    std::vector<double> b = getb(A, f);
    
    //Проверим норму B для сходимости
    double normB = cubeMatNorm(B);

    //Присваиваем значения веткора b начальному приближению
    std::vector<double> xk0 = b;
    std::vector<double> xk1 = { 0, 0, 0, 0, 0 };

    //Пока наша кубическая норма больше 10E-5, проводим итерационный процесс
    while (cubeNorm(xk1, xk0) >= EPS) {
        if (count != 0) {
            xk0 = xk1;
        }
        count++;

        xk1 = multVec(B, xk0);
        for (int i = 0; i < 5; ++i) {
            xk1[i] += b[i];
        }
    }
    std::vector<double> neuralVector = multVec(A, xk1);
    for (int i = 0; i < 5; ++i) {
        neuralVector[i] -= f[i];
    }

    system("color F0");

    std::cout << "===================================================Default matrix=================================================\n\n";
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << A[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "==================================================Default vector f================================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << f[i] << ' ';
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << "=================================================Resulting matrix B===============================================\n\n";
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << B[i][j] << ' ';
        }
        std::cout << std::endl;
    } 
    std::cout << std::endl;

    std::cout << "===================================================Norm matrix B=================================================\n\n";
    std::cout << "                                                    " << normB;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "=================================================Resulting vector b===============================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << b[i] << ' ';
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "===========================================Check diagonal domination in A=========================================\n\n";
    if (diagDomination(A) == true) {
        std::cout << "                                                        true";
    }
    else {
        std::cout << "                                                        false";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "===================================================Decision vector==================================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << xk1[i] << ' ';
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "===================================================Iteration count==================================================\n\n";
    std::cout << "                                                          " << count;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "====================================================Neural vector==================================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << neuralVector[i] << ' ';
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "====================================================Number of iteration==================================================\n\n";
    std::cout << "\t\t\t\t\t\t\t" << (log(EPS*(1 - normB)/((1 - normB) * normVec(f) + normVec(b))))/(log(normB));
    std::cout << std::endl;
    std::cout << std::endl;
}
