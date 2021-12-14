#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

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

double errorEstimation(std::vector<double> x2, std::vector<double> x1) {
    double cubeNorm = 0.;

    for (int i = 0; i < 5; ++i) {
        double res = abs(x2[i] - x1[i]);
        if (res > cubeNorm) {
            cubeNorm = res;
        }
    }

    return cubeNorm;
}

int main()
{
    const double EPS = 10E-5;
    const double w = 1.02;
    int count = 0;

    std::vector<std::vector<double>> A =
    { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423},
      {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323},
      {0.0162, -0.1454, 0.9367, 0.0178, 0.0565},
      {0.0485, 0.0000, -0.1179, 0.9367, 0.0000},
      {0.0323, -0.0485, 0.2342, -0.0194, 0.6783} };

    std::vector<double> f = { 1.5569, 2.0656, -2.9054, -8.0282, 3.4819 };

    //Для получения транспонированной положительно определённой матрицы системы
    std::vector<std::vector<double>> A_T = multMatrix(transpositionMatrix(A), A);
    std::vector<double> f_T = multVec(transpositionMatrix(A), f);

    std::vector<double> x0 = f;
    std::vector<double> xk = { 0, 0, 0, 0, 0 };

    while (errorEstimation(xk, x0) > EPS) {
        if (count != 0) {
            x0 = xk;
        }
        ++count;

        for (int i = 0; i < 5; ++i) {
            xk[i] = 0;
            xk[i] += (1 - w) * x0[i] + w * f[i] / A[i][i];

            double sum = 0.;
            for (int j = 0; j < i; ++j) {
                sum += A[i][j] * xk[j];
            }
            for (int j = i + 1; j < 5; ++j) {
                sum += A[i][j] * x0[j];
            }

            xk[i] -= w * sum / A[i][i];
        }
    }

    std::vector<double> neuralVector = multVec(A, xk);
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

    std::cout << "=================================================Resulting matrix A===============================================\n\n";
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << A_T[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "=================================================Resulting vector f===============================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << f_T[i] << ' ';
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "===================================================Decision vector==================================================\n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << xk[i] << ' ';
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
}
