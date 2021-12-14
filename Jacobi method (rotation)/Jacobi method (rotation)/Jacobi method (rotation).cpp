#include <iostream>
#include <vector>
#include <iomanip>

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

double getEps(std::vector<std::vector<double>> M) {
    double eps = 0.;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (i == j) {
                continue;
            }
            eps += pow(M[i][j], 2);
        }
    }

    return sqrt(eps);
}

int main() {
    std::vector<std::vector<double>> A =
    { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423},
      {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323},
      {0.0162, -0.1454, 0.9367, 0.0178, 0.0565},
      {0.0485, 0.0000, -0.1179, 0.9367, 0.0000},
      {0.0323, -0.0485, 0.2342, -0.0194, 0.6783} };

    std::vector<std::vector<double>> M = multMatrix(A, transpositionMatrix(A));

    const double EPS = 0.00001;

    int counter = 0;

    std::vector<std::vector<double>> resultT = {
      {1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0} };

    while (getEps(M) > EPS) {
        counter++;
        double max = 0;
        int maxI = 0;
        int maxJ = 0;
        
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == j) {
                    continue;
                }
                else {
                    if (abs(M[i][j]) > abs(max)) {
                        max = M[i][j];
                        maxI = i;
                        maxJ = j;
                    }
                }
            }
        }

        double tg = 2 * max / (M[maxI][maxI] - M[maxJ][maxJ]);
        double cos = sqrt((1 + 1 / sqrt(1 + pow(tg, 2))) / 2);
        double sin = sqrt(1 - pow(cos, 2));
        
        std::vector<std::vector<double>> T = {
      {1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0} };
        T[maxI][maxI] = cos;
        T[maxJ][maxJ] = cos;
        T[maxI][maxJ] = -sin;
        T[maxJ][maxI] = sin;

        M = multMatrix(transpositionMatrix(T), multMatrix(M, T));
        resultT = multMatrix(resultT, T);
    }
    system("color F0");
    std::cout << "================================Result Matrix (diag - lambda)================================" << std::endl;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(12) << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "=============================================Neural=============================================";
    std::cout << std::endl;
    std::cout << std::endl << "              ";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(12) << pow(M[i][i], 5) - 3.34804 * pow(M[i][i], 4) + 4.1574 * pow(M[i][i], 3) -
            2.35346 * pow(M[i][i], 2) + 0.596922 * pow(M[i][i], 1) - 0.0546702 << " ";
    }  
    std::cout << std::endl;
    std::cout << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << "=============================Eigenvector " << i + 1 << " =============================" << std::endl;
        double max = 0;
        for (int j = 0; j < 5; ++j) {
            if (abs(resultT[j][i]) > abs(max)) {
                max = resultT[j][i];
            }
        }
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << resultT[j][i] / max;
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "=============================Neural vector " << i + 1 << " =============================" << std::endl;
        std::vector<double> r = { 0., 0., 0., 0., 0. };
        std::vector<double> s = { resultT[0][i], resultT[1][i], resultT[2][i], resultT[3][i], resultT[4][i] };
        for (int j = 0; j < 5; ++j) {
            r[j] = multVec(multMatrix(A, transpositionMatrix(A)), s)[j] - M[i][i] * s[j];
            std::cout << std::setw(15) << r[j];
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
    std::cout << "=====================================Number of iterations=======================================\n";
    std::cout << "                               " << counter;
    std::cout << std::endl;
    std::cout << std::endl;
}