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
int main() {
    std::vector<std::vector<double>> A =
    { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423},
      {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323},
      {0.0162, -0.1454, 0.9367, 0.0178, 0.0565},
      {0.0485, 0.0000, -0.1179, 0.9367, 0.0000},
      {0.0323, -0.0485, 0.2342, -0.0194, 0.6783} };

    std::vector<std::vector<double>> F = multMatrix(A, transpositionMatrix(A));
    std::vector<std::vector<double>> S = {
      {1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0} };

    for (int i = 3; i >= 0; --i) {
        std::vector<std::vector<double>> M = {
      {1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0} };
        std::vector<std::vector<double>> M1 = {
      {1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0} };
        for (int j = 0; j < 5; ++j) {
            M1[i][j] = F[i + 1][j];
            if (i != j) {
                M[i][j] = -F[i + 1][j] / F[i + 1][i];
            }
        }
        M[i][i] = 1 / F[i + 1][i];

        F = multMatrix(M1, multMatrix(F, M));
        S = multMatrix(S, M);
    }

    std::vector<double> lambda = { 0.23075, 0.31026, 0.691624, 0.936767, 1.17864 };
    system("color F0");
    for (int i = 0; i < 5; ++i) {
        std::vector<double> y = {0., 0., 0., 0., 0.};
        for (int j = 4; j >= 0; --j) {
            y[4 - j] = pow(lambda[i], j);
        }
        std::vector<double> s = multVec(S, y);
        double max = 0;
        for (int j = 0; j < 5; ++j) {
            if (abs(s[j]) > abs(max)) {
                max = s[j];
            }
        }
        std::cout << std::endl;
        std::cout << "=============================Eigenvector " << i + 1 << " =============================" << std::endl;
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << s[j] / max;
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "=============================Neural vector " << i + 1 << " =============================" << std::endl;
        std::vector<double> r = { 0., 0., 0., 0., 0. };
        for (int j = 0; j < 5; ++j) {
            r[j] = multVec(multMatrix(A, transpositionMatrix(A)), s)[j] - lambda[i] * s[j];
            std::cout << std::setw(15) << r[j];
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(15) << F[i][j];
        }
        std::cout << std::endl;
    }
}
