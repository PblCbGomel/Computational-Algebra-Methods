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

    std::vector<std::vector<double>> A_Result = multMatrix(A, transpositionMatrix(A));
    
    std::vector<std::vector<double>> C;
    C.push_back({1., 0., 0., 0., 0.});

    for (int i = 0; i < 5; ++i) {
        std::vector<double> Ck = multVec(A_Result, C[i]);
        C.push_back(Ck);
    }

    std::vector<std::vector<double>> C_result = { 
       {C[4][0], C[3][0], C[2][0], C[1][0], C[0][0], C[5][0]},
      {C[4][1], C[3][1], C[2][1], C[1][1], C[0][1], C[5][1]},
      {C[4][2], C[3][2], C[2][2], C[1][2], C[0][2], C[5][2]},
      {C[4][3], C[3][3], C[2][3], C[1][3], C[0][3], C[5][3]},
      {C[4][4], C[3][4], C[2][4], C[1][4], C[0][4], C[5][4]} };

    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            for (int k = i + 1; k < 6; ++k) {
                C_result[j][k] -= C_result[j][i] / C_result[i][i] * C_result[i][k];
            }
        }
        for (int j = i + 1; j < 5; ++j) {
            C_result[j][i] = 0.;
        }
    }

    std::vector<double> answer = { 0., 0., 0., 0., 0. };

    for (int k = 4; k >= 0; --k) {
        double sum = 0;
        for (int j = k + 1; j < 5; ++j) {
            sum += C_result[k][j] * answer[j];
        }
        answer[k] = (C_result[k][5] - sum) / C_result[k][k];
    }
    system("color F0");
    std::cout << "====================================================qk====================================================";
    std::cout << std::endl;
    std::cout << std::endl << "              ";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(12) << answer[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::vector<double> lambda = {0.23075, 0.31026, 0.691624, 0.936767, 1.17864};
    std::cout << "=============================================Lambda (Wolfram)=============================================";
    std::cout << std::endl;
    std::cout << std::endl << "              ";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(12) << lambda[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "=============================================Neural=============================================";
    std::cout << std::endl;
    std::cout << std::endl << "              ";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(12) << pow(lambda[i], 5) -answer[0] * pow(lambda[i], 4) - answer[1] * pow(lambda[i], 3) - 
            answer[2] * pow(lambda[i], 2) - answer[3] * pow(lambda[i], 1) - answer[4] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}