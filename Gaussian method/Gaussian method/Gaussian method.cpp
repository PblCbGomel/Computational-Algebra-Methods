#include <iostream>
#include <iomanip>
#include <vector>

int searchMaxI(std::vector<std::vector<double>> M, int startI, int startJ) {
    int maxI = startI;
    int maxJ = startJ;
    for (int i = startI; i < M.size(); ++i) {
        for (int j = startJ; j < M.size(); ++j) {
            if (abs(M[i][j]) > abs(M[maxI][maxJ])) {
                maxI = i;
                maxJ = j;
            }
        }
    }

    return maxI;
}

int searchMaxJ(std::vector<std::vector<double>> M, int startI, int startJ) {
    int maxI = startI;
    int maxJ = startJ;
    for (int i = startI; i < M.size(); ++i) {
        for (int j = startJ; j < M.size(); ++j) {
            if (abs(M[i][j]) > abs(M[maxI][maxJ])) {
                maxI = i;
                maxJ = j;
            }
        }
    }

    return maxJ;
}

int main() {
    //Исходные данные
    std::vector<std::vector<double>> systemM = 
      { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423},
        {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323},
        {0.0162, -0.1454, 0.9367, 0.0178, 0.0565},
        {0.0485, 0.0000, -0.1179, 0.9367, 0.0000},
        {0.0323, -0.0485, 0.2342, -0.0194, 0.6783} };

    std::vector<double> columnN = { 1.5569, 2.0656, -2.9054, -8.0282, 3.4819 };

    std::vector<std::vector<double>> extendedSystemMat =
    { {0.7941, 0.0000, -0.2067, 0.1454, 0.2423, 1.5569},
      {-0.0485, 0.5168, 0.0000, -0.0985, 0.0323, 2.0656},
      {0.0162, -0.1454, 0.9367, 0.0178, 0.0565, -2.9054},
      {0.0485, 0.0000, -0.1179, 0.9367, 0.0000, -8.0282},
      {0.0323, -0.0485, 0.2342, -0.0194, 0.6783, 3.4819} };

    //Для обратной матрицы
    std::vector<std::vector<double>> inverseMatrix =
    { {1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1} };

    std::vector<int> permutations = { 0, 1, 2, 3, 4 };
    std::vector<int> inverseRowPermutations;
    std::vector<int> inverseColPermutations = { 0, 1, 2, 3, 4 };
    std::vector<double> answer = { 0, 0, 0, 0, 0 };

    double det = 1;
    
    //Прямой ход
    for (int k = 0; k < systemM.size(); ++k) {
        int maxI = searchMaxI(extendedSystemMat, k, k);
        int maxJ = searchMaxJ(extendedSystemMat, k, k);

        std::swap(permutations[k], permutations[maxJ]);
        std::swap(inverseColPermutations[k], inverseColPermutations[maxI]);

        //Найдя максимальный, свапаем строку, потом столбец
        for (int j = k; j < systemM.size() + 1; ++j) {
            std::swap(extendedSystemMat[k][j], extendedSystemMat[maxI][j]);
        }       
        for (int i = 0; i < systemM.size(); ++i) {
            std::swap(extendedSystemMat[i][k], extendedSystemMat[i][maxJ]);
        }

        //Так же меняем строки и столцы для подсчёта обратной матрицы
        for (int j = 0; j < systemM.size(); ++j) {
            std::swap(inverseMatrix[k][j], inverseMatrix[maxI][j]);
        }
        for (int i = 0; i < systemM.size(); ++i) {
            std::swap(inverseMatrix[i][k], inverseMatrix[i][maxJ]);
        }
       
        //Параллельно считаем определитель
        det *= extendedSystemMat[k][k] * pow(-1, maxI + maxJ - 2 * k);

        //Для обратной матрицы 
        for (int i = k + 1; i < systemM.size(); ++i) {
            for (int j = 0; j < systemM.size(); ++j) {
                inverseMatrix[i][j] -= inverseMatrix[k][j] * extendedSystemMat[i][k] / extendedSystemMat[k][k];
            }
        }

        //Далее k-й шаг метода, не трогая k-й столбец
        for (int i = k + 1; i < systemM.size(); ++i) {
            for (int j = k + 1;  j < systemM.size() + 1; ++j) {
                extendedSystemMat[i][j] -= extendedSystemMat[i][k] * extendedSystemMat[k][j] / extendedSystemMat[k][k];
            }
        }

        //Зануляем k-й столбец за исключением X k-й
        for (int i = k + 1; i < systemM.size(); ++i) {
            extendedSystemMat[i][k] = 0;
        }

    } //Прямой ход метода Гаусса завершен, переход к обратному ходу
    inverseRowPermutations = permutations; 

    for (int k = systemM.size() - 1; k >= 0; --k) {
        double sum = 0;
        for (int j = k + 1; j < systemM.size(); ++j) {
            sum += extendedSystemMat[k][j] * answer[j];
        }
        answer[k] = (extendedSystemMat[k][systemM.size()] - sum) / extendedSystemMat[k][k];
    }
    
    //Расставим переменные по своим местам
    int num = 0;
    for (int i = 0; i < answer.size(); ++i) {
        if (permutations[i] == num) {
            std::swap(answer[num], answer[i]);
            std::swap(permutations[num], permutations[i]);
            i = num;
            ++num;
        }
    }//Обратный ход метода Гаусса завершён

    //Посчитаем невязку
    std::vector<double> discrepancy = { 0, 0, 0, 0, 0 };
    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            discrepancy[i] += systemM[i][j] * answer[j];
        }
        discrepancy[i] -= columnN[i];
    }

    //Досчитаем обратную матрицу
    for (int k = systemM.size() - 1; k >= 0; --k) {
        for (int j = 0; j < systemM.size(); ++j) {
            inverseMatrix[k][j] /= extendedSystemMat[k][k];
        }
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < systemM.size(); ++j) {
                inverseMatrix[i][j] -= inverseMatrix[k][j] * extendedSystemMat[i][k];
            }
        }
    }

    //Вернём строки и столбцы обратной матрицы на своё место
    num = 0;
    for (int k = 0; k < inverseRowPermutations.size(); ++k) {
        if (inverseRowPermutations[k] == num) {
            for (int j = 0; j < systemM.size(); ++j) {
                std::swap(inverseMatrix[num][j], inverseMatrix[k][j]);
            }
            std::swap(inverseRowPermutations[num], inverseRowPermutations[k]);
            k = num;
            ++num;
        }
    }
    num = 0;
    for (int k = 0; k < inverseColPermutations.size(); ++k) {
        if (inverseColPermutations[k] == num) {
            for (int i = 0; i < systemM.size(); ++i) {
                std::swap(inverseMatrix[i][num], inverseMatrix[i][k]);
            }
            std::swap(inverseColPermutations[num], inverseColPermutations[k]);
            k = num;
            ++num;
        }
    }

    //Невязка для обратной матрицы
    std::vector<std::vector<double>> neuralMatrix =
    { {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0} };
    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            double sum = 0;
            for (int k = 0; k < systemM.size(); ++k) {
                sum += systemM[i][k] * inverseMatrix[k][j];
            }
            neuralMatrix[i][j] = sum;
            if (i == j) {
                neuralMatrix[i][j] -= 1;
            }
        }
    }

    //Найдём кубические нормы матриц и число обусловленности
    double normSystemM = 0;
    double maxSum = 0;
    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            maxSum += abs(systemM[i][j]);
        }
        if (maxSum > normSystemM) {
            normSystemM = maxSum;
        }
        maxSum = 0;
    }

    double normInverseM = 0;
    for (int i = 0; i < inverseMatrix.size(); ++i) {
        for (int j = 0; j < inverseMatrix.size(); ++j) {
           maxSum += abs(inverseMatrix[i][j]);
        }
        if (maxSum > normInverseM) {
            normInverseM = maxSum;
        }
        maxSum = 0;
    }

    double normNeuralM = 0;
    for (int i = 0; i < neuralMatrix.size(); ++i) {
        for (int j = 0; j < neuralMatrix.size(); ++j) {
            maxSum += abs(neuralMatrix[i][j]);
        }
        if (maxSum > normNeuralM) {
            normNeuralM = maxSum;
        }
        maxSum = 0;
    }

    double conditionNumber = normInverseM * normSystemM;

    std::cout << "===================================================Resulting matrix=================================================\n\n";

    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size() + 1; ++j) {
            if (j == 5) {
                std::cout << " | ";
            }
            std::cout << std::setw(15) << extendedSystemMat[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "===================================================Decision vector==================================================\n\n";

    std::cout << "                                ( ";
    for (int i = 0; i < answer.size(); ++i) {
        std::cout << answer[i];
        if (i + 1 != answer.size()) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl << std::endl;

    std::cout << "===================================================Neural vector====================================================\n\n";

    std::cout << "                                            ( ";
    for (int i = 0; i < discrepancy.size(); ++i) {
        std::cout << discrepancy[i];
        if (i + 1 != discrepancy.size()) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl << std::endl;

    std::cout << "===================================================Determinant======================================================\n\n";
    std::cout << "                                                    " << det << std::endl << std::endl;
    
    std::cout << "===================================================Inverce Matrix======================================================\n\n";

    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            if (j == 5) {
                std::cout << " | ";
            }
            std::cout << std::setw(15) << inverseMatrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "===================================================Neural Matrix======================================================\n\n";

    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            std::cout << std::setw(15) << neuralMatrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "===================================================Condition number=====================================================\n\n";
    std::cout << "                                                       " << conditionNumber << std::endl << std::endl;

    return 0;
}
