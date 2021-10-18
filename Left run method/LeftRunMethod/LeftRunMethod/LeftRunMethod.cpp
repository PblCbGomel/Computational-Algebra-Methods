#include <iostream>
#include <iomanip>
#include <vector>

int main()
{
    std::vector<std::vector<double>> systemM =
    { {0.7941, 0., 0., 0., 0.},
      {-0.0485, 0.5168, 0.0, 0.0, 0.0},
      {0., -0.1454, 0.9367, 0.0178, 0.0},
      {0., 0., -0.1179, 0.9367, 0.0},
      {0., 0., 0., -0.0194, 0.6783} };

    std::vector<double> columnN = { 1.5569, 2.0656, -2.9054, -8.0282, 3.4819 };

    std::vector<double> answer = { 0, 0, 0, 0, 0 };
    std::vector<double> xi = { 0, 0, 0, 0, 0 };
    std::vector<double> nu = { 0, 0, 0, 0, 0 };

    //Вычисление прогоночных коэффициентов (прямой ход)
    xi[4] = -systemM[4][3] / systemM[4][4];
    nu[4] = columnN[4] / systemM[4][4];

    for (int i = 3; i > 0; --i) {
        xi[i] = -systemM[i][i - 1] / (systemM[i][i] + xi[i + 1] * systemM[i][i + 1]);
        nu[i] = (columnN[i] - systemM[i][i + 1] * nu[i + 1]) / (systemM[i][i] + xi[i + 1] * systemM[i][i + 1]);
    }

    nu[0] = (columnN[0] - systemM[0][1] * nu[1]) / (systemM[0][0] + xi[1] * systemM[0][1]);

    //Нахождение решения (обратный ход)
    answer[0] = nu[0];
    
    for (int i = 0; i < 4; ++i) {
        answer[i + 1] = xi[i + 1] * answer[i] + nu[i + 1];
    }

    //Посчитаем невязку
    std::vector<double> discrepancy = { 0, 0, 0, 0, 0 };
    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            discrepancy[i] += systemM[i][j] * answer[j];
        }
        discrepancy[i] -= columnN[i];
    }
    system("color F0");

    std::cout << "===================================================System matrix=================================================\n\n";

    for (int i = 0; i < systemM.size(); ++i) {
        for (int j = 0; j < systemM.size(); ++j) {
            std::cout << std::setw(15) << systemM[i][j] << ' ';
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

    std::cout << "                                        ( ";
    for (int i = 0; i < discrepancy.size(); ++i) {
        std::cout << discrepancy[i];
        if (i + 1 != discrepancy.size()) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl << std::endl;

    //Кубическая норма невязки
    double cubicNorm = 0;
    for (double item : discrepancy) {
        if (cubicNorm < abs(item)) {
            cubicNorm = abs(item);
        }
    }
    std::cout << "=====================================================Cubic norm======================================================\n\n";
    std::cout << "                                                     " << cubicNorm << std::endl << std::endl;
    
    std::cout << "=================================================Run coefficients xi==================================================\n\n";
    std::cout << "                                         ";
    for (int i = 1; i < 5; ++i) {
        std::cout << xi[i] << " ";
    }

    std::cout << std::endl << std::endl;
}
