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

double getLambda(const std::vector<double> f1, const std::vector<double> f2) {
    double res = 0.;

    for (int i = 0; i < 5; ++i) {
        res += f1[i]/ f2[i];
    }
    
    return res / 5.;
}

std::vector<double> normV(std::vector<double> f) {
    double max = 0;
    for (int i = 0; i < 5; ++i) {
        if (abs(f[i]) > abs(max)) {
            max = f[i];
        }
    }
    for (int i = 0; i < 5; ++i) {
        f[i] /= max;
    }
    return f;
}

double getEps(double lk, double lk1) {
    return abs(lk - lk1);
}

int main()
{
    const double EPS = 0.00001;

    std::vector<std::vector<double>> A1 = {
        {1.67953, 0.291252, 0.513653,  0.32638, -0.737916},
        {0.291252, 3.87842, 0.500924, 0.384012, -0.280008 },
        {0.513653, 0.500924, 1.5147, 0.051476, -0.903941},
        {-0.32638, 0.384012, 0.051476, 1.24609, 0.169331},
        {-0.737916, -0.280008, -0.903941, 0.169331, 2.59987} };

    std::vector<double> y = {1, 1, 1, 1, 1};

    double l = 1;
    double l1 = 0;
    int counter = 0;

    while (getEps(l, l1) > EPS) {
        counter++;
        std::vector<double> y1 = multVec(A1, y);
        l = l1;
        l1 = getLambda(y1, y);
        y = normV(y1);

    }

    y = { 0.315612, 1, 0.395623, 0.070598, -0.4947832 };

    l1 = 1 / l1;
    system("color F0");
    std::cout << "=======================================Min lambda=======================================" << std::endl;
    std::cout << "                                            " << l1 << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "=====================================Eigenvector=======================================" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setw(15) << y[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "=====================================Number of iterations=======================================" << std::endl;
    std::cout << "                                               " << counter;
    std::cout << std::endl;
    std::cout << std::endl;

}

