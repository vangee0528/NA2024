#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "interpolation.h"  

// 目标函数 f(x) = 1 / (1 + 25x^2)
double f(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}


std::vector<double> chebyshev_points(int n) {
    std::vector<double> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = cos(M_PI * (2 * i + 1) / (2 * n));
    }
    return points;
}

int main() {

    // 参数设置
    std::vector<int> n_values = {5, 10, 15, 20};
    double a = -1.0, b = 1.0;
    std::ofstream outfile("data/ProblemC.txt");

    // 生成数据
    for (int n : n_values) {
        std::vector<double> x_vals = chebyshev_points(n), f_vals(n);

        // 计算 Chebyshev 插值点的函数值
        for (int i = 0; i < n; ++i) {
            f_vals[i] = f(x_vals[i]);
        }

        outfile << "n = " << n << std::endl;

        NewtonInterpolator interpolator(x_vals, f_vals);

        Polynomial interpolating_polynomial = interpolator.interpolate();

        for (double x = a; x <= b; x += 0.05) {
            double interpolated_value = interpolating_polynomial.evaluate(x);  
            outfile << x << " " << interpolated_value << std::endl;
        }

        outfile << "====" << std::endl;
    }

    outfile.close();

    return 0;
}