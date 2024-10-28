#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include "interpolation.h"  

// 目标函数 f(x) = 1 / (1 + x^2)
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {

    // 参数设置
    std::vector<int> n_values = {2, 4, 6, 8};
    std::cout << "============Problem B=============" << std::endl;
    double a = -5.0, b = 5.0;
    std::ofstream outfile("data/ProblemB.txt");

    // 生成数据
    for (int n : n_values) {
        std::vector<double> x_vals(n + 1), f_vals(n + 1);

        for (int i = 0; i <= n; ++i) {
            x_vals[i] = a + (b - a) * i / n;
            f_vals[i] = f(x_vals[i]);
        }

        outfile << "n = " << n << std::endl;

        NewtonInterpolator interpolator(x_vals, f_vals);

        std::cout << "n = " << n << " : ";
        Polynomial interpolating_polynomial = interpolator.interpolate();

        interpolating_polynomial.print();

        for (double x = a; x <= b; x += 0.1) {
            double interpolated_value = interpolating_polynomial.evaluate(x); 
            outfile << x << " " << interpolated_value << std::endl;
        }

        outfile << "====" << std::endl;
    }

    outfile.close();

    return 0;
}