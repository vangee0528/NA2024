#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "newton_interpolation.h"  // 牛顿插值法头文件

// 目标函数 f(x) = 1 / (1 + 25x^2)
double f(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

// 生成 Chebyshev 插值点
std::vector<double> chebyshev_points(int n) {
    std::vector<double> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = cos(M_PI * (2 * i + 1) / (2 * n));
    }
    return points;
}

int main() {
    // 定义 n 值
    std::vector<int> n_values = {5, 10, 15, 20};

    // 插值点区间
    double a = -1.0, b = 1.0;

    // 生成不同 n 值的插值数据并写入文件
    std::ofstream outfile("data/interpolation_chebyshev_data.txt");

    // 生成数据
    for (int n : n_values) {
        std::vector<double> x_vals = chebyshev_points(n), f_vals(n);

        // 计算 Chebyshev 插值点的函数值
        for (int i = 0; i < n; ++i) {
            f_vals[i] = f(x_vals[i]);
        }

        // 输出 n 值
        outfile << "n = " << n << std::endl;

        // 对 x 范围内的点进行插值计算
        for (double x = a; x <= b; x += 0.05) {
            double interpolated_value = newtonInterpolation(x_vals, f_vals, x);
            outfile << x << " " << interpolated_value << std::endl;
        }

        // 插入分隔符以便Python解析
        outfile << "====" << std::endl;
    }

    outfile.close();

    return 0;
}
