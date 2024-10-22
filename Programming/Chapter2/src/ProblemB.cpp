#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include "newton_interpolation.h"  // 引入牛顿插值法的头文件

// 定义目标函数 f(x) = 1 / (1 + x^2)
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {
    // 定义 n 值
    std::vector<int> n_values = {2, 4, 6, 8};

    // 插值点区间
    double a = -5.0, b = 5.0;

    // 生成不同 n 值的插值数据并写入文件
    std::ofstream outfile("data/interpolation_data.txt");

    // 生成数据
    for (int n : n_values) {
        std::vector<double> x_vals(n + 1), f_vals(n + 1);

        // 计算插值点和对应的函数值
        for (int i = 0; i <= n; ++i) {
            x_vals[i] = a + (b - a) * i / n;
            f_vals[i] = f(x_vals[i]);
        }

        // 输出 n 值
        outfile << "n = " << n << std::endl;

        // 对 x 范围内的点进行插值计算
        for (double x = a; x <= b; x += 0.1) {
            double interpolated_value = newtonInterpolation(x_vals, f_vals, x);
            outfile << x << " " << interpolated_value << std::endl;
        }

        // 插入分隔符以便Python解析
        outfile << "====" << std::endl;
    }

    outfile.close();

    return 0;
}
