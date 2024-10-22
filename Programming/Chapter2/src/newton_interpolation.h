#ifndef NEWTON_INTERPOLATION_H
#define NEWTON_INTERPOLATION_H

#include <vector>

// 差商的计算
double dividedDifference(const std::vector<double>& x, const std::vector<double>& f, int i, int j);

// 牛顿插值计算 p_n(f; x0, x1, ..., xn; x)
double newtonInterpolation(const std::vector<double>& x, const std::vector<double>& f, double x_val);


// 差商的计算
double dividedDifference(const std::vector<double>& x, const std::vector<double>& f, int i, int j) {
    if (i == j) {
        return f[i];
    }
    return (dividedDifference(x, f, i + 1, j) - dividedDifference(x, f, i, j - 1)) / (x[j] - x[i]);
}

// 牛顿插值计算 
double newtonInterpolation(const std::vector<double>& x, const std::vector<double>& f, double x_val) {
    int n = x.size();
    std::vector<double> coeffs(n);

    // 计算差商表的第一行（即牛顿插值公式中的系数）
    for (int i = 0; i < n; ++i) {
        coeffs[i] = dividedDifference(x, f, 0, i);
    }

    // 插值多项式的初始值
    double result = coeffs[0];
    double product = 1.0;

    // 计算插值多项式在 x_val 处的值
    for (int i = 1; i < n; ++i) {
        product *= (x_val - x[i - 1]);
        result += coeffs[i] * product;
    }

    return result;
}

#endif // NEWTON_INTERPOLATION_H
