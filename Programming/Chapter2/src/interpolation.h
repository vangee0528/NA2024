#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

class Polynomial {
public:
    std::vector<double> coefficients; 

    // 默认构造函数
    Polynomial();

    // 构造函数
    Polynomial(int degree, const std::vector<double>& coeffs);

    // 求值
    double evaluate(double x) const;

    // 求导
    Polynomial derivative() const;
    
    // 打印多项式
    void print() const;

    // 乘法和加法运算
    friend Polynomial operator*(const Polynomial& p1, const Polynomial& p2);
    friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
};

// 牛顿插值
class NewtonInterpolator  {
public:
    // 构造函数
    NewtonInterpolator(const std::vector<double>& x, const std::vector<double>& f);

    // 牛顿插值计算
    Polynomial interpolate() const;

private:
    std::vector<double> x_;    // 插值点
    std::vector<double> f_;    // 对应函数值
    std::vector<double> coeffs_; // 差商系数

    double dividedDifference(int i, int j) const;
};

class HermiteInterpolator {
public:
    // 构造函数
    HermiteInterpolator(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& df);

    // Hermite插值计算
    Polynomial interpolate() const;

private:
    std::vector<double> x_;    // 插值点
    std::vector<double> f_;    // 对应函数值
    std::vector<double> df_;   // 对应导数值
    std::vector<double> z_;    // Hermite节点
    std::vector<std::vector<double>> q_; // 差商表

    void computeDividedDifferences();
};


#endif // INTERPOLATION_H