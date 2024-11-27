#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <vector>
#include <iostream>
#include <cmath>
#include "lapack.h"

// 定义通用函数类
class MathFunction {
private:
    double (*function_ptr)(double x);
public:
    MathFunction() {}
    MathFunction(double (*func)(double x));
    virtual double evaluate(double x) const;
    virtual ~MathFunction() {}
};

// 定义多项式类
class Polynomial : public MathFunction {
private:
    std::vector<double> coefficients; // 存储多项式的系数（从低次到高次）
public:
    Polynomial() {}


    // 使用系数构造多项式
    Polynomial(const std::vector<double>& coef);

    // 使用 Newton 插值法构造多项式
    Polynomial(const std::vector<double>& x_values, const std::vector<double>& y_values);

    // 重载多项式的加法、减法、乘法
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;

    // 计算多项式在 x 点的值
    virtual double evaluate(double x) const;

    // 输出多项式公式
    void print() const;

    // 返回多项式的导数
    Polynomial derivative();

    virtual ~Polynomial() {}
};

// 分段多项式类
class PiecewisePolynomial : public MathFunction {
private:
    std::vector<double> points;       // 各分段的起点
    std::vector<Polynomial> polynomials;   // 各分段的多项式
public:
    PiecewisePolynomial() {}
    // 通过各段的表达式和分段点构造分段多项式
    PiecewisePolynomial(const std::vector<Polynomial>& p, const std::vector<double>& x): polynomials(p), points(x) {}

    // 计算分段多项式在 x 点的值
    virtual double evaluate(double x) const;

    // 输出分段多项式公式
    void print() const;

    // 返回分段多项式的阶数
    int segment_count() const { return polynomials.size(); }

    // 返回分段多项式的导数
    PiecewisePolynomial derivative();

    ~PiecewisePolynomial() {}
};

// 样条边界条件枚举
enum SplineBoundaryCondition {
    NO_CONDITION,              // 无边界条件
    CLAMPED,                   // 完全三次样条: s'(f;a) = f'(a) 且 s'(f;b) = f'(b)
    SECOND_DERIVATIVE_FIXED,   // 二阶导数固定: s''(f;a) = f''(a) 且 s''(f;b) = f''(b)
    NATURAL_SPLINE,            // 自然样条: s''(f;a) = 0 且 s''(f;b) = 0
    NOT_A_KNOT_CONDITION,      // 非节点条件: s'''(f;x) 在 x = x_2 和 x = x_{n-1} 处存在
    PERIODIC_CONDITION         // 周期边界条件: s(f;b) = s(f;a) 且 s'(f;b) = s'(f;a) 且 s''(f;b) = s''(f;a)
};

// 样条曲线基类
class Spline{
private:
    // 计算样条各段多项式
    virtual PiecewisePolynomial compute_spline_segments(SplineBoundaryCondition condition, 
                                                         const std::vector<double>& function_values, 
                                                         const std::vector<double>& nodes, 
                                                         double first_derivative_start, 
                                                         double first_derivative_end) = 0;
    
public:
    int dimensions;                             // 空间维度
    int spline_order;                           // 样条函数的次数即n, k = n - 1
    std::vector<PiecewisePolynomial> segments; // 每段的样条表达式
    Spline(int dim, int order) : dimensions(dim), spline_order(order) {}
    Spline(int dim, int order, const std::vector<PiecewisePolynomial>& spline_segments) 
        : dimensions(dim), spline_order(order), segments(spline_segments) {}
    std::vector<double> operator()(double t) const;
    void print() const;
    virtual ~Spline() {}
};

// 分段样条类
class PPSpline : public Spline {
private:
    PiecewisePolynomial compute_spline_segments(SplineBoundaryCondition bc, 
                                                 const std::vector<double>& f, 
                                                 const std::vector<double>& t, 
                                                 double da, 
                                                 double db);
public:
    PPSpline(int dim, int order, const MathFunction& f, 
                               double a, double d, 
                               SplineBoundaryCondition bc = NO_CONDITION, 
                               int num_interval = 100, 
                               double da = 0.0, 
                               double db = 0.0); // 均匀节点

    PPSpline(int dim, int order, const MathFunction& f, 
                               const std::vector<double>& t, 
                               SplineBoundaryCondition bc = NO_CONDITION, 
                               double da = 0.0, 
                               double db = 0.0); // 不均匀节点

    PPSpline(int dim, int order, const std::vector<std::vector<double>>& points, 
                               SplineBoundaryCondition bc = NO_CONDITION);
    virtual ~PPSpline() {}
};

// B 样条类
class BSpline : public Spline {
private:
    PiecewisePolynomial compute_spline_segments(SplineBoundaryCondition bc, 
                                                 const std::vector<double>& f,
                                                 const std::vector<double>& t, 
                                                 double da, 
                                                 double db);
    std::vector<double> knot_vector; // 节点向量

    // 计算 B 样条基函数 B_i^k 的值
    double evaluate_basis(int i, int k, double x) const;

    // 计算 B 样条基函数 B_i^k 的导数
    double evaluate_basis_derivative(int i, int k, double x) const;

    // 计算 B 样条基函数 B_i^k 的二阶导数
    double evaluate_basis_second_derivative(int i, int k, double x) const;

    // 计算 B 样条基函数 B_i^k 的三阶导数
    double evaluate_basis_third_derivative(int i, int k, double x) const;
public:
    BSpline(int dim, int order, const MathFunction& f, 
            double a, double b, int N = 100 ,SplineBoundaryCondition bc = NO_CONDITION); // 均匀节点

    BSpline(int dim, int order, const MathFunction& f, const std::vector<double>& time_points,SplineBoundaryCondition bc = NO_CONDITION); // 不均匀节点

    BSpline(int dim, int order, const std::vector<std::vector<double>>& points, SplineBoundaryCondition bc = NO_CONDITION); // 高维散点拟合
    virtual ~BSpline() {}
};

#endif // _SPLINE_H_