#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @class Polynomial
 * @brief 表示一个多项式。
 */
class Polynomial {
public:
    std::vector<double> coefficients; ///< 多项式的系数

    /** 
     * @brief 默认构造函数 
     */
    Polynomial();

    /** 
     * @brief 构造函数 
     * @param degree 多项式的次数
     * @param coeffs 多项式系数的向量
     */
    Polynomial(int degree, const std::vector<double>& coeffs);

    /** 
     * @brief 在给定点求多项式值 
     * @param x 自变量值
     * @return 在 x 处的多项式值
     */
    double evaluate(double x) const;

    /** 
     * @brief 计算多项式的导数 
     * @return 导数多项式
     */
    Polynomial derivative() const;
    
    /** 
     * @brief 打印多项式的表达式 
     */
    void print() const;

    /** 
     * @brief 将多项式输出到文件 
     * @param outfile 输出文件
     */
    void FileOut(std::ofstream& outfile) const;

    /** 
     * @brief 重载乘法运算符 
     * @param p1 第一个多项式
     * @param p2 第二个多项式
     * @return 乘积多项式
     */
    friend Polynomial operator*(const Polynomial& p1, const Polynomial& p2);

    /** 
     * @brief 重载加法运算符 
     * @param p1 第一个多项式
     * @param p2 第二个多项式
     * @return 和多项式
     */
    friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
};

/**
 * @class NewtonInterpolator
 * @brief 实现牛顿插值法。
 */
class NewtonInterpolator {
public:
    /** 
     * @brief 构造函数 
     * @param x 自变量值
     * @param f 因变量值
     */
    NewtonInterpolator(const std::vector<double>& x, const std::vector<double>& f);

    /** 
     * @brief 进行牛顿插值计算 
     * @return 插值多项式
     */
    Polynomial interpolate() const;

private:
    std::vector<double> x_;    ///< 插值点
    std::vector<double> f_;    ///< 对应函数值
    std::vector<double> coeffs_; ///< 差商系数

    /** 
     * @brief 计算分割差 
     * @param i 第一个索引
     * @param j 第二个索引
     * @return 分割差的值
     */
    double dividedDifference(int i, int j) const;
};

/**
 * @class HermiteInterpolator
 * @brief 实现埃尔米特插值法。
 */
class HermiteInterpolator {
public:
    /** 
     * @brief 构造函数 
     * @param x 自变量值
     * @param f 因变量值
     * @param df 导数值
     */
    HermiteInterpolator(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& df);

    /** 
     * @brief 进行 Hermite 插值计算 
     * @return 插值多项式
     */
    Polynomial interpolate() const;

private:
    std::vector<double> x_;    ///< 插值点
    std::vector<double> f_;    ///< 对应函数值
    std::vector<double> df_;   ///< 对应导数值
    std::vector<double> z_;    ///< Hermite节点
    std::vector<std::vector<double>> q_; ///< 差商表

    /** 
     * @brief 计算差商 
     */
    void computeDividedDifferences();
};

/**
 * @struct Point
 * @brief 表示二维空间中的一个点。
 */
struct Point {
    double x;       ///< x 坐标
    double y;       ///< y 坐标
    double theta;   ///< 与 x 轴正方向的夹角

    /** 
     * @brief 默认构造函数 
     */
    Point();
    
    /** 
     * @brief 构造函数 
     * @param x_val x 坐标值
     * @param y_val y 坐标值
     */
    Point(double x_val, double y_val);

    /** 
     * @brief 打印点的坐标 
     */
    void print() const; 

    /** 
     * @brief 重载加法运算符 
     * @param other 另一个点
     * @return 相加后的点
     */
    Point operator+(const Point& other) const;

    /** 
     * @brief 重载减法运算符 
     * @param other 另一个点
     * @return 相减后的点
     */
    Point operator-(const Point& other) const;

    /** 
     * @brief 重载乘法运算符 
     * @param scalar 缩放因子
     * @return 缩放后的点
     */
    Point operator*(double scalar) const;

    /** 
     * @brief 重载除法运算符 
     * @param scalar 缩放因子
     * @return 缩放后的点
     */
    Point operator/(double scalar) const;
};

/**
 * @class Bezier
 * @brief 表示由控制点定义的贝塞尔曲线。
 */
class Bezier {
public:
    /** 
     * @brief 构造函数 
     * @param control_points 控制点
     */
    Bezier(const std::vector<Point>& control_points);

    /** 
     * @brief 打印控制点信息 
     */
    void printOut() const;

    /** 
     * @brief 将参数方程输出到文件 
     * @param outfile 输出文件
     * @param x 自变量值
     */
    void FileOut(std::ofstream& outfile, double x) const; 

private:
    std::vector<Point> control_points_; ///< 控制点
};

/** 
 * @brief 计算阶乘 
 * @param n 输入值
 * @return n 的阶乘
 */
int Factorial(int n);

/** 
 * @brief 计算组合数 
 * @param n 总数
 * @param k 选择数
 * @return 组合数 C(n, k)
 */
int Combination(int n, int k);

/** 
 * @brief 生成字符串形式的 Bernstein 多项式 
 * @param n 多项式的次数
 * @param k 选项
 * @return Bernstein 多项式的字符串形式
 */
std::string Bernstein(int n, int k);

/** 
 * @brief 计算心形函数的 y 值 
 * @param x 自变量值
 * @return 对应的 y 值
 */
std::vector<double> heart_function(double x);

/** 
 * @brief 计算心形函数的切向量 
 * @param x1 自变量值
 * @param i 索引
 * @return 切向量
 */
Point heart_tangent(double x1, int i);

/** 
 * @brief 寻找心形函数上的 m+1 个点 
 * @param m 点的数量
 * @return 点的向量
 */
std::vector<Point> find_points(int m);

/** 
 * @brief 将点列转换为控制点列 
 * @param points 点的向量
 * @return 控制点的向量
 */
std::vector<std::vector<Point>> convert_to_control_points(const std::vector<Point>& points);

#endif // INTERPOLATION_H
