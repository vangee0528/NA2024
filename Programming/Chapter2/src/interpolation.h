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

// Point 结构体
struct Point {
    double x;       // x坐标
    double y;       // y坐标
    double theta;   //与x轴正方向的夹角
    Point();
    Point(double x_val, double y_val);


    void print() const; 

    // 运算符重载
    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double scalar) const;
    Point operator/(double scalar) const;

};

// Bezier 类
class Bezier {
public:
    Bezier(const std::vector<Point>& control_points);

    void printOut() const;
    void FileOut(std::ofstream& outfile,double x) const; // 输出参数方程表达式到文件

private:
    std::vector<Point> control_points_;
};

// 辅助函数
int Factorial(int n);
int Combination(int n, int k);
std::string Bernstein(int n, int k); // 生成字符串形式的Bernstein多项式

// 心形函数
std::vector<double> heart_function(double x);   // 心形函数x对应的y值
Point heart_tangent(double x1,int i);           // 心形函数x对应的切向量

// 寻找心形函数上的 m+1 个点，并按顺序排列
std::vector<Point> find_points(int m);

// 将点列转换为控制点列
std::vector<std::vector<Point>> convert_to_control_points(const std::vector<Point>& points);


#endif // INTERPOLATION_H