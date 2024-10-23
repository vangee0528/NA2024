#include <vector>
#include <iostream>
#include <fstream>
#include "interpolation.h"
#include <cmath>

// 计算阶乘
int Factorial(int n)  {
     int result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

// 计算 C(n, k)
int Combination(int n, int k) {
    return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

//计算b_n,k(t),输出字符串
std::string Bernstein(int n, int k) {
    std::string result = "";
    int c = Combination(n, k);
    if (c != 1) {
        result += std::to_string(c)+"*";
    }
    if (n - k != 0) {
        result += "t**(" + std::to_string(n - k) + ")";
        if (k != 0) {
            result += "*";
        }
    }
    if (k != 0) {
        result += "(1-t)**(" + std::to_string(k) + ")";
    }

    return result;
}

struct Point {
    double x;  // x 坐标
    double y;  // y 坐标

    // 默认构造函数
    Point() : x(0), y(0) {}

    // 带参数的构造函数
    Point(double x_val, double y_val) : x(x_val), y(y_val) {}

    // 打印点的坐标
    void print() const {
        std::cout << "Point(" << x << ", " << y << ")" << std::endl;
    }

    // 两点相加（用于生成控制点等操作）
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y);
    }

    // 两点相减
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }

    // 缩放点（与常数相乘）
    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar);
    }

    // 缩放点（与常数相除）
    Point operator/(double scalar) const {
        return Point(x / scalar, y / scalar);
    }
};


// Bezier 曲线
class Bezier {
public:
    // 构造函数
    Bezier(const std::vector<Point>& control_points) : control_points_(control_points) {}

    // 计算 Bezier 曲线上的点
    void printOut() const{
        std::cout << "x(t) =";
        for (int i = 0; i < control_points_.size(); i++) {
             std::cout << control_points_[i].x <<"*("<< Bernstein(control_points_.size() - 1, i) << ")";
             if (i != control_points_.size() - 1) {
                 std::cout << "+";
             }
        }
        std::cout << std::endl;
        std::cout << "y(t) =";
        for (int i = 0; i < control_points_.size(); i++) {
             std::cout << control_points_[i].y <<"*("<< Bernstein(control_points_.size() - 1, i) << ")";
                if (i != control_points_.size() - 1) {
                    std::cout << "+";
                }
        }
        std::cout << ")" << std::endl;
    }

    void FileOut(std::ofstream& outfile) const{
        outfile << "x(t) =";
        for (int i = 0; i < control_points_.size(); i++) {
             outfile << control_points_[i].x <<"*("<< Bernstein(control_points_.size() - 1, i) << ")";
             if (i != control_points_.size() - 1) {
                 outfile << "+";
             }
        }
        outfile << std::endl;
        outfile << "y(t) =";
        for (int i = 0; i < control_points_.size(); i++) {
             outfile << control_points_[i].y <<"*("<< Bernstein(control_points_.size() - 1, i) << ")";
                if (i != control_points_.size() - 1) {
                    outfile << "+";
                }
        }
        outfile << ")" << std::endl;
    }


private:
    std::vector<Point> control_points_;  // 控制点

};

// 心形函数
std::vector<double> heart_function(double x) {
    double y1 = 2.0/3.0*(sqrt(3-x*x)+sqrt(abs(x)));
    double y2 = 2.0/3.0*(sqrt(3+x*x)-sqrt(abs(x)));
    std::vector<double> y = {y1, y2};
    return y;
}

// 心形函数固定点的切向量
Point heart_tangent(double x1) {
    double y1 = heart_function(x1)[0];
    double x2 = x1 + 0.0001;
    double y2 = heart_function(x2)[0];
    double slope = (y2 - y1) / (x2 - x1);
    double angle = atan(slope);
    Point tangent(cos(angle), sin(angle));
    return tangent;
}


//在心形函数上寻找m个点
std::vector<Point> find_points(int m) {
    std::vector<Point> points;
    for(int j = 0 ; j < m / 2 ; j++){
        double x = -1.7 + 3.4 * j / (m / 2);
        double y1 = heart_function(x)[0];
        double y2 = heart_function(x)[1];
        points.push_back(Point(x, y1));
        points.push_back(Point(x, y2));
    }
    return points;
}

// 将一个点列转化为四个控制点的点列
std::vector<std::vector<Point>> convert_to_control_points(const std::vector<Point>& points) {
    std::vector<std::vector<Point>> control_points;
    for (int i = 0; i < points.size()-1; i++) {
        Point p_i = points[i];
        Point p_i_plus_1 = points[(i + 1)];
        Point q0, q1, q2, q3;
        q0 = p_i;
        q1 = p_i + heart_tangent(p_i.x) / 3.0;
        q2 = p_i_plus_1 - heart_tangent(p_i_plus_1.x) / 3.0;
        q3 = p_i_plus_1;
        control_points.push_back({q0, q1, q2, q3});
    }
    return control_points;
}




