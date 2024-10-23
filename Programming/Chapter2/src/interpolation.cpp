#include "interpolation.h"
#include <iostream>
#include <fstream>

// POLYNOMIAL CLASS
Polynomial::Polynomial() {}

Polynomial::Polynomial(int degree, const std::vector<double>& coeffs) : coefficients(coeffs) {
    coefficients.resize(degree + 1, 0);
}

double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

Polynomial Polynomial::derivative() const {
    std::vector<double> newCoeffs(coefficients.size() - 1, 0.0);
    for (int i = 1; i < coefficients.size(); ++i) {
        newCoeffs[i - 1] = coefficients[i] * i;
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

void Polynomial::print() const {
    for (int i = coefficients.size() - 1; i >= 0; --i) {
        if (i < coefficients.size() - 1 && coefficients[i] >= 0) {
            std::cout << "+";
        }
        std::cout << std::fixed << std::setprecision(2) << coefficients[i] << "x^" << i << " ";
    }
    std::cout << std::endl;
}

Polynomial operator*(const Polynomial& p1, const Polynomial& p2) {
    std::vector<double> newCoeffs(p1.coefficients.size() + p2.coefficients.size() - 1, 0.0);
    for (int i = 0; i < p1.coefficients.size(); ++i) {
        for (int j = 0; j < p2.coefficients.size(); ++j) {
            newCoeffs[i + j] += p1.coefficients[i] * p2.coefficients[j];
        }
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

Polynomial operator+(const Polynomial& p1, const Polynomial& p2) {
    int maxSize = std::max(p1.coefficients.size(), p2.coefficients.size());
    std::vector<double> newCoeffs(maxSize, 0.0);

    for (int i = 0; i < maxSize; ++i) {
        if (i < p1.coefficients.size()) {
            newCoeffs[i] += p1.coefficients[i];
        }
        if (i < p2.coefficients.size()) {
            newCoeffs[i] += p2.coefficients[i];
        }
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

// NEWTON INTERPOLATOR CLASS
double NewtonInterpolator::dividedDifference(int i, int j) const {
    if (i == j) {
        return f_[i];
    }
    return (dividedDifference(i + 1, j) - dividedDifference(i, j - 1)) / (x_[j] - x_[i]);
}

NewtonInterpolator::NewtonInterpolator(const std::vector<double>& x, const std::vector<double>& f)
    : x_(x), f_(f), coeffs_(x.size()) {
    for (size_t i = 0; i < x.size(); ++i) {
        coeffs_[i] = dividedDifference(0, i);
    }
}

Polynomial NewtonInterpolator::interpolate() const {
    int n = x_.size();
    Polynomial result(0, {coeffs_[0]});
    Polynomial term(0, {1});

    for (int i = 1; i < n; ++i) {
        term = term * Polynomial(1, {-x_[i - 1], 1});
        result = result + term * Polynomial(0, {coeffs_[i]});
    }

    return result;
}

// HERMITE INTERPOLATOR CLASS
HermiteInterpolator::HermiteInterpolator(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& df)
    : x_(x), f_(f), df_(df), z_(2 * x.size()), q_(2 * x.size(), std::vector<double>(2 * x.size(), 0.0)) {
    computeDividedDifferences();
}

void HermiteInterpolator::computeDividedDifferences() {
    int n = x_.size();

    for (int i = 0; i < n; ++i) {
        z_[2 * i] = x_[i];
        z_[2 * i + 1] = x_[i];
        q_[2 * i][0] = f_[i];
        q_[2 * i + 1][0] = f_[i];
        q_[2 * i + 1][1] = df_[i];
        if (i != 0) {
            q_[2 * i][1] = (q_[2 * i][0] - q_[2 * i - 1][0]) / (z_[2 * i] - z_[2 * i - 1]);
        }
    }

    for (int i = 2; i < 2 * n; ++i) {
        for (int j = 2; j <= i; ++j) {
            q_[i][j] = (q_[i][j - 1] - q_[i - 1][j - 1]) / (z_[i] - z_[i - j]);
        }
    }
}

Polynomial HermiteInterpolator::interpolate() const {
    int n = x_.size();
    Polynomial result(0, {q_[0][0]});
    Polynomial term(0, {1});

    for (int i = 1; i < 2 * n; ++i) {
        term = term * Polynomial(1, {-z_[i - 1], 1});
        result = result + term * Polynomial(0, {q_[i][i]});
    }

    return result;
}


// Point 结构体实现
Point::Point() : x(0), y(0) {}

Point::Point(double x_val, double y_val) : x(x_val), y(y_val) {}

void Point::print() const {
    std::cout << "Point(" << x << ", " << y << ")" << std::endl;
}

Point Point::operator+(const Point& other) const {
    return Point(x + other.x, y + other.y);
}

Point Point::operator-(const Point& other) const {
    return Point(x - other.x, y - other.y);
}

Point Point::operator*(double scalar) const {
    return Point(x * scalar, y * scalar);
}

Point Point::operator/(double scalar) const {
    return Point(x / scalar, y / scalar);
}

// Bezier 类实现
Bezier::Bezier(const std::vector<Point>& control_points) : control_points_(control_points) {}
void Bezier::printOut() const {
    std::cout << "x(t) =";
    for (int i = 0; i < control_points_.size(); i++) {
        std::cout << control_points_[i].x << "*(" << Bernstein(control_points_.size() - 1, i) << ")";
        if (i != control_points_.size() - 1) {
            std::cout << "+";
        }
    }
    std::cout << std::endl;

    std::cout << "y(t) =";
    for (int i = 0; i < control_points_.size(); i++) {
        std::cout << control_points_[i].y << "*(" << Bernstein(control_points_.size() - 1, i) << ")";
        if (i != control_points_.size() - 1) {
            std::cout << "+";
        }
    }
    std::cout << std::endl;
}

void Bezier::FileOut(std::ofstream& outfile,double x) const {
    outfile << "x = " << x << std::endl;
    outfile << "x(t) =";
    for (int i = 0; i < control_points_.size(); i++) {
        outfile << control_points_[i].x << "*(" << Bernstein(control_points_.size() - 1, i) << ")";
        if (i != control_points_.size() - 1) {
            outfile << "+";
        }
    }
    outfile << std::endl;

    outfile << "y(t) =" ;
    for (int i = 0; i < control_points_.size(); i++) {
        outfile << control_points_[i].y << "*(" << Bernstein(control_points_.size() - 1, i) << ")";
        if (i != control_points_.size() - 1) {
            outfile << "+";
        }
    }
    outfile << std::endl;
}

// 辅助函数实现
int Factorial(int n) {
    int result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

int Combination(int n, int k) {
    return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

std::string Bernstein(int n, int k) {
    std::string result = "";
    int c = Combination(n, k);
    if (c != 1) {
        result += std::to_string(c) + "*";
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

// 心形函数实现
std::vector<double> heart_function(double x) {
    double y1 = 2.0 * (sqrt(3 - x * x) + sqrt(sqrt(x*x))) / 3.0;
    double y2 = 2.0 * (-sqrt(3 - x * x) + sqrt(sqrt(x*x))) / 3.0;
    std::vector<double> y = {y1, y2};
    return y;
}

// 计算心形函数的切线
Point heart_tangent(double x,int i) {
    double d1, d2;
    if(x >= 0){
        if(i==0) return Point(1.0, (-x/sqrt(3.0 - x * x) + 1.0/(2.0 * sqrt(x))) * 2.0 / 3.0);
        if(i==1) return Point(1.0, (x/sqrt(3.0 - x * x) + 1.0/(2.0 * sqrt(x))) * 2.0 / 3.0);
    }else{
        if(i==0) return Point(1.0, (-x/sqrt(3.0 - x * x) - 1.0/(2.0 * sqrt(-x))) * 2.0 / 3.0);
        if(i==1) return Point(1.0, (x/sqrt(3.0 - x * x) - 1.0/(2.0 * sqrt(-x))) * 2.0 / 3.0);
    }
    return Point(0,0);
}

// 寻找心形函数上的 m + 1 个点
std::vector<Point> find_points(int m) {
    std::vector<Point> points;
    int half_m = m / 2;
    for (int j = 0; j < half_m; j++) {
        double x1 = -1.7 + 1.6 * j / half_m;
        double x2 = 0.01 + 1.6 * j / half_m;
        double y1_1 = heart_function(x1)[0];
        double y1_2 = heart_function(x1)[1];
        double y2_1 = heart_function(x2)[0];
        double y2_2 = heart_function(x2)[1];
        points.push_back(Point(x1, y1_1));
        points.push_back(Point(x1, y1_2));
        points.push_back(Point(x2, y2_1));
        points.push_back(Point(x2, y2_2));
    }
    points.push_back(Point(-1.71, 1.09033));
    return points;
}


// 将点列转换为控制点
std::vector<std::vector<Point>> convert_to_control_points(const std::vector<Point>& points) {
    std::vector<std::vector<Point>> control_points;
    for (int i = 0; i < points.size() - 1; i++) {
        Point p_i = points[i];
        Point p_i_plus_1 = points[i + 1];
        Point q0 = p_i;
        Point q1 = p_i + heart_tangent(p_i.x,(i) % 2) / 3.0;
        Point q2 = p_i_plus_1 - heart_tangent(p_i_plus_1.x, (i) % 2) / 3.0;
        Point q3 = p_i_plus_1;
        control_points.push_back({q0, q1, q2, q3});
    }
    return control_points;
}