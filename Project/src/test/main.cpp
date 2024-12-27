#include "../Packages/spline.h"
#include <iostream>

// 线性函数
double f1 (double x) {
    return 2*x+5;
}

// 二次函数
double f2 (double x) {
    return x*x+2*x+1;
}

// 高次函数
double f3 (double x) {
    return x*x*x+2*x*x+3*x+1;
}

// 指数函数，对数函数的组合
double f4 (double x) {
    return exp(x)*x-log(x);
}

// 三角函数的组合
double f5 (double x) {
    return sin(x)*cos(x)/x;
}


MathFunction f_func1 (f1);
MathFunction f_func2 (f2);
MathFunction f_func3 (f3);
MathFunction f_func4 (f4);
MathFunction f_func5 (f5);


int main () {
    freopen ("output/test/func.txt", "w", stdout);
    std::vector<MathFunction>  f_v = {f_func1};
    BSpline cubicSpline1 (1, 3, f_v, -1.0, 1.0, 21, NATURAL_SPLINE);
    cubicSpline1.print ();
    std::cout<<"===="<<std::endl;
    f_v = {f_func2};
    BSpline cubicSpline2 (1, 3, f_v, -1.0, 1.0, 21, NATURAL_SPLINE);
    cubicSpline2.print ();
    std::cout<<"===="<<std::endl;
    f_v = {f_func3};
    BSpline cubicSpline3 (1, 3, f_v, -1.0, 1.0, 21, NATURAL_SPLINE);
    cubicSpline3.print ();
    std::cout<<"===="<<std::endl;
    f_v = {f_func4};
    BSpline cubicSpline4 (1, 3, f_v, 0.1, 1.0, 21, NATURAL_SPLINE);
    cubicSpline4.print ();
    std::cout<<"===="<<std::endl;
    f_v = {f_func5};
    BSpline cubicSpline5 (1, 3, f_v, 0.1, 1.0, 21, NATURAL_SPLINE);
    cubicSpline5.print ();
    fclose (stdout);
    return 0;
}