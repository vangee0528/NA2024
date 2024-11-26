#include "../spline.h"
#include <iostream>

// 定义MathFunction f(x) = 1/(1 + 25 * x^2)
double f(double x) {
    return 1.0 / (1 + 25 * x * x);
}

MathFunction f_func(f);

int main() {
    std::cout << f_func.evaluate(1.0) <<std::endl;
}