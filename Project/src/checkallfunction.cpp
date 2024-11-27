#include "Packages/spline.h"
#include <iostream>
#include <vector>
#include <functional>
#include <string>

double f(double x) {
    return 1.0 / (1 + 25 * x * x);
}

MathFunction f_func(f);

void check_1(std::string output_file) {
    std::cout << "[1] Now start to check Point1: pp-Form and Bspline of linear SPLINE S^0_1 (20pts)." << std::endl;
    std::cout << "We try to interpolate the function f(x) = 1/(1 + 25 * x^2) on the interval [-1, 1] with 20 points with both pp-Form and BSpline." << std::endl;
    int N = 20;
    PPSpline spline_pp(1, 1, f_func, -1.0, 1.0, CLAMPED, N);
    BSpline spline_bs(1, 1, f_func, -1.0, 1.0, N);
    freopen(output_file.c_str(), "w", stdout);
    std::cout << "CheckPoint 1" <<std::endl;
    spline_pp.print();
    std::cout << "===" <<std::endl;
    spline_bs.print();
    fclose(stdout);
}

void check_2() {
    // check_2 的实现
}

void check_3() {
    // check_3 的实现
}

int main() {
    std::vector<std::function<void(std::string)>> checks = {check_1};

    for (size_t i = 0; i < checks.size(); ++i) {
        std::string output_file = "output/check/check" + std::to_string(i + 1) + ".txt";
        
        checks[i](output_file);
        
        std::string command = "python3 src/plotcheck.py " + output_file;
        system(command.c_str());
    }

    return 0;
}