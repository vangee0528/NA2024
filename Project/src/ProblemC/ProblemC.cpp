#include "../Packages/spline.h"
#include <iostream>

double f(double x) {
    return 1 / (1 + x * x);
}

MathFunction f_func(f);

int main() {
    std::vector<MathFunction> f_v = {f_func};
    std::vector<double> t1;
    for (int i = 1; i <= 11; i++) {
        t1.push_back(i - 6);
    }
    freopen("output/problemC/s23.txt", "w", stdout);
    BSpline spline1(1, 3, f_v, t1, NATURAL_SPLINE);
    spline1.print();
    fclose(stdout);

    std::vector<double> t2;
    for (int i = 1; i <= 10; i++) {
        t2.push_back(i - 5.5);
    }
    freopen("output/problemC/s12.txt", "w", stdout);
    BSpline spline2(1, 2, f_v, t2);
    spline2.print();
    fclose(stdout);

    return 0;
}