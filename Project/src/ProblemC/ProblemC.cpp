#include "../Packages/spline.h"
#include <iostream>

double f(double x) {
    return 1 / (1 + x * x);
}

MathFunction f_func(f);

int main() {
    try {
        std::vector<double> t1;
        for (int i = 1; i <= 11; i++) {
            t1.push_back(i - 6);
        }
        freopen("output/problemC/s23.txt", "w", stdout);
        BSpline spline1(1, 3, f_func, t1);
        spline1.print();
        fclose(stdout);

        std::vector<double> t2;
        for (int i = 1; i <= 10; i++) {
            t2.push_back(i - 5.5);
        }
        freopen("output/problemC/s12.txt", "w", stdout);
        BSpline spline2(1, 2, f_func, t2);
        spline2.print();
        fclose(stdout);
    } catch (const std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    } catch (const char *msg) {
        std::cerr << "Exception: " << msg << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception occurred" << std::endl;
    }

    return 0;
}