#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Function.hpp"          // 引入 Function 虚类
#include "EquationSolver.hpp"    // 引入 EquationSolver 虚类

class Function2 : public Function {
public:
    virtual double operator() (double x) const override {
        return x - tan(x); // f(x) = x - tan(x)
    }

    virtual double derivative(double x) const override {
        return 1 - 1 / cos(x) / cos(x); // f'(x) = 1 - sec^2(x)
    }
};

void checkRoot(const Function &f, double root) {
    if(fabs(f(root)) < 1e-7) {
        std::cout << "\033[32mRoot is correct.\033[0m" << std::endl << std::endl;
        std::cout << "--------------------------------" << std::endl << std::endl;
    } else {
        std::cout << "\033[31mRoot is incorrect.\033[0m" << std::endl << std::endl;
        std::cout << "--------------------------------" << std::endl << std::endl;
    }
}


int main() {
    try {
        Function2 f2;
        std::cout << "--------------------------------" << std::endl << std::endl;

        // 寻找接近 4.5 的根
        std::cout << "Testing f(x) = x - tan(x) near 4.5:" << std::endl;
        double x1 = 4.5;
        Newton_Method solver1(f2, x1);
        double root1 = solver1.solve();
        std::cout << "Root near 4.5: " << root1 << std::endl;
        std::cout << "f(root) = " << f2(root1) << std::endl;
        checkRoot(f2, root1);

        // 寻找接近 7.7 的根
        std::cout << "Testing f(x) = x - tan(x) near 7.7:" << std::endl;
        double x2 = 7.7;
        Newton_Method solver2(f2, x2);
        double root2 = solver2.solve();
        std::cout << "Root near 7.7: " << root2 << std::endl;
        std::cout << "f(root) = " << f2(root2) << std::endl;
        checkRoot(f2, root2);

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
