#include "Function.hpp"
#include "EquationSolver.hpp"
#include <cmath>
#include <iostream>

// 定义 F 类继承自 Function，用于表示方程
class F : public Function {
public:
    // 构造函数，初始化参数 l, h, D, beta_1
    F(double l, double h, double D, double beta_1) : l(l), h(h), D(D), beta_1(beta_1) {}

    // 重载 operator()，用于计算方程 f(alpha)
    double operator()(double alpha) const override {
        double sin_alpha = std::sin(alpha);
        double cos_alpha = std::cos(alpha);
        double tan_beta_1 = std::tan(beta_1);
        double sin_beta_1 = std::sin(beta_1);
        double cos_beta_1 = std::cos(beta_1);

        double A = l * sin_beta_1;
        double B = l * cos_beta_1;
        double C = (h + 0.5 * D) * sin_beta_1 - 0.5 * D * tan_beta_1;
        double E = (h + 0.5 * D) * cos_beta_1 - 0.5 * D;

        return A * sin_alpha * cos_alpha + B * sin_alpha * sin_alpha - C * cos_alpha - E * sin_alpha;
    }

private:
    double l, h, D, beta_1; // 定义参数 l, h, D, beta_1
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

// 检查解是否符合精度要求

int main() {
    // 题目要求的常量
    double l = 89.0 * 0.0254;   // 转换为米
    double h = 49.0 * 0.0254;   // 转换为米
    double D = 55.0 * 0.0254;   // 转换为米
    double beta_1 = 11.5 * M_PI / 180.0;  // 转换为弧度

    // 构造函数对象
    F f(l, h, D, beta_1);
    
    std::cout << "-----------PROBLEM F---------------" << std::endl << std::endl;

    // (a) 二分法求解大概值
    std::cout << "(a) Use the bisection method to find approximate values of alpha when l = 89 in., h = 49 in., D = 55 in., and beta_1 = 11.5 degrees:\n";
    Bisection_Method solver1(f, 0, M_PI, 1e-3, 1, 1000); //将epsilon和delta设置的更宽泛
    double root1 = solver1.solve();
    std::cout << "\033[35mApproximate root:" << root1 * 180 / M_PI << " degrees" << "\033[0m" << std::endl << std::endl; // 弧度转换角度


    // (b) 牛顿法求解，初值 33°
    std::cout << "--------------------------------" << std::endl << std::endl;
    std::cout << "(b) Use Newton's method to find alpha with the initial guess 33 degrees when l,h,D, and beta_1 are as given above:\n";
    double initial_guess = 33.0 * M_PI / 180.0; 
    Newton_Method solver2(f, initial_guess);
    double root2 = solver2.solve();
    std::cout << "\033[33mRoot near 33 degrees: " << root2 * 180 / M_PI << "\033[0m" << std::endl;
    std::cout << "f(root) = " << f(root2) << std::endl;
    checkRoot(f, root2);

    // (c) 割线法，选择远离33°的初值，寻找其他根
    std::cout << "Use secant method with initial guesses far away from 33 to find alpha:\n" << std::endl;;
    Secant_Method solver3(f, 100, 200);
    double root3 = solver3.solve();
    std::cout << "Initial guesses: 100, 200" << std::endl;
    std::cout << "\033[34mRoot: " << root3 * 180 / M_PI << "\033[0m" << std::endl;
    std::cout << "f(root) = " << f(root3) << std::endl;
    checkRoot(f, root3);

    Secant_Method solver5(f, -1, -2);
    double root4 = solver5.solve();
    std::cout << "Initial guesses: -1, -2" << std::endl;
    std::cout << "\033[34mRoot: " << root4 * 180 / M_PI << "\033[0m" << std::endl;
    std::cout << "f(root) = " << f(root4) << std::endl;
    checkRoot(f, root4);

    Secant_Method solver4(f, 10, 20);
    double root5 = solver4.solve();
    std::cout << "Initial guesses: 10, 20" << std::endl;
    std::cout << "\033[34mRoot: " << root5 * 180 / M_PI << "\033[0m" << std::endl;
    std::cout << "f(root) = " << f(root5) << std::endl;
    checkRoot(f, root5);



    return 0;
}
