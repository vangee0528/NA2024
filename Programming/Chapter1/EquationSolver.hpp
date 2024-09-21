#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <cmath>
#include <iostream>

// 方程求解器基类
class EquationSolver {
protected:
    const Function &F; 
public:
    EquationSolver(const Function &F) : F(F) {}
    virtual double solve() = 0;  
    virtual ~EquationSolver() {} 
};

// 二分法求解器类
class Bisection_Method : public EquationSolver {
private:
    double a, b; // 区间端点
    double eps, delta;// 收敛精度和区间精度
    int Maxiter; // 最大迭代次数
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-7, int Maxiter = 1000) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}

    virtual double solve() {
        double fa = F(a);
        double fb = F(b);

        // 预处理：检查区间两端符号是否相反
        if (fa * fb > 0) {
            throw std::runtime_error("Error: The function values at the endpoints of the interval must have opposite signs.");
        }


        // 二分迭代
        double c, fc;
        for (int i = 0; i < Maxiter; ++i) {
            // 中点
            c = (a + b) / 2;
            fc = F(c);

            // // 当前区间小于delta或函数值接近0
            // if (std::fabs(fc) < eps || (b - a) / 2 < delta) {
            //     std::cout << "Converged after " << iter << " iterations." << std::endl;
            //     return c;  // 找到近似根
            // }

            // 当前区间小于delta且函数值接近0
            if (std::fabs(fc) < eps && (b - a) / 2 < delta) { //为了获得更好的近似根，同时满足区间精度和收敛精度
                std::cout << "Converged after " << i << " iterations." << std::endl;
                return c;  // 找到近似根
            }

            //选择下一步区间
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }

        }
        std::cout << "\033[1;31mError: Maximum iterations reached without convergence.\033[0m" << std::endl;
        return NAN;

    }
};




// 牛顿法求解器类
class Newton_Method : public EquationSolver {
private:
    double x0;  // 初始猜测值
    double eps;  // 收敛精度
    int Maxiter;  // 最大迭代次数
public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-7, int Maxiter = 8) :
        EquationSolver(F), x0(x0), Maxiter(Maxiter), eps(eps) {}
    
    virtual double solve() {
        double x = x0, fx, dfx;
        for (int i = 0; i < Maxiter+1; ++i) {
            fx = F(x);
            dfx = F.derivative(x);  
            if (std::fabs(fx) < eps) {
                std::cout << "Converged after " << i << " iterations." << std::endl;
                return x;  // 如果已经满足精度要求，返回根
            }
            if (dfx == 0) throw std::runtime_error("Error: Derivative is zero."); // 防止除零
            x = x - fx / dfx;  
        }

        // 超过最大迭代次数
        std::cout << "\033[1;31mError: Maximum iterations reached without convergence.\033[0m" << std::endl;
        return NAN;
        
        return x;  // 返回最终结果
    }
};


// 割线法求解器类
class Secant_Method : public EquationSolver {
private:
    double x0, x1;  // 初始猜测值
    double eps;  // 收敛精度
    int Maxiter;  // 最大迭代次数
public:
    Secant_Method(const Function &F, double x0, double x1, 
        double eps = 1e-7, int Maxiter = 100) :
        EquationSolver(F), x0(x0), x1(x1), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() {
        double x2, fx0 = F(x0), fx1 = F(x1);
        for (int i = 0; i < Maxiter+1; ++i) {
            if (std::fabs(fx1 - fx0) < eps) {
                std::cout << "\033[1;31mError: Function values are too close.\033[0m" << std::endl;
                return NAN;
            }
            x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);  
            if (std::fabs(F(x2)) < eps) {
                std::cout << "Converged after " << i << " iterations." << std::endl;
                return x2;  // 找到近似根
            }
            x0 = x1;
            fx0 = fx1;
            x1 = x2;
            fx1 = F(x1);
        }
        // 超过最大迭代次数
        std::cout << "\033[1;31mError: Maximum iterations reached without convergence.\033[0m" << std::endl;
        return NAN;
    }
};


#endif
