#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <cmath>
#include <iostream>

// 抽象基类 EquationSolver
class EquationSolver {
protected:
    const Function &F;  // 引用Function对象
public:
    EquationSolver(const Function &F) : F(F) {}
    virtual double solve() = 0;  // 纯虚函数
    virtual ~EquationSolver() {} // 虚析构函数
};

class Bisection_Method : public EquationSolver {
private:
    double a, b;
    double eps, delta;
    int Maxiter;
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-6, int Maxiter = 50) :
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
        int iter = 0;
        while (iter < Maxiter) {
            // 中点
            c = (a + b) / 2;
            fc = F(c);

            // 检查收敛条件：当前区间小于delta或函数值接近0
            if (std::fabs(fc) < eps || (b - a) / 2 < delta) {
                std::cout << "Converged after " << iter << " iterations." << std::endl;
                return c;  // 找到近似根
            }

            // 二分法步骤：选择下一步区间
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }

            iter++;
        }

        // 如果达到最大迭代次数，返回当前近似解并告知未达到精度
        throw std::runtime_error("Error: Maximum number of iterations reached without sufficient precision.");
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
        for (int i = 0; i < Maxiter; ++i) {
            fx = F(x);
            dfx = F.derivative(x);  // 假设 Function 提供导数函数
            if (std::fabs(fx) < eps) return x;  // 满足精度条件
            x = x - fx / dfx;  // 牛顿迭代公式
        }
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
        double eps = 1e-7, int Maxiter = 8) :
        EquationSolver(F), x0(x0), x1(x1), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() {
        double x2, fx0 = F(x0), fx1 = F(x1);
        for (int i = 0; i < Maxiter; ++i) {
            if (std::fabs(fx1 - fx0) < eps) return x1;  // 防止除零
            x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);  // 割线公式
            if (std::fabs(F(x2)) < eps) return x2;  // 满足收敛条件
            x0 = x1;
            fx0 = fx1;
            x1 = x2;
            fx1 = F(x1);
        }
        return x1;  // 返回最终结果
    }
};

#endif
