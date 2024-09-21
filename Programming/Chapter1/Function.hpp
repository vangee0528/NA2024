#ifndef FUNCTION
#define FUNCTION
#include <cmath>
#include <limits>

class Function {
public:
    // 纯虚函数，用于计算函数值
    virtual double operator() (double x) const = 0;

    // 虚函数，用于计算导数（可以被重载）
    virtual double derivative(double x) const {
        // 数值导数近似法，使用有限差分法
        double h = 1e-5;  // 取一个非常小的值
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }


    // 虚析构函数
    virtual ~Function() {}

    // 判断函数是否在区间 [a, b] 上连续
    virtual bool isContinuous(double a, double b, int samples = 100) const {
        double step = (b - a) / samples;
        for (int i = 0; i < samples; ++i) {
            double x1 = a + i * step;
            double x2 = a + (i + 1) * step;
            if (std::fabs((*this)(x2) - (*this)(x1)) > std::numeric_limits<double>::epsilon()) {
                return false; // 存在跳跃或不连续
            }
        }
        return true; // 连续
    }
};

#endif
