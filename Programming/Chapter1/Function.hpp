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

};

#endif
