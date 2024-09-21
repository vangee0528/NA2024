#ifndef FUNCTION
#define FUNCTION
#include <cmath>
#include <limits>

class Function {
public:
    //用于计算函数值
    virtual double operator() (double x) const = 0;

    //用于计算导数（可以被重载）
    virtual double derivative(double x) const {
        //使用有限差分法
        double h = 1e-5;
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }

    // 析构函数
    virtual ~Function() {}

};

#endif
