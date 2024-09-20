#ifndef FUNCTION
#define FUNCTION

class Function {
public:
    virtual double operator() (double x) const = 0;
    virtual double derivative(double x) const {
        /* Type your code here. */
        return 0;   // to be replaced
    }
};

#endif