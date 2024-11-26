#include "../Packages/spline.h"
#include <iostream>

double f (double x) {
    return 1/(1 + 25 * x * x);
}

MathFunction f_func (f);

int main () {
    freopen ("output/problemB/cubic.txt", "w", stdout);
    BSpline cubicSpline (1, 3, f, -1.0, 1.0, 21);
    cubicSpline.print ();
    fclose (stdout);
    freopen ("output/problemB/quadratic.txt", "w", stdout);
    BSpline quadSpline (1, 2, f, -1.0, 1.0, 21);
    quadSpline.print ();
    fclose (stdout);
    return 0;
}