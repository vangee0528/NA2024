#include "../Packages/spline.h"
#include <iostream>

// 定义MathFunction f(x) = 1/(1 + 25 * x^2)
double f(double x) {
    return 1.0 / (1 + 25 * x * x);
}

MathFunction f_func(f);

std::vector<int> Ns = {6, 11, 21, 41, 81} ;

int main () {
    for (int i = 0; i < Ns.size (); ++ i) {
        freopen (("output/problemA/N_" + std :: to_string (Ns [i]) + ".txt").c_str (), "w", stdout);
        PPSpline spline (1, 3, f_func, -1.0, 1.0, CLAMPED, Ns [i]); // Build a complete cubic spline
        spline.print ();
        double maxError = 0.0;
        for (int j = 0; j < Ns [i] - 1; ++ j) {
            
            double x = -1.0 + j * 2.0 / (Ns [i] - 1) + 1.0 / (Ns [i] - 1);
            double error = fabs (spline (x)[0] - f (x));
            if (error > maxError)
                maxError = error;
        }
        std :: cerr << "Error (N = " << Ns [i] << "): " << maxError << std :: endl;
        fclose (stdout);
    }
    return 0;
}