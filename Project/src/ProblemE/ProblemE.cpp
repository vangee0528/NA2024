#include "../Packages/spline.h"
#include <iostream>

// 心形曲线r1(t) = (sqrt(3) * cos(t), 2 / 3 * (sqrt(sqrt(3) * |cos(t)|) + sqrt(3) * sin(t)))
double r1x (double t) {
    return sqrt (3) * cos (t);
}

double r1y (double t) {
    return 2.0 / 3.0 * (sqrt (sqrt (3) * fabs (cos (t))) + sqrt (3) * sin (t));
}

int main () {
    std :: vector <int> n = {10, 40, 160};
    for (int i = 0; i < n.size (); ++ i) {
        PPSpline splineX (1, 3, r1x, -M_PI / 2, M_PI / 2, NATURAL_SPLINE, n [i]);
        PPSpline splineY (1, 3, r1y, -M_PI / 2, M_PI / 2, NATURAL_SPLINE, n [i]);
        freopen (("output/problemE/r1x_" + std :: to_string (n [i]) + ".txt").c_str (), "w", stdout);
        splineX.print ();
        fclose (stdout);
        freopen (("output/problemE/r1y_" + std :: to_string (n [i]) + ".txt").c_str (), "w", stdout);
        splineY.print ();
        fclose (stdout);
    }
    return 0;
}