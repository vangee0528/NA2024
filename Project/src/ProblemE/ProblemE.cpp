#include "../Packages/spline.h"
#include <iostream>

// 心形曲线r1(t) = (sqrt(3) * cos(t), 2 / 3 * (sqrt(sqrt(3) * |cos(t)|) + sqrt(3) * sin(t)))
double r1x (double t) {
    return sqrt (3) * cos (t);
}

double r1y (double t) {
    return 2.0 / 3.0 * (sqrt (sqrt (3) * fabs (cos (t))) + sqrt (3) * sin (t));
}

MathFunction r1x_func (r1x);
MathFunction r1y_func (r1y);

int main () {
    std :: vector <int> n = {10, 40, 160};
    for (int i = 0; i < n.size (); ++ i) {
        std::vector<MathFunction> f_v = {r1x_func, r1y_func};
        PPSpline spline_unit (2, 3, f_v, -M_PI , M_PI, n [i], NATURAL_SPLINE);
        freopen (("output/problemE/r1_unit_N" + std :: to_string (n [i]) + ".txt").c_str (), "w", stdout);
        spline_unit.print ();
        fclose (stdout);

        PPSpline spline_chord (2, 3, f_v, -M_PI , M_PI, n [i], NATURAL_SPLINE, 0, 0, "chordal");
        freopen (("output/problemE/r1_chord_N" + std :: to_string (n [i]) + ".txt").c_str (), "w", stdout);
        spline_chord.print ();
        fclose (stdout);

    }
    return 0;
}