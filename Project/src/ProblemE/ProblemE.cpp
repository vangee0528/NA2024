#include "../Packages/spline.h"
#include <iostream>
#include <vector>
#include <cmath>

std::vector<std::vector<MathFunction>> fvs;

// 心形曲线r1(t) = (sqrt(3) * cos(t), 2 / 3 * (sqrt(sqrt(3) * |cos(t)|) + sqrt(3) * sin(t)))
double r1x(double t) {
    return sqrt(3) * cos(t);
}

double r1y(double t) {
    return 2.0 / 3.0 * (sqrt(sqrt(3) * fabs(cos(t))) + sqrt(3) * sin(t));
}

MathFunction r1x_func(r1x);
MathFunction r1y_func(r1y);


// 曲线r2(t) = (sin(t) + tcos(t), cos(t) - tsin(t))
double r2x(double t) {
    return sin(t) + t * cos(t);
}

double r2y(double t) {
    return cos(t) - t * sin(t);
}

MathFunction r2x_func(r2x);
MathFunction r2y_func(r2y);



// 曲线r3(t) = [sin(cos(t))*cos(sin(t)), sin(cos(t))*sin(sin(t)), cos(cos(t))]
double r3x(double t) {
    return sin(cos(t)) * cos(sin(t));
}

double r3y(double t) {
    return sin(cos(t)) * sin(sin(t));
}

double r3z(double t) {
    return cos(cos(t));
}

MathFunction r3x_func(r3x);
MathFunction r3y_func(r3y);
MathFunction r3z_func(r3z);

int main() {
    std::vector<int> n = {10, 40, 160};
        std::vector<std::vector<double>> segs = {{-M_PI, M_PI}, {0, 6 * M_PI}, {0, 2 * M_PI}};

        fvs.push_back({r1x_func, r1y_func});
        fvs.push_back({r2x_func, r2y_func});
        fvs.push_back({r3x_func, r3y_func, r3z_func});

        for (int j = 0; j < fvs.size(); ++j) {
            std::vector<MathFunction> f_v = fvs[j];
            for (int i = 0; i < n.size(); ++i) {
                PPSpline spline_unit(f_v.size(), 3, f_v, segs[j][0], segs[j][1], n[i], NATURAL_SPLINE);
                freopen(("output/problemE/r" + std::to_string(j+1) + "_unit_N" + std::to_string(n[i]) + ".txt").c_str(), "w", stdout);
                spline_unit.print();
                fclose(stdout);

                PPSpline spline_chord(f_v.size(), 3, f_v, segs[j][0], segs[j][1], n[i], NATURAL_SPLINE, 0, 0, "chordal");
                freopen(("output/problemE/r" + std::to_string(j+1) + "_chord_N" + std::to_string(n[i]) + ".txt").c_str(), "w", stdout);
                spline_chord.print();
                fclose(stdout);
            }
        }

        // 使用一个B样条来验证在相同的边界条件和节点时，PPSpline 和 BSpline 是否会得到相同的结果
        BSpline spline_compared(fvs[0].size(), 3, fvs[0], segs[0][0], segs[0][1], n[0], NATURAL_SPLINE);
        freopen("output/problemE/r1_unit_N10_Bspline.txt", "w", stdout);
        spline_compared.print();
        fclose(stdout);

    return 0;
}