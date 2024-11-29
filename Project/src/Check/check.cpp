#include "../Packages/spline.h"
#include <iostream>
#include <cstdlib>
#include <algorithm>

// Point 1 : PP Spline & BSpline S^0_1 (20pts)

// 测试函数 f(x) = e^x - x^2 + 1

double f1(double x) {
    return exp(x) - x * x + 1;
}

MathFunction f1_func(f1);
std::vector<MathFunction> f_v = {f1_func};

void check_P1() {
    
    PPSpline ppspline(1, 1, f_v, -1, 1, 40);
    BSpline bspline(1, 1, f_v, -1, 1, 40);
    
    FILE* original_stdout = stdout;
    
    freopen("output/check/P1_ppspline.txt", "w", stdout);
    ppspline.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出

    freopen("output/check/P1_bspline.txt", "w", stdout);
    bspline.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出
    system("python3 src/Check/plot.py output/check/P1_ppspline.txt");
    system("python3 src/Check/plot.py output/check/P1_bspline.txt");
    stdout = original_stdout;
}

// Point 2 : PP Spline S^2_3 with 3 different boundary conditions of any knots (75pts)
void check_P2() {

    // 在 [-1, 1] 上不均匀的随机选取 11 个节点
    std::vector<double> t1;
    for (int i = 1; i <= 11; i++) {
        t1.push_back(-1 + 2.0 * rand() / RAND_MAX);
    }
    std::sort(t1.begin(), t1.end()); // 使用 std::sort 进行排序
    
    FILE* original_stdout = stdout;

    freopen("output/check/P2_s23_natural.txt", "w", stdout);
    PPSpline spline1(1, 3, f_v, t1, NATURAL_SPLINE);
    spline1.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出

    freopen("output/check/P2_s23_clamped.txt", "w", stdout);
    PPSpline spline2(1, 3, f_v, t1, CLAMPED);
    spline2.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出

    freopen("output/check/P2_s23_periodic.txt", "w", stdout);
    PPSpline spline3(1, 3, f_v, t1, PERIODIC_CONDITION);
    spline3.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出
    system("python3 src/Check/plot.py output/check/P2_s23_natural.txt");
    system("python3 src/Check/plot.py output/check/P2_s23_clamped.txt");
    system("python3 src/Check/plot.py output/check/P2_s23_periodic.txt");
    stdout = original_stdout;
}

// Point 3 : B Spline S^2_3 with 3 different boundary conditions of any knots (75pts)

void check_P3() {
    // 在 [-1, 1] 上不均匀的随机选取 11 个节点
    std::vector<double> t1;
    for (int i = 1; i <= 11; i++) {
        t1.push_back(-1 + 2.0 * rand() / RAND_MAX);
    }
    std::sort(t1.begin(), t1.end()); // 使用 std::sort 进行排序
    
    FILE* original_stdout = stdout;

    freopen("output/check/P3_s23_natural.txt", "w", stdout);
    BSpline spline1(1, 3, {f_v}, t1, NATURAL_SPLINE);
    spline1.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出

    freopen("output/check/P3_s23_clamped.txt", "w", stdout);
    BSpline spline2(1, 3, f_v, t1, CLAMPED);
    spline2.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出

    freopen("output/check/P3_s23_periodic.txt", "w", stdout);
    BSpline spline3(1, 3, f_v, t1, PERIODIC_CONDITION);
    spline3.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出
    system("python3 src/Check/plot.py output/check/P3_s23_natural.txt");
    system("python3 src/Check/plot.py output/check/P3_s23_clamped.txt");
    system("python3 src/Check/plot.py output/check/P3_s23_periodic.txt");
    stdout = original_stdout;

}

// Point 4 : Check PP Spline & BSpline get the same curve with the same boundary conditions and knots 

// 通过绘制output/ProbelmE/r1_unit_N10.txt 和 output/ProbelmE/r1_unit_N10_Bspline.txt 进行比较
void check_P4(){
    system("python3 src/Check/compare.py output/problemE/r1_unit_N10.txt output/problemE/r1_unit_N10_Bspline.txt");
}


// Point 5 : BSpline in any order 
void check_P5(){
    //随机生成节点序列N=11
    std::vector<double> t1;
    for (int i = 1; i <= 11; i++) {
        t1.push_back(-0.5 + 2.0 * rand() / RAND_MAX);
    }
    std::sort(t1.begin(), t1.end()); // 使用 std::sort 进行排序

    std::vector<double> coefficients;
    for (int i = 1; i <= 14; i++) {
        coefficients.push_back(-1.0 + 2.0 * rand() / RAND_MAX);
    }
    BSpline spline1(1, 4, coefficients, t1);
    FILE* original_stdout = stdout;
    freopen("output/check/P5_bspline.txt", "w", stdout);
    spline1.print();
    fflush(stdout);
    freopen("/dev/tty", "a", stdout); // 恢复标准输出
    system("python3 src/Check/plot.py output/check/P5_bspline.txt");
    stdout = original_stdout;
}

int main() {
    std::cout << "Checking P1 : PP Spline & BSpline S^0_1" << std::endl;
    check_P1();
    std::cout << "figures are saved in output/check" << std::endl;
    std::cout << "Checking P2 : S^2_3 SPLINE OF PP-Form with 3 different boundary conditions" << std::endl;
    check_P2();
    std::cout << "figures are saved in output/check" << std::endl;
    std::cout << "Checking P3 : S^2_3 SPLINE OF B-Form with 3 different boundary conditions" << std::endl;
    check_P3();
    std::cout << "figures are saved in output/check" << std::endl;
    std::cout << "Checking P4 : Check PP Spline & BSpline get the same curve with the same boundary conditions and knots" << std::endl;
    check_P4();
    std::cout << "figures are saved in output/check" << std::endl;
    std::cout << "Checking P5 : BSpline in any order" << std::endl;
    check_P5();
    std::cout << "figures are saved in output/check" << std::endl;
    
    return 0;
}