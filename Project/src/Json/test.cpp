#include <iostream>
#include "../Packages/spline.h"

// 测试函数
double test_function1(double x) {
    return 1.0 / (1 + 25 * x * x);
}

double test_function2(double x) {
    return std::sin(x);
}

int main() {
    try {
        std::vector<MathFunction> functions = {MathFunction(test_function1), MathFunction(test_function2)};
        
        // 输出 BSpline 结果到文件
        FILE* bspline_file = freopen("output/jsontest/output_bspline.txt", "w", stdout);
        if (bspline_file == nullptr) {
            std::cerr << "Error: Unable to open output_bspline.txt for writing" << std::endl;
            return 1;
        }
        BSpline bspline("src/Json/input/input_bspline.json", functions);
        bspline.print();
        fclose(bspline_file);
        
        // 输出 PPSpline 结果到文件
        FILE* ppspline_file = freopen("output/jsontest/output_ppspline.txt", "w", stdout);
        if (ppspline_file == nullptr) {
            std::cerr << "Error: Unable to open output_ppspline.txt for writing" << std::endl;
            return 1;
        }
        PPSpline ppspline("src/Json/input/input_ppspline.json", functions);
        ppspline.print();
        fclose(ppspline_file);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}