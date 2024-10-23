#include <iostream>
#include <fstream>
#include <vector>
#include "interpolation.h"

int main() {
    // 输入数据
    std::vector<double> days = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> sp1_weights = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> sp2_weights = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    // 计算牛顿插值多项式
    NewtonInterpolator sp1_interpolator(days, sp1_weights);
    NewtonInterpolator sp2_interpolator(days, sp2_weights);

    Polynomial sp1_polynomial = sp1_interpolator.interpolate();
    Polynomial sp2_polynomial = sp2_interpolator.interpolate();

    // 打开文件输出
    std::ofstream outfile("data/ProblemE.txt");

    // 输出样本1的多项式
    outfile << "Sample 1 Weight-Time Polynomial:\n";
    for (int i = sp1_polynomial.coefficients.size() - 1; i >= 0; --i) {
        if (i < sp1_polynomial.coefficients.size() - 1 && sp1_polynomial.coefficients[i] >= 0) {
            outfile << "+";
        }
        outfile << sp1_polynomial.coefficients[i] << "x**" << i << " ";
    }
    outfile << "\n";

    // 输出样本2的多项式
    outfile << "Sample 2 Weight-Time Polynomial:\n";
    for (int i = sp2_polynomial.coefficients.size() - 1; i >= 0; --i) {
        if (i < sp2_polynomial.coefficients.size() - 1 && sp2_polynomial.coefficients[i] >= 0) {
            outfile << "+";
        }
        outfile << sp2_polynomial.coefficients[i] << "x**" << i << " ";
    }
    outfile << "\n";

    // 输出给定10个时间点的重量
    outfile << "Day Sample1_Weight Sample2_Weight\n";
    for (int i = 0; i < days.size(); ++i) {
        outfile << days[i] << " " << sp1_weights[i] << " " << sp2_weights[i] << "\n";
    }

    // (a) 在 t=43 时预测重量
    double futureTime = 43.0;
    double sp1_future_weight = sp1_polynomial.evaluate(futureTime);
    double sp2_future_weight = sp2_polynomial.evaluate(futureTime);

    // 输出 t=43 时的预测重量
    std::cout << "Predicted Weight of Sample 1 at t = " << futureTime << ": " << sp1_future_weight << " grams" << std::endl;
    std::cout << "Predicted Weight of Sample 2 at t = " << futureTime << ": " << sp2_future_weight << " grams" << std::endl;

    // (b) 预测幼虫是否会死亡
    if (sp1_future_weight <= 0) {
        std::cout << "\033[1;31mSample 1 larvae will die after another 15 days.\033[0m" << std::endl;
    } else {
        std::cout << "\033[1;32mSample 1 larvae will survive after another 15 days.\033[0m" << std::endl;
    }

    if (sp2_future_weight <= 0) {
        std::cout << "\033[1;31mSample 2 larvae will die after another 15 days.\033[0m" << std::endl;
    } else {
        std::cout << "\033[1;32mSample 2 larvae will survive after another 15 days.\033[0m" << std::endl;
    }

    std::cout << "----------------------------------------" << std::endl;

    // 关闭文件
    outfile.close();

    return 0;
}