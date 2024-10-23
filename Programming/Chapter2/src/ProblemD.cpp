#include <iostream>
#include <fstream>
#include <vector>
#include "interpolation.h"

int main() {
    
    // 参数设置
    std::vector<double> time = {0, 3, 5, 8, 13}; 
    std::vector<double> displacement = {0, 225, 383, 623, 993}; 
    std::vector<double> speed = {75, 77, 80, 74, 72}; 

    // 计算Hermite插值多项式
    HermiteInterpolator hermiteInterpolation(time, displacement, speed);
    Polynomial p = hermiteInterpolation.interpolate(); 
    Polynomial dp = p.derivative(); 

    // 打开文件输出
    std::ofstream outfile("data/ProblemD.txt");

    // 输出距离-时间多项式
    outfile << "Displacement-Time Polynomial:\n";
    for (int i = p.coefficients.size() - 1; i >= 0; --i) {
        if (i < p.coefficients.size() - 1 && p.coefficients[i] >= 0) {
            outfile << "+";
        }
        outfile << p.coefficients[i] << "x**" << i << " ";
    }
    outfile << "\n";

    // 输出速度-时间多项式
    outfile << "Speed-Time Polynomial:\n";
    for (int i = dp.coefficients.size() - 1; i >= 0; --i) {
        if (i < dp.coefficients.size() - 1 && dp.coefficients[i] >= 0) {
            outfile << "+";
        }
        outfile << dp.coefficients[i] << "x**" << i << " ";
    }
    outfile << "\n";

    // 输出给定的时间点的位置和速度，用于比较拟合结果
    outfile << "Time Displacement Speed\n";
    for (int i = 0; i < time.size(); ++i) {
        outfile << time[i] << " " << displacement[i] << " " << speed[i] << "\n";
    }

    // (a) 在 t=10 时预测位置和速度
    double testTime = 10.0;
    double estimatedDisplacement = p.evaluate(testTime);
    double estimatedSpeed = dp.evaluate(testTime);

    // 输出 t=10 时的位置和速度
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "\033[1;32mEstimated Displacement at t = " << testTime << ": " << estimatedDisplacement << " feet\033[0m" << std::endl;
    std::cout << "\033[1;32mEstimated Speed at t = " << testTime << ": " << estimatedSpeed << " feet/second\033[0m" << std::endl;

    // (b) 确定汽车是否曾经超过速度限制
    bool exceedsSpeedLimit = false;
    double exceedTime = 0.0;
    double exceedSpeed = 0.0;
    for (double t = time.front(); t <= time.back(); t += 0.1) {
        double spd = dp.evaluate(t);
        if (spd > 81) {
            exceedsSpeedLimit = true;
            exceedTime = t;
            exceedSpeed = spd;
            break;
        }
    }

    if (exceedsSpeedLimit) {
        std::cout << "\033[1;31mThe car exceeds the speed limit of 81 feet/second at t = " << exceedTime << " seconds.\033[0m" << std::endl;
        std::cout << "\033[1;31mSpeed at that time: " << exceedSpeed << " feet/second\033[0m" << std::endl;
    } else {
        std::cout << "\033[1;32mThe car never exceeds the speed limit of 81 feet/second.\033[0m" << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;

    outfile.close();

    return 0;
}