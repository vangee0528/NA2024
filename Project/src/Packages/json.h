#ifndef JSON_H
#define JSON_H

#include <string>
#include <vector>

struct SplineParameters {   
    std::string spline_type;                        // 样条类型
    int dimension;                                  // 维数(参数方程数量)
    int order;                                      // 阶数
    std::string method;                             // 选点方法
    std::vector<double> interval;                   // 区间
    int num_intervals;                              // 区间数
    std::vector<double> time_points;                // 自定义选点
    std::vector<double> coefficients;               // 系数
    std::string boundary_condition;                 // 边界条件
    double da;                                      // 左边界导数
    double db;                                      // 右边界导数
};  

SplineParameters read_json(const std::string& file_path);

#endif // JSON_H