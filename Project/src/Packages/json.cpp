#include "json.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

// 辅助函数：去除字符串两端的空白字符
std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

// 辅助函数：解析字符串值
std::string parse_string(const std::string& s) {
    size_t first_quote = s.find("\"");
    size_t last_quote = s.rfind("\"");
    if (first_quote == std::string::npos || last_quote == std::string::npos || last_quote <= first_quote) {
        throw std::runtime_error("Invalid string format: " + s);
    }
    return s.substr(first_quote + 1, last_quote - first_quote - 1);
}

// 辅助函数：解析数组
std::vector<double> parse_array(const std::string& s) {
    std::vector<double> array;
    size_t start = s.find("[");
    size_t end = s.find("]", start);
    if (start == std::string::npos || end == std::string::npos) {
        throw std::runtime_error("Invalid array format: " + s);
    }
    std::string array_content = s.substr(start + 1, end - start - 1);
    std::stringstream ss(array_content);
    std::string item;
    while (std::getline(ss, item, ',')) {
        array.push_back(std::stod(trim(item)));
    }
    return array;
}

SplineParameters read_json(const std::string& file_path) {
    std::ifstream input_file(file_path);
    if (!input_file.is_open()) {
        throw std::runtime_error("无法打开文件: " + file_path);
    }

    SplineParameters params;
    std::string line;
    while (std::getline(input_file, line)) {
        // 移除注释
        size_t comment_pos = line.find("//");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 查找键值对
        size_t colon_pos = line.find(":");
        if (colon_pos == std::string::npos) continue;

        std::string key = trim(line.substr(0, colon_pos));
        key = parse_string(key);  // 提取键名

        std::string value = trim(line.substr(colon_pos + 1));
        // 移除可能的逗号
        if (!value.empty() && value.back() == ',') {
            value.pop_back();
            value = trim(value);
        }

        if (value.empty()) continue;

        if (value[0] == '\"') {  // 字符串
            if (key == "spline_type") params.spline_type = parse_string(value);
            else if (key == "method") params.method = parse_string(value);
            else if (key == "boundary_condition") params.boundary_condition = parse_string(value);
        }
        else if (value[0] == '[') {  // 数组
            if (key == "interval") {
                params.interval = parse_array(value);
            }
            else if (key == "time_points") {
                params.time_points = parse_array(value);
            }
            else if (key == "coefficients") {
                params.coefficients = parse_array(value);
            }
        }
        else {  // 数字
            if (key == "dimension") {
                params.dimension = std::stoi(value);
            }
            else if (key == "order") {
                params.order = std::stoi(value);
            }
            else if (key == "num_intervals") {
                params.num_intervals = std::stoi(value);
            }
            else if (key == "da") {
                params.da = std::stod(value);
            }
            else if (key == "db") {
                params.db = std::stod(value);
            }
        }
    }

    return params;
}