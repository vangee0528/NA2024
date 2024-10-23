#include "interpolation.h" 
#include <iostream>
#include <fstream>
#include <string>

int main() {
    // 参数设置
    int ms[] = {10, 40, 160};
    
    
    for (int i = 0; i < 3; i++) {

        std::ofstream outfile;
        std::string filename = "data/ProblemF_" + std::to_string(i + 1) + ".txt";
        outfile.open(filename);

        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return 1; 
        }

        std::vector<Point> Points = find_points(ms[i]);
        std::vector<std::vector<Point>> control_points = convert_to_control_points(Points);

        for (int j = 0; j < control_points.size(); j++) {
            Bezier bezier(control_points[j]);
            bezier.FileOut(outfile, Points[j].x);
        }

        outfile.close();
    }

    return 0;
}
