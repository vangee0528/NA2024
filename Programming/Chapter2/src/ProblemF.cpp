#include "Bezier.h"
#include <iostream>

int main(){
    int ms[] = {10};
    for(int i = 0; i < 1; i++){
        std::vector<Point> Points = find_points(ms[i]);
        std::vector<std::vector<Point>> control_points = convert_to_control_points(Points);

        for(int j = 0; j < control_points.size(); j++){
            Bezier bezier(control_points[j]);
            
            bezier.printOut();
        }
    }

}