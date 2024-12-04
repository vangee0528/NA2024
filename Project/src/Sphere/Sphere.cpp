#include "../Packages/spline.h"
#include <iostream>
#include <vector>
#include <cmath>

int main(){
    // 在球面x^2+y^2+(z-1)^2=1上随机选取10个点
    try{
    std::vector<std::vector<double>> spherical_points(3);
    system("mkdir -p output/Sphere");
    freopen("output/Sphere/output.txt","w",stdout);
    for (int i = 0; i < 10; ++i) {
        double theta = 2 * M_PI * rand() / RAND_MAX;
        double phi = M_PI * rand() / RAND_MAX;
        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi) + 1;
        spherical_points[0].push_back(x);
        spherical_points[1].push_back(y);
        spherical_points[2].push_back(z);
    }
    SplineOnSphere spline_on_sphere(spherical_points, 2, CLAMPED, 0, 0);
    spline_on_sphere.print();
    fclose(stdout);
    }
    catch(const char* msg){
        std::cerr << msg << std::endl;
    }

}