#ifndef binarycollision_info_H
#define binarycollision_info_H


#include <fstream>
#include <cstring>

using namespace std;

class binarycollision_info
{
    private:
        int grid_nx, grid_ny;
        double grid_dx, grid_dy;
        double grid_x_min, grid_y_min;
        double **rho_binary;
        double rho_binary_max;
        string filename;
        int initialization_status;

    public:
        binarycollision_info();
        ~binarycollision_info();

        void init(string filename, int kind);
        double get_binary_collision_density(double x, double y);
        void get_sample(double *x_picked, double *y_picked);
        void check_samples();
        void print_info();

};

#endif
