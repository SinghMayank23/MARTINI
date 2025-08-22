#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "binarycollision_info.h"

using namespace std;

binarycollision_info::binarycollision_info()
{
    srand48(time(NULL));  // set random seed to system time
    initialization_status = 0;
}

binarycollision_info::~binarycollision_info()
{
    if(initialization_status == 1)
    {
        // clean up the table
        for(int i = 0; i < grid_nx; i++)
            delete [] rho_binary[i];
        delete [] rho_binary;
    }
}

void binarycollision_info::init(string filename, int kind)
{
    initialization_status = 1;

    double sum = 0.0;   // integrated value
    double temp;
    rho_binary_max = 0.0;

    string dummy;
    ifstream fin(filename.c_str());

    // Chun's initial binary collision profile
    if(kind == 2)
    {
        double Ncoll;
	// read in profile
	// first read in information line
        // # ncoll = 234 nx = 301 ny = 301 dx = 0.05 dy = 0.05 xmin = -7.5 ymin = -7.5
        fin >> dummy >> dummy >> dummy >> Ncoll 
            >> dummy >> dummy >> grid_nx
            >> dummy >> dummy >> grid_ny
            >> dummy >> dummy >> grid_dx
            >> dummy >> dummy >> grid_dy
            >> dummy >> dummy >> grid_x_min 
            >> dummy >> dummy >> grid_y_min;
    
        rho_binary = new double* [grid_nx];
        for(int i = 0; i < grid_nx; i++)
            rho_binary[i] = new double [grid_ny];
        
        for(int i = 0; i < grid_nx; i++)
        {
            for(int j = 0; j < grid_ny; j++)
            {
                fin >> temp;
                if(temp < 0.0) temp = 0.0;
                if(rho_binary_max < temp) rho_binary_max = temp;
                rho_binary[i][j] = temp;
                sum += rho_binary[i][j];
            }
        }
    }
    // Scott's initial energy density profile
    else if(kind == 3)
    {
        // .0  0.0  0.0   n_eta=  4  nx= 200  ny= 200
        // deta=   5.0    dx=  0.098650     dy=  0.098650
        fin >> dummy >> dummy >> dummy
            >> dummy >> dummy >> dummy
            >> grid_nx >> dummy >> grid_ny
            >> dummy >> dummy >> dummy
            >> grid_dx >> dummy >> grid_dy;

        rho_binary = new double* [grid_nx];
        for(int i = 0; i < grid_nx; i++)
            rho_binary[i] = new double [grid_ny];

        for(int i = 0; i < grid_nx; i++)
        {
            for(int j = 0; j < grid_ny; j++)
            {
                fin >> dummy >> dummy >> dummy
                    >> temp >> dummy >> dummy
                    >> dummy >> dummy >> dummy
                    >> dummy >> dummy;
                if(temp < 0.0) temp = 0.0;
                if(rho_binary_max < temp) rho_binary_max = temp;
                rho_binary[i][j] = temp;
                sum += rho_binary[i][j];
            }
        }
        grid_x_min = -1.*grid_nx*grid_dx/2.;
        grid_y_min = -1.*grid_ny*grid_dy/2.;
    }
    // IP-GLASMA: u_field with pi , RMY Oct 27, 2021
    // Difference between this and option #3: more dummies
    // to read after we get our rho.
    else if (kind == 4)
    {
        // .0  0.0  0.0   n_eta=  4  nx= 200  ny= 200
        // deta=   5.0    dx=  0.098650     dy=  0.098650
        fin >> dummy >> dummy >> dummy
            >> dummy >> dummy >> dummy
            >> grid_nx >> dummy >> grid_ny
            >> dummy >> dummy >> dummy
            >> grid_dx >> dummy >> grid_dy;

        rho_binary = new double* [grid_nx];
        for(int i = 0; i < grid_nx; i++)
            rho_binary[i] = new double [grid_ny];

        for(int i = 0; i < grid_nx; i++)
        {
            for(int j = 0; j < grid_ny; j++)
            {
                fin >> dummy >> dummy >> dummy
                    >> temp  >> dummy >> dummy
                    >> dummy >> dummy >> dummy
                    >> dummy >> dummy >> dummy
                    >> dummy >> dummy >> dummy
                    >> dummy >> dummy >> dummy;
                if(temp < 0.0) temp = 0.0;
                if(rho_binary_max < temp) rho_binary_max = temp;
                rho_binary[i][j] = temp;
                sum += rho_binary[i][j];
            }
        }
        grid_x_min = -1.*grid_nx*grid_dx/2.;
        grid_y_min = -1.*grid_ny*grid_dy/2.;
    }

    sum = sum*grid_dx*grid_dy;

    // normalized to 1
    for(int i = 0; i < grid_nx; i++)
    {
        for(int j = 0; j < grid_ny; j++)
        {
            rho_binary[i][j] = rho_binary[i][j]/sum;
        }
    }
    rho_binary_max = rho_binary_max/sum;

    print_info();
}

void binarycollision_info::print_info()
{
    cout << "===============================================" << endl;
    cout << "grid information for binary collision profile:" << endl;
    cout << "nx = " << grid_nx << ", ny = " << grid_ny << endl;
    cout << "dx = " << grid_dx << ", dy = " << grid_dy << endl;
    cout << "x_min = " << grid_x_min << ", y_min = " << grid_y_min << endl;
    cout << "===============================================" << endl;
}

double binarycollision_info::get_binary_collision_density(double x, double y)
{
    int idx_x = (int)((x - grid_x_min)/grid_dx);
    int idx_y = (int)((y - grid_y_min)/grid_dy);

    // avoid underflow
    if(idx_x < 0) return(0.0);
    if(idx_y < 0) return(0.0);

    // avoid overflow
    if(idx_x > grid_nx - 2) return(0.0);
    if(idx_y > grid_ny - 2) return(0.0);

    double frac_x = (x - grid_x_min - idx_x*grid_dx)/grid_dx;
    double frac_y = (y - grid_y_min - idx_y*grid_dy)/grid_dy;
    double f1 = rho_binary[idx_x][idx_y];
    double f2 = rho_binary[idx_x+1][idx_y];
    double f3 = rho_binary[idx_x+1][idx_y+1];
    double f4 = rho_binary[idx_x][idx_y+1];

    double rho_binary_interp = (f1*(1. - frac_x)*(1. - frac_y) 
        + f2*frac_x*(1. - frac_y) + f3*frac_x*frac_y + f4*(1. - frac_x)*frac_y);
    return(rho_binary_interp);
}

void binarycollision_info::get_sample(double *x_picked, double *y_picked)
{
    // perform sampling according to rejection method

    double sample_x, sample_y;
    double rho_binary_sampled, random_test;
    do
    {
        sample_x = (drand48() - 0.5)*2.*fabs(grid_x_min);
        sample_y = (drand48() - 0.5)*2.*fabs(grid_y_min);
        rho_binary_sampled = get_binary_collision_density(sample_x, sample_y);
        random_test = drand48()*rho_binary_max;
    }while(random_test > rho_binary_sampled);

    *x_picked = sample_x;
    *y_picked = sample_y;
}

void binarycollision_info::check_samples()
{
    // generate 10,000 samples for checking purpose
    ofstream check_file("check_samples.dat");
    for(int i = 0; i < 10000; i++)
    {
        double x, y;
        get_sample(&x, &y);
        check_file << scientific << setprecision(8) << setw(18)
                   << x << "   " << y << endl;
    }
    check_file.close();
    exit(0);
}
