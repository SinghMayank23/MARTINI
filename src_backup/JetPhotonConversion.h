// JetPhotonConversion.h is a part of the MARTINI event generator.
// Copyright (C) 2020 Rouzbeh Yazdi (?)
#ifndef JETPHOTONCONVERSION_H
#define JETPHOTONCONVERSION_H
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <Basics.h>
/*
*   @author Rouz
*   ConversionPhotonGrid:
*       - Class to hold information on the conversion photon grid
*           as well as the grid itself. During a MARTINI event run,
*           the ConversionPhotonGrid object can be used to interpolate
*           for a 
*/
class ConversionPhotonGrid
{
    private:
        double*** grid;
        int jet_energy_grid_size;
        int temperature_grid_size;
        int photon_energy_grid_size;
        double jet_energy_min;
        double jet_energy_max;
        double dEjet;
        double temperature_min;
        double temperature_max;
        double dT;
        double photon_energy_min;
        double photon_energy_max;
        double dEgamma;
    public:
       ConversionPhotonGrid(std::string gridinfo, std::string gridFile);
       ~ConversionPhotonGrid();
       double getConvRate_Trilinear(double eJet, double T, double eGamma);
};
/*
**
*   ConversionAngularGrid:
*       - Class to hold the information for conversion photon rate
*         differential in photon kT and kz as well as jet momentum and 
*         local temperature.
**
*/
class ConversionAngularGrid
{
    private: 
        double**** grid;
        int jet_momentum_grid_size;
        int temperature_grid_size;
        int transverse_k_grid_size;
        int longitudinal_k_grid_size;
        double jet_mom_min;
        double jet_mom_max;
        double dPjet;
        double temperature_min;
        double temperature_max;
        double dT;
        double kT_min;
        double kT_max;
        double dkT;
        double kZ_min;
        double kZ_max;
        double dkZ;
        
        double integrand_kz(double kz, double phiGamma, double jetx, double jety, double jetz, double kT, double vx, double vy, double vz, double T);
        double integrate_kz(double phiGamma, double eta, double jetx, double jety, double jetz, double kT, double vx, double vy, double vz, double T);
    public:
        ConversionAngularGrid(std::string grid_info_file, std::string grid_data_file);
        ~ConversionAngularGrid();
        double interpolate(double jet_momentum, double T, double kT, double kZ);
        double convert(double jetx, double jety, double jetz, double kT, double eta, double vx, double vy, double vz, double T);
};
#endif 
