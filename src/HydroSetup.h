// HydroSetup.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions that return interpolated data at a given space-time point

#ifndef HydroSetup_h
#define HydroSetup_h

#include<vector>
#include<cstring>
#include "HydroCell.h"
#include "Basics.h"

using namespace std;

class HydroSetup
{
    private:
        double hydroTau0;        // tau_0 in the hydro data files
        double hydroTauMax;      // tau_max in the hydro data files
        double hydroDtau;        // step dtau in fm/c in the hydro data files
        double hydroXmax;        // maximum x in fm in the hydro data files 
                                 // [-xmax, +xmax] for both x and y
        double hydroZmax;        // maximum z in fm in the hydro data files 
                                 // [-zmax, +zmax] for 3D hydro
        double hydroDx;          // step dx in fm in the hydro data files
        double hydroDz;          // step dz in fm in the hydro data files in 
                                 // the z-direction for 3D hydro
                                 
        double hydroTfinal;      // temperature at which jet energy loss stops

        int hydroWhichHydro;     // choose a hydro evolution model to use
        int use_tau_eta_coordinate; 

        bool boost_invariant;

        int itaumax, ixmax, ietamax;

        vector<HydroCell> *lattice; // array to store hydro information

    public:
        HydroSetup();   //constructor
        ~HydroSetup();  //destructor
        
        double get_hydro_tau_max() {return hydroTauMax;};
        // RMY Additiosn begin: 
        // Usage aim: MARTINI::fragmentation method 2 
        // Add helper functions for Hydro Info retrieval
        double get_hydro_Dx() { return hydroDx;};
        double get_hydro_Dz() { return hydroDz;};
        double get_hydro_Dtau() { return hydroDtau;};
        int get_ixmax() { return ixmax;};
        int get_itaumax() { return itaumax;};
        int get_izmax() { return ietamax;}; // z is actually eta. 
        // end RMY additions
        void readHydroData(double tau0, double taumax, double dtau, 
		     double xmax, double zmax, double dx, double dz, 
                 int nskip_tau, int nskip_x, int nskip_z,
                 int whichHydro, double Tfinal, int taueta_coord,
                 string evolution_name);
        
        HydroInfo getHydroValues(double x, double y, double z, double t);
        
        void output_temperature_evolution(string filename_base);
        void update_grid_info(
                double tau0, double tau_max, double dtau,
                double x_max, double dx, double z_max, double dz);
        void output_temperature_evolution_list_format(string filename_base);
};

#endif
