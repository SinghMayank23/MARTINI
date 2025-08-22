
// HydroSetup.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions that return interpolated data at a given space-time point

#ifndef HydroSetup_h
#define HydroSetup_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <list>
#include "HydroCell.h"
#include "Basics.h"
#include "Pythia.h"

using namespace Pythia8;

class HydroSetup
{
 private:
  int a;
 public:
  HydroSetup(){};//constructor
  ~HydroSetup(){};//destructor
  void readHydroData(double A, double tau0, double taumax, double dtau, 
		     double xmax, double zmax, double dx, 
		     double dz, int whichHydro, int subset,  bool viscous, double b, vector<HydroCell> *lattice, string evolution_name);
  //readHydroData is overloaded for FIC. -CFY
  void readHydroData(double A, double tau0, double taumax, double dtau, 
		     double xmax, double zmax, double dx, 
		     double dz, int whichHydro, int file_number, int subset,  bool viscous, double b, vector<HydroCell> *lattice, string evolution_name);
  HydroInfo getHydroValues(double x, double y, double z, double t, double hydroXmax, double hydroZmax, double hydroTauMax, 
			   double hydroTau0, double hydroDx, double hydroDz, double hydroDtau, int hydroWhichHydro, int fixedDistribution, 
			   vector<HydroCell>* lattice, bool trackHistory=0);
  // getHydroValues2 is the old routine for the full 3+1D that was not optimized for memory usage using its symmetries
  HydroInfo getHydroValues2(double x, double y, double z, double t, double hydroXmax, double hydroZmax, double hydroTau0, 
			    double hydroDx, double hydroDz, double hydroDtau, int hydroWhichHydro, vector<HydroCell>* lattice);
};

#endif
