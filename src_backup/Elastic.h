// Elastic.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the elastic processes for energy loss and momentum broadening

#ifndef Elastic_h
#define Elastic_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include <gsl/gsl_sf_lambert.h> // includes the Lambert function from the gnu scientific library
#include "Basics.h"
#include "Random.h"
#include "Import.h"
#include "Pythia8/Pythia.h"


using namespace Pythia8;

class Elastic
{
 private:

 public:
  Elastic(){};//constructor
  ~Elastic(){};//destructor
  double totalRate(double p, double T, double alpha_s, int Nf, int process);
  double totalRatePos(double p, double T, double alpha_s, int Nf, int process);
  double totalRateNeg(double p, double T, double alpha_s, int Nf, int process);
  double function(double p, double omega, double alpha_s, Import *import, int process);
  double functionOmegaQ(double p, double omega, double q, double alpha_s, Import *import, int process);
  double area (double y, double alpha_s, int posNegSwitch, int process);
  double areaOmegaQ (double y, double omega, double alpha_s, int process);
  double areaOmegaQ2 (double y, double omega, double alpha_s, int process);
  double findValuesRejection(double u, double T, double alpha_s, int Nf,
		Random *random, Import *import, int process);
  double findValuesRejectionOmegaQ(double p, double omega, double T, double alpha_s, int Nf,
				   Random *random, Import *import, int process);
  Vec4 getNewMomentum(Vec4 vecpRest, double omega, double q, Random *random);
  // sample thermal parton which is going to be recoiled by momentum transfer
  // this thermal parton should be constrained in a way that 
  // new recoil particle become on shell. added by CP
  Vec4 getRecoilMomentum(Vec4 vecq, double temp, int kind, Random *random);
  Vec4 getThermalMomentum(Vec4 vecq, double temp, int kind, Random *random);
};

#endif
