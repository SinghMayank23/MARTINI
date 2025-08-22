
// Rates.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains sampling routines for the radiative transition rates

#ifndef Rates_h
#define Rates_h

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

class Rates
{
 private:
  int rateSelector;
 public:
  Rates(int rateSelectorIn);//constructor
  ~Rates(){};//destructor
  /// calculates mq^2 using T:                                                                                                   
  double mqs(double T, double alpha_s){return (1./6.)*alpha_s*4.*PI*T*T;};
  // return the rate
  double function(double x, double y, Import *import, int process);
  // area under the exact rates
  Norms integratePos(double p, double T, double alpha_s, int Nf);
  Norms integrateNeg(double p, double T, double alpha_s, int Nf);
  Norms integrate(double p, double T, double alpha_s, int Nf);
  double integratePhoton(double p, double T, double alpha_s);
  // conversion rates
  double Gammaqg(double En, double T, double alpha_s);
  double Gammaqgamma(double En, double T, double alpha_s);
  double Gammagq(double En, double T, double alpha_s, int Nf);
  // area under the envelope
  double area (double y, double u, int posNegSwitch, int process, Import * import);
  // sample the rates
  ReturnValue findValuesRejection(double u, double T, double alpha_s, int Nf, 
				  Random *random, Import *import, int process);
};

#endif
