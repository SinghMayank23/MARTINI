// ImportLRates.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in transition rates for radiative and elastic processes from files

#ifndef ImportLRates_h
#define ImportLRates_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#define ABS(x)      ( (x) > 0 ? (x) : -(x) )

class ImportLRates
{
 private:
  
  //FILE *openfile ( char filename[] , char mode[] );
  
  /* the driver produces the data line per line.
     p and k are in GeV, L is in fermis, and gam is the dimensionless
     rate  dGamma/dk / g^4. */
  typedef struct {
    double p;
    double L;
    double *k;
    double *gam[4];
  } lineofdata;
  
  /* data structure I use: the k of every point is stored  */
  /* data is ordered as "dat[p][L].gam[k]" */
  typedef struct {
    int nlines, nk;
    double T;  /* This is some order-of-magnitude estimate of the temperature.
		  It is used to strip Bose factors in the interpolation routine and
		  thereby improve fits.
		  The actual temperature profile is time-dependent in fact but this T
		  is constant over all one profile.
	       */
    int profile;
    lineofdata *dat;
  } rawdata;

  
  typedef struct
  {
    double df;
    double da;
    double cf;
    double ca;
    int    Nc;
    int    Nf;
    int    Bethe_Heitler;
    int    BDMPS;
    int    include_gluons;
    
    double alpha_s;
    double alpha;
    double delta_x;
    double dx_max;
    double dp;
    double p_min;
    double p_max;
    int   n_p;
    long   n_pmin;
    double k_min;
    double k_max;
    int   n_k;
    long   n_kmin;
    
  } Gamma_info;

  Gamma_info   dat;
  rawdata      Gam;

  // nlines should be an upper bound on the number of possible data strings
  void write_header(FILE *wfile, Gamma_info *dat, int profile, double T, int nlines);
  
  void write_one_line(FILE *out, int nk, lineofdata line);

  // Perform the interpolation with respect to k in a single data line.
  // This is only one part of the interpolation w/r to p,L,k.
  double interpolate_k(rawdata *raw, int p, double k, int process);
  

 public:
  ImportLRates();//constructor

  ~ImportLRates();//destructor

  int read_raw(Gamma_info *dat, rawdata *raw);
  void free_raw(rawdata *raw);
  double use_raw ( double p , double k , double L, rawdata *raw, int process);
  void init();
  double getRate(double p, double k, double L);
  double getRate_em(double p, double k, double L);
  double getRate_gqq(double p, double k, double L);
  double getRate_ggg(double p, double k, double L);

};
#endif
