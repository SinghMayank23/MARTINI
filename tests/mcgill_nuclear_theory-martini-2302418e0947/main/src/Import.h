
// Import.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in transition rates for radiative and elastic processes from files

#ifndef Import_h
#define Import_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#define ABS(x)      ( (x) > 0 ? (x) : -(x) )
#define PauliBlock(x) ( (x)>0 ? 1/(1+exp(-(x))) : exp(x)/(exp(x)+1) )
#define BoseStim(x)   ( (x)>0 ? 1/(1-exp(-(x))) : exp(x)/(exp(x)-1) )
#define READ_LETTER(x,y) while ( (( (x) = getc (y) ) < 'a' ||       \

#define NNM 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

class Import
{
 private:

  static const int NP = 230;
  static const int NK = 381;
  static const int Nf = 3;
  static const int Nalphas = 11; //for elastic only
  static const int Nomega = 120; //for elastic only
  static const int Nq = 60;      //for elastic only
  static const double omegaStep = 0.2; //for elastic only
  static const double qStep = 0.2; //for elastic only
  static const double alphaMin = 0.15;
  static const double alphaStep = 0.03;


  typedef struct
  {
    double ddf;
    double dda;
    double dcf;
    double dca;
    int    include_gluons;
    int Nc;
    int Nf;
    int BetheHeitler;
    int BDMPS;
    double alpha_s;
    double alpha;
    double delta_x;
    double dx_max;
    double dp;
    double p_min;
    double p_max;
    long   n_p;
    long   n_pmin;
    double k_min;
    double k_max;
    long   n_k;
    long   n_kmin;
    
  } Gamma_info;
  
/// Structure to store information about splitting functions, to interpolate from
  typedef struct
  {
    double dGamma[NP][NK];
    double dGamma_gqq[NP][NK];
    double dGamma_ggg[NP][NK];
    double dGamma_em[NP][NK];
    double tau[NP][NK];
    double tau_gqq[NP][NK];
    double tau_ggg[NP][NK];
    double tau_em[NP][NK];
  } dGammas;
  
  
  Gamma_info	        dat;
  dGammas		Gam;

  double dGamma[NP*NK];
  double dGamma_em[NP*NK];
  double dGamma_gqq[NP*NK];
  double dGamma_ggg[NP*NK];

  double *dGamma_qq;
  double *dGamma_qg;
  double *dGamma_qq_q;
  double *dGamma_qg_q;
  
 public:
  Import();//constructor

  ~Import();//destructor
  double use_table ( double p , double k , double dGamma[NP][NK] , int which_kind );
  double use_elastic_table ( double omega , double alpha_s, int which_kind );
  double use_elastic_table_omega_q ( double omega , double q, double alpha_s, int which_kind );
  void prep_dGamma ( Gamma_info *dat , dGammas *gam , 
		     double T , double beta , double cos_phi );
  void list_dGamma ( Gamma_info *dat , dGammas *gam , 
		     double T , double beta , double cos_phi );
  void read_table  ( Gamma_info *dat , dGammas *Gam, int rateSelector);
  void init(int rateSelector);
  void show_dGamma();
  void readElasticRate();
  void readElasticRateOmegaQ();

  double getRate(double p, double k);
  double getRate_gqq(double p, double k);
  double getRate_ggg(double p, double k);
  double getRate_em(double p, double k);
  double getElasticRate(double p, double omega, double alpha_s, int process);
  double getElasticRateOmegaQ(double p, double omega, double q, double alpha_s, int process);

  int get_NP(){return(NP);};
  int get_NK(){return(NK);};

};
#endif
