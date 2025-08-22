#ifndef Import_h
#define Import_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;

#define ABS(x)      ( (x) > 0 ? (x) : -(x) )
#define PauliBlock(x) ( (x)>0 ? 1/(1+exp(-(x))) : exp(x)/(exp(x)+1) )
#define BoseStim(x)   ( (x)>0 ? 1/(1-exp(-(x))) : exp(x)/(exp(x)-1) )
#define READ_LETTER(x,y) while ( (( (x) = getc (y) ) < 'a' ||       \

#define NN 312
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
  
 public:
  Import(){};//constructor
  ~Import(){};//destructor
  double use_table ( double p , double k , double dGamma[NP][NK] , int which_kind );
  void prep_dGamma ( Gamma_info *dat , dGammas *gam , 
		     double T , double beta , double cos_phi );
  void list_dGamma ( Gamma_info *dat , dGammas *gam , 
		     double T , double beta , double cos_phi );
  void read_table  ( Gamma_info *dat , dGammas *Gam );
  void init();
  void show_dGamma();

  double getRate(double p, double k);
  double getRate_gqq(double p, double k);
  double getRate_ggg(double p, double k);

  int get_NP(){return(NP);};
  int get_NK(){return(NK);};

};

#endif
