#ifndef Evolution_h //avoids multiple inclusions of the header file
#define Evolution_h
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>


#include "Constants.h"
#include "tools/FileNames.h"

#define ABS(x)      ( (x) > 0 ? (x) : -(x) )
#define PauliBlock(x) ( (x)>0 ? 1/(1+exp(-(x))) : exp(x)/(exp(x)+1) )
#define BoseStim(x)   ( (x)>0 ? 1/(1-exp(-(x))) : exp(x)/(exp(x)-1) )
#define READ_LETTER(x,y) while ( (( (x) = getc (y) ) < 'a' ||       \
                                    (x) > 'z' ) && ( (x) < 'A' || (x) > 'Z' ) )    

using namespace std;

/// Includes all the time evolution functions
class Evolution
{
 private:
  
  double LogEmax;
  double LogEmin;
  double LogStepE;
  double LogOmegaMin;
  double LogStepOmega;
  double Tfile;
  double stepT;

  double PI;
  double hbarc;
  double EulerGamma;
  /// cb = -EulerGamma + Zeta'(2)/Zeta(2) + Log(2)
  double cb; 
  /// cf = -EulerGamma + Zeta'(2)/Zeta(2)
  double cf;
  double cs;
  
  static const int NP = 230;
  static const int NK = 381;
  
  int Nc;


/// Structure to store information about couplings, momentum range and discretization, group theoretic factors, etc 
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
  
  double T;
  double dt;
  double alphas;
  double dE;
  int Nf;
  double initialE;
  double width;
  double maxEnergy;
  int BetheHeitler;
  int BDMPS;
  int doRad;
  int doCol;
  int newRadGamma;
  int counter;
  
  double		*dGamma, *dGamma_gqq, *dGamma_ggg, *dGamma_em;
  Gamma_info	        dat;
  dGammas		Gam;
  double	        *P_spare;
  double		AA, tau, x; 
  double		pinitial, pginitial, pfinal, Rperd;
  double		taui, Ti, tauc, Tc, rd, tauh, fqgp;
  double 		p;
  
  int			num_step, n_dx; 
  int			ip, np, np_g;
  char		        ch, flag;
  char		        filename[100], filename_q[100], filename_g[100], filename_e[100];

  const FileNames* evolutionsFileNameList;

 public:
  Evolution(const FileNames* fileList, const Constants* constants, double dt, double dE, double T, const double alphas, const int Nf, const double initialE, const double width, const double maxEnergy, const int BetheHeitler, const int BDMPS, const int doRad, const int doCol, const int newRadGamma, const int size=1, const int rank=0);
  ~Evolution();
/// determines K(z) = K_0(z) + ln(z/2) + EulerGamma
/// Sangyong: This is modified to include Peter's interpolation.
/// Original K_0 + ln + gamma_E is renamed to BK

  double K(double z, int just_K0);
  
  double BK(double z, int just_K0);

/// determines the BDMPS approximation to the above: 0.5*(1-z*K_1(z)), with K_1 the mofified Bessel fct
  double K_BDMPS(double z);
/// perform a 4th order Runge Kutta step to get dGamma:
  void RK_step (double f[4], double A, double D, double *b, double *k,  
		double C[3], double r[3] ); 
/// solves differential equation for f(b->0), which is the integral over F we need:
  double find_f_atzero ( double A, double D, double C[3], double r[3],
			 double *one_over_dE ); 
/// Solves the integral equations to give the production rate at p,k:
  double solve_int_eqn ( double p, double k,  
			 double df, double da, double cf, double ca, 
			 int n_quark_flavors, int p_is_gluon, 
			 int all_are_gluons, int emit_photon, 
			 double *one_over_deltaE); 
/// make table with radiative dGamma and save temporarily to a file:
  void build_table (Gamma_info * dat, dGammas *Gam);
/// write the table in binary to a file:
  void write_table ( Gamma_info *dat , dGammas Gam );   
/// uses looup table and simple interpolation to get the value of dGamma/dkdt at some value of p,k
  double use_table ( double p, double k, double dGamma[NP][NK], 
		     double tau[NP][NK], double *inv_dE, int which_kind );
/// find Gamma:
  void prep_dGamma ( Gamma_info *dat , dGammas *gam , 
		     double T , double beta , double cos_phi , 
		     double *dGamma , double *dGamma_gqq , 
		     double *dGamma_ggg , double *dGamma_em );
/// find dP for radiative loss:
double find_dP ( double *P , double *Pg , double *Pem , double *dGamma , 
		 double *dGamma_gqq , double *dGamma_ggg , 
		 double *dGamma_em , double *dP , 
		 double *dPg , double *dPem , Gamma_info dat );
/// evolve one time step:
double evolve_one_step ( double *P , double *Pg , double *Pem , 
			 double *dGamma , double *dGamma_gqq , 
			 double *dGamma_ggg , double *dGamma_em ,  
			 Gamma_info dat , double xmax );
/// this takes the input as I ave it and converts so radiative machinery can handle it:
void evolve_in_medium ( double *P , double *Pg , double *Pem , 
			double *dGamma , double *dGamma_gqq , 
			double *dGamma_ggg , double *dGamma_em ,  
			Gamma_info *dat , dGammas *gam , double T , 
			double beta , double cos_phi , double x_inFermi );
/// read im table:
  void read_table ( Gamma_info *dat , dGammas *Gam );
/// prepare table:
  void prepare_table ( Gamma_info *dat , dGammas *Gam );
/// convert parameters for radiative part:
  void prep_equipment ( Gamma_info *dat, dGammas *Gam, 
			double ** dGam1, double ** dGam2, 
			double ** dGam3, double ** dGam4 );
/// calculates mg^2 using T:
  double mgs(double T){return 0.5*(1.+static_cast<double>(Nf)/6.)*alphas*4.*PI*T*T;};
/// calculates mq^2 using T:
  double mqs(double T){return (1./6.)*alphas*4.*PI*T*T;};
/// energy loss rate dE/dt_qq for a hard quark scattering with a thermal quark:
  double dEdtqq(double T, double En);
/// energy loss rate dE/dt_qg for a hard quark scattering with a thermal gluon:
  double dEdtqg(double T, double En);
/// energy loss rate dE/dt_gq for a hard gluon scattering with a thermal quark:
  double dEdtgq(double T, double En);
/// energy loss rate dE/dt_gg for a hard gluon scattering with a thermal gluon:
  double dEdtgg(double T, double En);
/// Bose distribution:
  double fb(double omega, double t){return 1./(exp(omega/T)-1);};

/// transition rates:


  double Gammaqg(double T, double En);
  double Gammagq(double T, double En);

/*
  double Gammaqg2(double T, double En){return (fb(dE,T))*dEdtqg(T,En);};
  double Gammagq2(double T, double En){return (fb(dE,T))*dEdtgq(T,En);};
*/

  double Gammaqq1(double T, double En){return (1+fb(dE,T))*(dEdtqq(T,En)+dEdtqg(T,En));};
  double Gammaqq2(double T, double En){return (fb(dE,T))*(dEdtqq(T,En)+dEdtqg(T,En));};
  double Gammagg1(double T, double En){return (1+fb(dE,T))*(dEdtgg(T,En)+dEdtgq(T,En));};
  double Gammagg2(double T, double En){return (fb(dE,T))*(dEdtgg(T,En)+dEdtgq(T,En));};

/// smooth transition rates
  double GammaqqS(double T, double En, double omega, double* trqq, double* trqg);
  double GammaggS(double T, double En, double omega, double* trgq, double* trgg);

/// time evolution of quark distribution dPq/dt for delta function transition rate:
  double dPqdt(double T, double En, double* PqV, double* PgV);

/// time evolution of gluon distribution dPg/dt for delta function transition rate:
  double dPgdt(double T, double En, double* PqV, double* PgV);

/// time evolution of quark distribution dPq/dt for smooth transition rate:
  double dPqdtS(double T, double En, double* PqV, double* PgV, double* trqq, double* trqg, double* trgq, double* trgg);

/// time evolution of gluon distribution dPg/dt for smooth transition rate:
  double dPgdtS(double T, double En, double* PqV, double* PgV, double* trqq, double* trqg, double* trgq, double* trgg);

/// output of dEdtab:
  void dEdtOutput();

  void initP(double* Pq, double* Pg, double* Pem);

  void run(double* Pq, double* Pg, double* Pem, double* PqA, double* PgA);

  void runsmooth(double* Pq, double* Pg, double* Pem, double* PqA, double* PgA, double* trqq, double* trqg, double* trgq, double* trgg);

};
#endif
