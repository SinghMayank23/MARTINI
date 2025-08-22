
// Glauber.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines for sampling the initial geometry according to the Glauber model

#ifndef Glauber_h //avoids multiple inclusions of the header file
#define Glauber_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <math.h>
#include "Pythia.h"
#include "Random.h"
#include "Basics.h"

#define TOL (1.0e-6)
#define tiny (1.0e-10)
#define limit 10000

class Glauber
{
 private:
  
  double cotan(double i) { return(1. / tan(i)); }
  double acotan(double i) { return( atan(1./i)); }

  typedef double (Glauber::*ptr_func)(double);

  typedef struct nucleus 
  {
    char *name;
    double A;
    double Z;
    int AnumFunc;
    int AnumFuncIntegrand;
    int DensityFunc;
    double w_WS;
    double a_WS;
    double R_WS;
    double rho_WS;
    
  } Nucleus;
    
  typedef struct data
  {
    double SigmaNN; 
    Nucleus Target;
    Nucleus Projectile;
    double SCutOff;
    int InterMax; 
    
    /* trap door */
  } Data;
  
  double AnumR, NuInS_S;
  Nucleus *Nuc_WS; 
  Data LexusData;  
  ptr_func tempFunc;
  double b; // impact parameter
  double gamma; // parameters of envelope function
  double height; // parameters of envelope function
  double currentTAB;
  double currentA;
  bool glauberEnvelope;
  
 public:
  double nucleusA() const {return currentA;} 
  int IsFile(char *file_name);
  void FindNucleusData(Nucleus *nucleus, char *target, char *file_name);
  void PrintLexusData();
  void PrintNucleusData(Nucleus *nucleus);
  int LinearFindXorg(double x, double *Vx, int ymax);
  double FourPtInterpolate(double x, double *Vx, double *Vy, double h, int x_org, int ymax);
  void MakeCoeff(double *a, double *b, double *c, double *d, 
		 double *Vy, double *Vx, double h, int x_org);
  double VInterpolate(double x, double *Vx, double *Vy, int ymax);
  int FindXorg(double x, double *Vx, int ymax);
  double *MakeVx(double down, double up, int maxi_num);
  double *MakeVy(char *, double *vx, int maxi_num);
  double *ReadInVx(char *, int maxi_num, int quiet);
  double *ReadInVy(char *, int maxi_num, int quiet);

  double InterNuPInSP(double s);
  double InterNuTInST(double s);
  void CalcRho(Nucleus *nucleus);
  double NuInS(double s);

  double Anum3Fermi(double R_WS);
  double Anum3FermiInt(double xi);
  double NuInt3Fermi(double xi);
  double Anum3Gauss(double R_WS);
  double Anum3GaussInt(double xi);
  double NuInt3Gauss(double xi);
  double Anum2HO(double R_WS);
  double Anum2HOInt(double xi);
  double NuInt2HO(double xi);
 
  char *StringFind(char *file_name, char *st);
  double DFind(char *file_name, char *st);
  int integer(double x);
  char *char_malloc(int n1);
  double *vector_malloc(int n1);
  void char_free(char *vec);
  double integral (int id, double down, double up, double tol, int *count);
  double qnc7(int id, double tol, double down, double dx, double *f_of, 
	      double pre_sum, double area, int *count);
  void ReWriteString(char *file_name, char *st, char *x);
  void FileCopy(char *in_file, char *out_file);
  void FileCat(char *in_file, char *out_file);
  double OLSIntegrand(double s);
  double TAB();
  double PAB(double x, double y);
  ReturnValue SamplePAB(Random *random);
  ReturnValue SamplePABRejection(Random *random);
  double areaPAB(double x, double A);
  double areaTA(double x, double A);
  ReturnValue SampleTARejection(Random *random);
  void preInit(double SigmaNN, string Target, string Projectile, double inb, int imax);
  void init(double SigmaNN, string Target, string Projectile, double b, int imax, bool glauberEnvelope);
};

#endif

