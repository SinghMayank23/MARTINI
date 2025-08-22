// Basics.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains basic structures and constants

#ifndef Basics_h //avoids multiple inclusions of the header file
#define Basics_h

using namespace std;

static const double PI = 3.141592653589793;
static const double hbarc = 0.197327053;
static const int Nc = 3;

struct ReturnValue 
{
  double x;
  double y;
  int rejections;
  int acceptances;
};

struct HydroInfo
{
  double T;
  double QGPfrac;
  double vx;
  double vy;
  double vz;
  double veta;
};
  
struct HydroInfoViscous
{
  double Wtautau;
  double Wtaux;
  double Wtauy;
  double Wtaueta;
  double Wxx;
  double Wxy;
  double Wxeta;
  double Wyy;
  double Wyeta;
  double Wetaeta;
  double BulkPi;
  double Entropy;
  double cs2;
};
  
struct Norms 
{
  double Gamma;
  double Gamma_gqq;
  double Gamma_ggg;
  double Gamma_em;
};

#endif
