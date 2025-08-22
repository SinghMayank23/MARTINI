// HydroCell.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the HydroCell class which contains information on the soft background in on cell

#ifndef HydroCell_h
#define HydroCell_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Basics.h"
#include "Pythia.h"

using namespace Pythia8;

class HydroCell
{
 private:
  double itsT;
  Vec4   itsV;
  double itsQGPfrac;
  double itsTauMixed;

 public:
  HydroCell(){};//constructor
  ~HydroCell(){};//destructor

  // flow velocity
  Vec4 v() const { return itsV; };
  void v(Vec4 value) { itsV=value; };

  double vx() const { return itsV.px(); };
  double vy() const { return itsV.py(); };
  double vz() const { return itsV.pz(); };
  void vx( double value ) { itsV.px(value); };
  void vy( double value ) { itsV.py(value); };
  void vz( double value ) { itsV.pz(value); };

  // temperature
  double T() const { return itsT; };
  void T(double value) { itsT=value; };

  // QGP fraction
  double QGPfrac() const { return itsQGPfrac; };
  void QGPfrac(double value) { itsQGPfrac = value; };
  double tauMixed() const { return itsTauMixed; };
  void tauMixed(double value) { itsTauMixed = value; };
};

#endif
