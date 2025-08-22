// Parton.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the Parton class which contains information on one parton

#ifndef Parton_H
#define Parton_H

#include "Pythia8/Pythia.h" // include Pythia for Vec4 class

using namespace Pythia8;

class Parton
{
 private:
  Vec4 itsP;         // momentum
  Vec4 itsOldP;      // momentum during the previous timestep
  Vec4 itsInitialP;  // initial momentum
  double itsX;       // X position
  double itsY;       // Y position
  double itsZ;       // Z position 
  double itsOldX;       // X position one timestep before
  double itsOldY;       // Y position one timestep before
  double itsOldZ;       // Z position one timestep before
  double itsInitialX;// X position
  double itsInitialY;// Y position
  double itsInitialZ;// Z position 
  double itsInitialT;// initial time when parton was created 
  double itsDistanceTraveled;//distance traveled by the parton in the medium
  double itstFinal;  // final t for fragmentation
  
  int itsID;         // gluon 21, quark < 4
  int itsCol;        // color
  int itsACol;       // anti-color
  int itsSplittings; // keeps track of a parton's splittings
  int itsElasticCollisions; // counts a parton's number of elastic collisions
  double itsMass;    // mass
  int itsFrozen;     // indicates whether parton is not evolving anymore
  int itsStatus;     // status from PYTHIA
  int itsEventNumber;// in a full event - save here which individual nucleon-nucleon collision the parton came from
  int itsSource;     // ONLY FOR PHOTONS: Origin of the photon: 0/1/2 Initial/conversion/AMY
  int itsAntiI;        // The position on the list of the heavy anti-quark created in the same pQCD event as this parton.

 public:
  // momentum
  Vec4 p() const { return itsP; }
  void p (Vec4 value) { itsP=value; }
  void p (double px, double py, double pz) { itsP.px(px);  itsP.py(py); itsP.pz(pz); itsP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
  Vec4 pOld() const { return itsOldP; }
  void pOld (Vec4 value) { itsOldP=value; }
  void pOld (double px, double py, double pz) { itsOldP.px(px);  itsOldP.py(py); itsOldP.pz(pz); itsOldP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }

  Vec4 pini() const { return itsInitialP; }
  void pini (Vec4 value) { itsInitialP=value; }

  // position
  double x() const { return itsX; }
  double y() const { return itsY; }
  double z() const { return itsZ; }
  double tFinal() const { return itstFinal; }
  void x (double value) { itsX=value; }
  void y (double value) { itsY=value; }
  void z (double value) { itsZ=value; }
  void tFinal (double value) { itstFinal=value; }

  double xOld() const { return itsOldX; }
  double yOld() const { return itsOldY; }
  double zOld() const { return itsOldZ; }
  void xOld (double value) { itsOldX=value; }
  void yOld (double value) { itsOldY=value; }
  void zOld (double value) { itsOldZ=value; }

  double xini() const { return itsInitialX; }
  double yini() const { return itsInitialY; }
  double zini() const { return itsInitialZ; }
  double tini() const { return itsInitialT; }
  void xini (double value) { itsInitialX=value; }
  void yini (double value) { itsInitialY=value; }
  void zini (double value) { itsInitialZ=value; }
  void tini (double value) { itsInitialT=value; }

  double distanceTraveled() const { return itsDistanceTraveled; }
  void distanceTraveled (double value) { itsDistanceTraveled=value; }
  
  // id
  int  id() const { return itsID; }
  void id(int value) { itsID=value; }

  // mass
  double  mass() const { return itsMass; }
  void mass(double value) { itsMass=value; }
 
  // color
  int col() const { return itsCol; }
  void col(int value) { itsCol=value; }

  // anti-color
  int acol() const { return itsACol; }
  void acol(int value) { itsACol=value; }

  // splittings
  int splits() const { return itsSplittings; }
  void splits(int value) { itsSplittings=value; }

  // elastic collisions
  int elasticCollisions () const { return itsElasticCollisions; }
  void elasticCollisions(int value) { itsElasticCollisions=value; }

  // frozen
  int frozen() const { return itsFrozen; }
  void frozen(int value) { itsFrozen=value; }

  // status
  int status() const { return itsStatus; }
  void status(int value) { itsStatus=value; }

  // eventNumber
  int eventNumber() const { return itsEventNumber; }
  void eventNumber(int value) { itsEventNumber=value; }

  //  vector<double> tempHist; // temperature history

  // source
  int source() const {return itsSource; }
  void source(int value) { itsSource=value; }

  // antiI
  int antiI() const {return itsAntiI; }
  void antiI(int value) { itsAntiI=value; }

};

#endif
