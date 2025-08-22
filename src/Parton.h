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
        Vec4 itsP;           // momentum
        Vec4 itsInitialP;    // Mayank: Initial momentum
        Vec4 itsEndOfPreeqP; // Mayank: Momentum at end of Preeq
        Vec4 itsAtSplitP;    // Sangyong: momentum at the split
        Vec4 itsOldP;        // momentum during the previous timestep
        
        double itsX;       // X position
        double itsY;       // Y position
        double itsZ;       // Z position 
        double itstFinal;  // final t for fragmentation
	double itsprevtau; //proper-time for previous tau step
        
        double itsCreationTau; // the proper time tau at which the photon was emitted
                               // RMY Nov. 4th 2021
        
        int itsID;         // gluon 21, quark < 4
        int itsCol;        // color
        int itsACol;       // anti-color
        int itsElasticCollisions; // counts a parton's number of elastic collisions
        double itsMass;    // mass
        int itsFrozen;     // indicates whether parton is not evolving anymore
        int itsFormed;     // indicates whether parton is formed or not
        int itsInHydro;    // indicates whether parton is formed or not
        int itsStatus;     // status from PYTHIA
        int itsSource;     // ONLY FOR PHOTONS: Origin of the photon: 0/1/2 Initial/conversion/AMY
                           // RY: extended for QCD partons, encoding what created them.
                           // 11 emitted gluon from fermion
                           // 12 emitted gluon from gluon
                           // 13 fermions from gluon splitting 
                           // 14 gluon from q-->g conversion
                           // 15 quark or anti-quark from g-->q conversion 
                         
        int itsisCharged;
        int itsRecoil;     // if it is a recoil particle or not
  int itsAntiI;        // The position on the list of the heavy anti-quark created in the same pQCD event as this parton.
        
  double itsOldX;       // X position one timestep before
  double itsOldY;       // Y position one timestep before
  double itsOldZ;       // Z position one timestep before
        //Sangyong's additions start
        int itsDaughter; // daughter's current position in the plist. initialize to -1 (means nobody's daughter)
        int itsMother; // mother's current position in the plist. initialize to -1 (means nobody's mother)
        //Sangyong's additions end
        // 8 public functions to access these values are added at the end of the class def.
        
        // Rouz's additions start
        int num_hard_radiations;     // number of hard radiation, defined as \Delta p > pCut in MARTINI
        int num_soft_radiations;     // number of soft radiations, defined as \Delta p < pCut in MARTINI
        int num_soft_elastic; // number of soft elastic scatterings 
        int num_hard_elastic; // number of hard elastic scatterings
        // Rouz's additions end

    public:
        // momentum
        Vec4 p() const { return itsP; }
        void p (Vec4 value) { itsP=value; }
        void p (double px, double py, double pz) { itsP.px(px);  itsP.py(py); itsP.pz(pz); itsP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
        Vec4 Initialp() const { return itsInitialP; }
        void Initialp (Vec4 value) { itsInitialP=value; }
        void Initialp (double px, double py, double pz) { itsInitialP.px(px);  itsInitialP.py(py); itsInitialP.pz(pz); itsInitialP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
        Vec4 EndOfPreeqp() const { return itsEndOfPreeqP; }
        void EndOfPreeqp (Vec4 value) { itsEndOfPreeqP=value; }
        void EndOfPreeqp (double px, double py, double pz) { itsEndOfPreeqP.px(px);  itsEndOfPreeqP.py(py); itsEndOfPreeqP.pz(pz); itsEndOfPreeqP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
  Vec4 pOld() const { return itsOldP; }
  void pOld (Vec4 value) { itsOldP=value; }
  void pOld (double px, double py, double pz) { itsOldP.px(px);  itsOldP.py(py); itsOldP.pz(pz); itsOldP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
         
        //Sangyong's additions start
        Vec4 p_at_split() const { return itsAtSplitP; }
        void p_at_split (Vec4 value) { itsAtSplitP=value; }
        void p_at_split (double px, double py, double pz) 
        { itsAtSplitP.px(px);  itsAtSplitP.py(py); itsAtSplitP.pz(pz); itsAtSplitP.e(sqrt(px*px+py*py+pz*pz+itsMass*itsMass)); }
        //Sangyong's additions end
        
        // position
        double x() const { return itsX; }
        double y() const { return itsY; }
        double z() const { return itsZ; }
        double tFinal() const { return itstFinal; }
	double prevtau() const { return itsprevtau; }
        double tauAtEmission() const {return itsCreationTau;}
        void x (double value) { itsX=value; }
        void y (double value) { itsY=value; }
        void z (double value) { itsZ=value; }
        void tFinal (double value) { itstFinal=value; }
        void prevtau (double value) { itsprevtau=value; }
        void tauAtEmission(double value) { itsCreationTau=value; }
        
  double xOld() const { return itsOldX; }
  double yOld() const { return itsOldY; }
  double zOld() const { return itsOldZ; }
  void xOld (double value) { itsOldX=value; }
  void yOld (double value) { itsOldY=value; }
  void zOld (double value) { itsOldZ=value; }

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
        
  // antiI
  int antiI() const {return itsAntiI; }
  void antiI(int value) { itsAntiI=value; }
  
        // frozen
        int frozen() const { return itsFrozen; }
        void frozen(int value) { itsFrozen=value; }
        
        // formed
        int formed() const { return itsFormed; }
        void formed(int value) { itsFormed=value; }
        
        // inhydro
        int inhydro() const { return itsInHydro; }
        void inhydro(int value) { itsInHydro=value; }
        
        // status
        int status() const { return itsStatus; }
        void status(int value) { itsStatus=value; }
        
        // source
        int source() const {return itsSource; }
        void source(int value) { itsSource=value; }
        
        // isCharged // only for thermal hadrons from Cooper-Frye
        int isCharged() const {return itsisCharged; }
        void isCharged(int value) { itsisCharged=value; }
        
        // recoil
        int recoil() const {return itsRecoil; }
        void recoil(int value) { itsRecoil=value; }
          
        //Sangyong's additions start
        // mother's plist position
        int mother() const {return itsMother; }
        void daughter_of(int value) { itsMother=value; }
          
        // daughter's plist position
        int daughter() const {return itsDaughter; }
        void mother_of(int value) { itsDaughter=value; }
        //Sangyong's additions end

        //Rouz's additions begin
        int number_soft_radiation() { return num_soft_radiations;}
        void increment_soft_radiation() {num_soft_radiations += 1;}

        int number_hard_radiation() { return num_hard_radiations;}
        void increment_hard_radiation() { num_hard_radiations += 1;}

        int number_soft_elastic() { return num_soft_elastic;}
        void increment_soft_elastic() { num_soft_elastic += 1;}

        int number_hard_elastic() { return num_hard_elastic;}
        void increment_hard_elastic() { num_hard_elastic += 1;}

        void init_counts() {num_soft_elastic = 0; num_hard_elastic = 0; num_soft_radiations = 0; num_hard_radiations = 0;}
        // Rouz's additions end
};

class Source
{
  private:
    int itsType;
    double itsTau, itsX, itsY, itsEta;
    double itsDE, itsDpx, itsDpy, itsDpz;
    double itsVx, itsVy, itsVz;

  public:
    int type() const {return itsType;}
    void type(int value) {itsType=value;}

    double tau() const {return itsTau;}
    void tau(double value) {itsTau=value;}

    double x() const {return itsX;}
    void x(double value) {itsX=value;}

    double y() const {return itsY;}
    void y(double value) {itsY=value;}

    double eta() const {return itsEta;}
    void eta(double value) {itsEta=value;}

    double dE() const {return itsDE;}
    void dE(double value) {itsDE=value;}

    double dpx() const {return itsDpx;}
    void dpx(double value) {itsDpx=value;}

    double dpy() const {return itsDpy;}
    void dpy(double value) {itsDpy=value;}

    double dpz() const {return itsDpz;}
    void dpz(double value) {itsDpz=value;}

    double vx() const {return itsVx;}
    void vx(double value) {itsVx=value;}

    double vy() const {return itsVy;}
    void vy(double value) {itsVy=value;}

    double vz() const {return itsVz;}
    void vz(double value) {itsVz=value;}
};

#endif
