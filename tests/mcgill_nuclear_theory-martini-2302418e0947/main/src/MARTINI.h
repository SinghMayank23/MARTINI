// MARTINI.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class
// MARTINI: provide the main user interface to everything else.

#ifndef MARTINI_H
#define MARTINI_H

#include <stdio.h>
#include <string>
#include <cmath>
#include <list>
#include "Pythia8/Pythia.h"            // include Pythia
#include "LHAPDF/LHAPDF.h"     // include Les Houches Accord for PDFs
#include "Basics.h"
#include "Setup.h"
#include "Import.h"
#include "ImportLRates.h"
#include "Random.h"
#include "Parton.h"
#include "Rates.h"
#include "Elastic.h"
#include "HydroCell.h"
#include "HydroSetup.h"
#include "Testing.h"
#include "Information.h"
#include "Glauber.h"
#include "Welcome.h"
#include "PauseEvolution.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>

//using namespace Pythia8;

class MARTINI
{
 private:
  Import  *import;
  ImportLRates  *importLRates;
  Random  *random;
  Testing *testing;
  Rates   *rates;
  Elastic *elastic;
  HydroSetup *hydroSetup;
  Glauber *glauber;
  Information info;
  vector<HydroCell> * lattice;     

  vector<ReturnValue> nucleusA;  // list of x and y coordinates of nucleons in nucleus A      
  vector<ReturnValue> nucleusB;  // list of x and y coordinates of nucleons in nucleus B 

  gsl_rng *gsl_rand ;

  int runs;                // number of events
  int evolution;           // do evolution or not
  double cmEnergy;         // center of mass energy
  double maxTime;          // maximal time in fm
  double fixedT;           // temperature in GeV
  double dtfm;             // time step in fm
  double alpha_s;          // strong coupling constant
  int fixedEnergy;         // switch between phythia generated initial partons (0) or fixed energy initial partons (1)
  int fixedTemperature;    // switch between fixed temperature and hydro evolution
  int fragmentationSwitch; // do fragmentation or not
  int photonSwitch;        // produce photons or not
  int Nf;                  // number of colors
  double pCut;             // low momentum cut for radiated partons (below this they are not taken into the list) in [GeV], fixed value
  bool pCutPropToT;      // if on low momentum cut for radiated partons will be pCut=4T, regardless of parameter pCut 
  double jetpTmin;         // the lower momentum cutoff on the jet cross section -> sent to PYTHIA
  double jetpTmax;         // the upper momentum cutoff on the jet cross section -> sent to PYTHIA -CP 14/1/2015
  double jetXSec;          // the jet cross section from PYTHIA computed with jetpTmin as infrared cutoff
  double totalXSec;        // the total cross section from PYTHIA (including elastic and diffractive parts)
  double elasticXSec;      // the elastic cross section from PYTHIA
  double inelasticXSec;    // the inelastic cross section (totalXSec-elasticXSec)
  int fullEvent;           // 0: just one pp event per total event (position from Glauber), 1: sample number of jet events from Glauber
  int moveBeforeTau0;      // 0: do not move partons before tau_0, 1: evolve position before tau_0
  int doElastic;           // 0: no elastic collisions, 1: include elastic collisions
  int doRadiative;         // 0: no radiation processes, 1: include radiation
  int transferTransverseMomentum; 
                           // 0: Eikonal approximation, 1: include momentum broadening
  int nuclearEffects;      // 0: no nuclear effects, 1: nuclear effects EKS98
  int NCol;
  double qtTot;
  double lengthTraveled;   // store the length traveled
  double mfp;              // store mean free path (used for testing)

  int tauEtaCoordinates;
  double hydroTau0;        // tau_0 in the hydro data files
  double hydroTauMax;      // tau_max in the hydro data files
  double hydroDtau;        // step dtau in fm/c in the hydro data files
  double hydroXmax;        // maximum x in fm in the hydro data files [-xmax, +xmax] for both x and y
  double hydroZmax;        // maximum z in fm in the hydro data files [-zmax, +zmax] for 3D hydro
  double hydroDx;          // step dx in fm in the hydro data files
  double hydroDz;          // step dz in fm in the hydro data files in the z-direction for 3D hydro
  double hydroTfinal;      // temperature at which evolution of the jet stops
  int hydroWhichHydro;     // choose a hydro evolution model to use
  int hydroSubset;         // choose a subset of parameters for a given evolution model (currently only for Jeon/Schenke hydro)
  bool hydroViscous;       // use viscous hydro or ideal if supported
  bool fullVacuumShower;   // decide whether to do a full vacuum shower (on) or to cut off at Q=\sqrt{E/\tau_0} (off)
  bool fixedDistribution;  // decide whether the hydro background should be fixed or evolving
  int Ldependence;         // decide whether to use the original AMY (0), my parametrized L dependence (1), or ...
  int fragmentationMethod; // determines the fragmentation method: 1=strings only among jet partons as in pp
                           //                                      2=conenct strings in the end including near-by medium partons
  int rateSelector;
  vector<double> tauMixed; // time when mixed phase starts (for Kolb hydro only)
  

  double glauberSigmaNN;   // total cross section
  string glauberTarget;    // name of target nucleus
  string glauberProjectile;// name of projectile nucleus
  double glauberImpactParam; // impact parameter
  int glauberIMax;         // max number of interpolation points
  bool glauberEnvelope;    // turn additional phenomenological envelope over x-y distribution on or off
  string dataPath;         // path to main directory 
  
  bool useLHAPDF;          // switch that holds whether or not to use LHAPDF (Les Houches Accord Parton Distribution Functions)
  string PDFname;          // name of the PDF to use
  int PDFmember;           // member function of the PDF set to use
  bool trackPartons;       // switch to turn on tracking of partons for pretty pictures
  bool trackHistory;       // switch to turn on tracking of partons trajectory and energy-momentum loss
  double initialXjet;      // initial parameters when putting a jet by hand (for setting FixedEnergy).
  double initialYjet;
  double initialPXjet;
  double initialPYjet;
  double multElasticRate;  // simple multiplicative factor to increase the elastic rate
  int totalNNs;            // total number of hard collisions in a full event
  bool allFromCenter;      // if on all partons are produced at x=y=0
  int quarksEqualGluons;   // 0: normal, 1: rates involving quarks are increased by a factor of 9/4, 2: rates involving gluons are decreased by 4/9
  bool runningCoupling;   // If on, coupling constant for hard process is running 

  //These settings were added for handling the coordinates of the collisions from MUSIC. -CFY 11/2/2010
  int nbinFromFile;
  int file_number;
  string evolution_name;
  //The array of the coordinates of the collisions.
  double **binary;
  int Nbin;
  
  //These settings were added for runs examining heavy quarks. -CFY 5/24/2011
  int examineHQ;
  int setInitialPToZero;
  double totalHQXSec;
  double charmWidth;
  double bottomWidth;
  double T_C_HQ;
  double TwoPiTD_HQ;

 public:
  MARTINI();
  ~MARTINI();
  /// declare Welcome object to display the welcome screen
  Welcome welcome;
  /// declare Pythia object as public - will be accessible by main program
  Pythia pythia;
  /// Make the setup class public for ease of use.
  Setup settings;
  /// Initialize MARTINI 
  bool init(int argc, char** argv);
  /// ReadString as adopted from PYTHIA will read in values for parm, flag, word, and mode
  bool readString(string line, bool warn=true);
  /// ReadFile as adopted from PYTHIA will read in values for parm, flag, word, and mode from a specified file
  bool readFile(string, bool warn = true);
  /// get T and the flow velocity at the current position usig interpolation in 3D resp. 4D (2+1D and 3+1D hydro)
  HydroInfo getHydroValues(double x, double y, double z, double t);
  /// for now this is where the whole evolution is done - this will change 
  void pythiaEvents();
  /// calls the next MARTINI event. Only good for events with initial parton content from PYTHIA
  vector<Parton>* next();
  /// Samples the momentum fraction z for Peterson fragmentation of heavy quarks:
  double SamplePetersonFunction(double M);
  // Checks whether or not a heavy quark pair would be bound in the color evaporation model:
  int isBound(int id1, int id2, double E_CM);
  /// routine for preparing and calling the fragmentation done by PYTHIA in a full event
  int fragmentation(vector<Parton> ** plist, int currentEvent = 1);
  // Hadronizes a single heavy quark into an open heavy flavor meson:
  void hadronizeOneHeavyQuark(int id, double *pmu, double M);
  /// This hadronizes heavy quarks in a full event according to a simple thermal model by C. Young.
  int hadronizeHeavyQuarks(vector<Parton> ** plist);
  void hadronizeByColorEvaporation(int id1, int id2, double E_CM, double *p1mu, double *p2mu, double M1, double M2, int status);  
  /// This simply removes quarkonium from the PYTHIA event and then decays the remaining open heavy flavor hadrons:
  int decayOpenHeavyFlavor();
  /// generate an event sampling the T_A's directly including a shower with a minimal p_T from PYTHIA 
  int generateEvent(vector<Parton> *plist);
  int generateEventHeavyQuarks(vector<Parton> *plist);
  int generateTestEvent_HQ(vector<Parton> *plist, double L);
  /// sets basic PYTHIA variables and calls its initialization
  void initPythia();
  /// evolve the system for one time step 
  int evolve(vector<Parton> **plist, int counter, int it);
  int evolveHeavyQuarks(vector<Parton>  **plist, int counter, int it);
  int evolveAndHadronizeHeavyQuarks(vector<Parton>  **plist, int counter, int it);
  double CornellPotential(double r);
  double CornellPotential(double r, int id1, int id2);
  double F(double r, double T);
  double FCornell(double r);
  double FCornell(double r, int id1, int id2);
  void stepMomentumInCMFrame(double *ptot, double *pmu, double dt, double *deltax, double *pnew, double T);
  double rInCMFrame(double *p1, double *p2, double *deltaX);
  /// functions for testing: these call the functions defined in the class Testing
  void generateSpectrum(double p, int process){testing->generateSpectrum(p, fixedT, alpha_s, Nf, 
									 random, import, rates, process);};
  void sampleElasticRate(double p, int process){testing->sampleElasticRate(p, fixedT, alpha_s, Nf, 
									   random, import, elastic, process);};
  void sampleElasticRateOmegaQ(double p, double omega, int process){testing->sampleElasticRateOmegaQ(p, omega, fixedT, alpha_s, Nf, 
										 random, import, elastic, process);};
  void quarkBrick(double p){testing->quarkBrick(p, fixedT, alpha_s, Nf, dtfm, maxTime, random, import, rates, runs);};
  void thermalPlot(){testing->thermalPlot(fixedT, random);};
  void sampleTA();
  void plotHydroData();
  void plotHydroDataXYZ();
  int returnFixedEnergy(){return fixedEnergy;};
  int returnTotalNNs(){return totalNNs;};
  double returnMaxTime(){return maxTime;};
  double returnDtfm(){return dtfm;};
  int returnRuns(){return runs;};
  double returnQtTot(){return qtTot;};
  double returnMfp(){return mfp;};
  void setQtTot(double value){qtTot = value;};
  int returnNCol(){return NCol;};
  void setNCol(int value){NCol = value;};
  int returnFullEvent(){return fullEvent;};
  int returnEvolution(){return evolution;};
  int returnFragmentationSwitch(){return fragmentationSwitch;};
  bool returnTrackHistory(){return trackHistory;};
  double returnInitialXjet(){return initialXjet;};
  double returnInitialYjet(){return initialYjet;};
  double returnInitialPXjet(){return initialPXjet;};
  double returnInitialPYjet(){return initialPYjet;};
  double returnHydroTau0(){return hydroTau0;};

  //Sangyong's additions
  int hasDaughter(vector<Parton> *plist, int i);
  int hasMother(vector<Parton> *plist, int i);

  // history
  vector<double> *Tt; // stores the history of the trajectory: t, x(t), y(t), z(t), dE/dt(t), dpx/dt(t), dpy/dt(t), dpz/dt(t) *for one particle*
  vector<double> *Tx;
  vector<double> *Ty;
  vector<double> *Tz;
  vector<double> *TE;
  vector<double> *Tpx;
  vector<double> *Tpy;
  vector<double> *TQGPfrac;
  vector<double> *TdEdt; 
  vector<double> *Tdpxdt; 
  vector<double> *Tdpydt; 
  vector<double> *Tdpzdt; 
};

#endif
