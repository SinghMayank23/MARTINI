#ifndef Constants_h //avoids multiple inclusions of the header file
#define Constants_h
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>

using namespace std;

/// saves all the important constants to distribute them to the classes that need them
class Constants
{
private:
  double dt;
  double alphas;
  double dE;
  int Nf;
  int strate;
  double initialE;
  double width;
  double maxTime;
  double maxEnergy;
  int BetheHeitler;
  int BDMPS;
  int doRad;
  int doCol;
  int newRadGamma;
  int mperfermi;

 public:
  static const double PI = 3.141592653589793;
  static const double hbarc = 0.197327053;
  static const double EulerGamma = 0.57721566490153286061l;
  /// cb = -EulerGamma + Zeta'(2)/Zeta(2) + Log(2)
  static const double cb = -1.147176657996066; 
  /// cf = -EulerGamma + Zeta'(2)/Zeta(2)
  static const double cf = -0.4540294774361204;
  static const double cs = -1.66246;
  
  static const int Nc = 3;

/// now these are the parameters as in the tabulated data files of collisional transition rates
  static const double LogEmax=5;//11.5; //5;
  static const double LogEmin=1; //0.91;
  static const double LogStepE=0.2; //0.3;
  static const double LogOmegaMin=-5; //-5;
  static const double LogStepOmega=0.2; //0.3;
  /// Tfile is the temperature at which the transition rate in the data file was computed (didnt use 1, to reduce the error from different matching points between hard and soft scale) 
  static const double Tfile=0.4;
  static const double stepT=1;
};
#endif
