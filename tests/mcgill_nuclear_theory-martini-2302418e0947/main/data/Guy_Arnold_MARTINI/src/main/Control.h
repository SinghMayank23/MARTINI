#ifndef Control_h //avoids multiple inclusions of the header file
#define Control_h
#include <cmath>
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>

#include "tools/FileNames.h"
#include "Evolution.h"

/// controls the evolution, loops over time steps in evolve
class Control
{
 private:

  static const double hbarc = 0.197327053;

  double dt;
  double dE;
  double T;
  double alphas;
  int Nf;
  double initialE;
  double width;
  double maxTime;
  double maxEnergy;
  int maxSteps;
  int strate;
  int BetheHeitler;
  int BDMPS;
  int doRad;
  int doCol;
  int newRadGamma;
  int mperfermi;
  int size;
  int rank;

  double LogEmax;
  double LogEmin;
  double LogStepE;
  double LogOmegaMin;
  double LogStepOmega;
  double Tfile;
  double stepT;

  double* Pq;
  double* Pg;
  double* Pem;
  double* PqA;
  double* PgA;

  ///arrays with transition rates:
  /**
     E=Emin+iE*stepE
     T=Tmin+iT*stepT
     omega=-Emax+iOmega*stepOmega
     then
     dGamma/(domega dt)_qq(E,T,omega)=trqq[iE*(Tsize*(omegaSize+1))+iT*(omegaSize+1)+iOmega]
  **/

  double* trqq;
  double* trqg;
  double* trgq;
  double* trgg;

  Evolution* timeEvolution;
  const FileNames* controlsFileNameList;
 public:
  Control(const FileNames* fileList, const Constants* constants, double dt, double dE, double T, const double alphas, int Nf, const double initialE, const double width, const double maxTime, const double maxEnergy, const int strate, const int BetheHeitler, const int BDMPS, const int doRad, const int doCol, const int newRadGamma, const int mperfermi, const int size=1, const int rank=0);
  ~Control(){delete [] Pq; delete [] Pg; delete [] Pem; delete [] PqA; delete [] PgA;};
  ///loops over time steps
  void evolve();

  double* getPq(){return Pq;};
  double* getPg(){return Pg;};

};
#endif

