#ifndef Setup_h //avoids multiple inclusions of the header file
#define Setup_h
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <cstdlib>

#include "tools/ReadParameter.h"
#include "Control.h"
#include "Constants.h"

using namespace std;

/// reads in parameters and sets up output files
class Setup
{
private:
  int setupsSize;
  int setupsRank;
  ifstream ifstr;
  string inputFile;
  string tempFileName[20];
  FileNames* fileNameList;
  ReadParameter* readparam;
  double T;
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
  Constants * constants;

 public:
  Setup(int argc, char** argv, string f, const int size=1, const int rank=0);
  ~Setup();
/// generates an instance of the Control class and passes all relevant parameters to it
  Control* createControl();
/// allows access to the private variable fileNameList
  FileNames* getFileNameList();
};
#endif
