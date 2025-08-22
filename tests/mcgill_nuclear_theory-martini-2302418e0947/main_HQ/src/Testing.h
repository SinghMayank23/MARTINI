// Testing.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains several testing routines

#ifndef Testing_h
#define Testing_h

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include "Basics.h"
#include "Random.h"
#include "Import.h"
#include "Rates.h"
#include "Elastic.h"
#include "Pythia.h"
#include "Glauber.h"


using namespace Pythia8;

class Testing
{
 private:

 public:
  Testing(){};//constructor
  ~Testing(){};//destructor
  void generateSpectrum(double p, double T, double alpha_s,
                        int Nf, Random *random, Import *import, Rates * rates, int process);
  void sampleElasticRate(double p, double T, double alpha_s,
                        int Nf, Random *random, Import *import, Elastic * elastic, int process);
  void sampleElasticRateOmegaQ(double p, double omega, double T, double alpha_s, int Nf,
			       Random *random, Import *import, Elastic * elastic, int process);
  void quarkBrick(double initp, double T, double alpha_s, 
		  int Nf, double dtfm, double maxTime, 
		  Random * random, Import * import, Rates * rates, int runs);
  void thermalPlot(double T, Random * random);
  void initialDist(Random *random, Glauber *glauber);
  void sampleTA(Random *random, Glauber *glauber);
};

#endif
