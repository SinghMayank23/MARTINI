#ifndef PauseEvolution_H
#define PauseEvolution_H

#include <stdio.h>
#include <string>
#include "Pythia8/Pythia.h"

//**************************************************************************

// PauseEvolution is a derived class for user access to program execution.
// It is a simple example, to kill events with > nMax partons at scale 
// pTcheck in the combined evolution, but only counting partons in the 
// hardest interaction and its associated ISR + FSR activity.

using namespace Pythia8;

class PauseEvolution : public UserHooks {

public:

  // Constructor.
  PauseEvolution( int nMaxIn, double pTIn ) : nMax(nMaxIn), 
    pT(pTIn){}

  // Possibility to veto combined MI + ISR + FSR evolution and
  // kill event, e.g. for MLM-style matching of matrix elements.
  virtual bool canVetoPT() {return true;}  

  // Transverse-momentum scale for veto test. 
  virtual double scaleVetoPT() {return pT;} 

  // Decide whether to veto current event or not.
  virtual bool doVetoPT( int , const Event& ) {return false;}  

  void setPT(double value){pT=value;}
  double getPT(){return pT;}

private:

  // Saved values from constructor.
  int    nMax;
  double pT;

};

//**************************************************************************

#endif
