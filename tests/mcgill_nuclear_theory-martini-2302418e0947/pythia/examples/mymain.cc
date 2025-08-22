#include <stdlib.h>
#include "Pythia.h"

using namespace Pythia8; 
int main() {
  Pythia pythia;
  pythia.readString("HardQCD:all = on");    
  pythia.readString("SigmaProcess:renormMultFac = 1."); 
  pythia.readString("SigmaProcess:factorMultFac = 6."); 
  pythia.readString("PhaseSpace:pTHatMin = 0."); 
  pythia.readString("PDF:useLHAPDF = on"); 
  pythia.readString("PDF:NuclearEffects = 1"); 
  pythia.readString("PDF:LHAPDFset = cteq5l.LHgrid"); 
  pythia.init( 2212, 2212, 200.);
  
  int runs = 40000;
  int workingruns = 0;
  int id;
  int bins = 10;
  double scale = 8.;

  Hist mult("pi0_pt", bins, 0., scale);

  double pt;
  double pl;
  double En;
  double y; // rapidity
  
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < runs; ++iEvent) {
    if (!pythia.next()) continue;
    workingruns+=1;
    for (int p_i = 0; p_i < pythia.event.size(); ++p_i) 
      {
	id = pythia.event[p_i].id();
	if ( id == 111)
	  {		
	    pt = pythia.event[p_i].pT(); //p_trans
	    pl = pythia.event[p_i].pz(); //p_long
	    En = pythia.event[p_i].e(); //energy
	    y = 0.5*log((En+pl)/(En-pl)); //rapidity
	    if(abs(y)<=1.)
	      {
		mult.fill(pt);
	      }
	  }
      }
  }  // End of event loop.
  pythia.statistics();
  mult/=workingruns;
  mult.table();
  return 0;
}
