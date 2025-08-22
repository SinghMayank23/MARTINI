#include "Pythia.h"
#define nbins 100
#define nEvent 1
#define ptmax 10.

using namespace Pythia8; 

int main() {

  int ievent=0;
  bool n1,n2;
  bool flag;
  Pythia pythia;

  pythia.readFile("setup_main_pp_LHC_pi-.dat");
  pythia.init();

  double pion_hist[nbins];          // The histogram for pi-

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    pion_hist[ih] = 0.;
  }

  n1=false;						//if neutron, n1=true
  n2=false;

  //Creating a file for event output
  string e_name = "./EventOutput/ALICE_pp.dat"; 
  fstream foutEvents( e_name.c_str() ,ios::out); 
  
  // The main loop
  while(ievent<nEvent){
    flag=pythia.next(n1,n2);
    if(flag==true){
      ievent+=1;
  
  for (int i=0; i<pythia.event.size(); ++i){
  	if(pythia.event[i].isFinal()){
	    foutEvents << pythia.event[i].id() << " " << pythia.event[i].pT() << " " 
		     << pythia.event[i].y() << " " << pythia.event[i].eta() << endl;
	  }
    if (pythia.event[i].isFinal() && pythia.event[i].id() == -211 ){
      int pbin = (int)((double)nbins*pythia.event[i].pT()/ptmax);


      if(pbin >= 0 && pbin < nbins){
        pion_hist[pbin] += (double)nbins/(double)nEvent;
      }
    }
  }
    }

  }
  foutEvents << "EndOfEvent" << endl;

  string charged_name = "pp_LHC_pi-_test.dat";
  fstream foutcharged( charged_name.c_str() ,ios::out);   

  for(int ib=0; ib<nbins; ib++){
    foutcharged << ptmax*(double)ib/(double)nbins << " " << pion_hist[ib]<< endl;
  }

  foutcharged.close();
  foutEvents.close();

  return 0;
}  
