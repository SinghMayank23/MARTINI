#include "Pythia.h"
#define nbins 500
#define nEvent 20
#define ptmax 100.
#define nSubrun 5

using namespace Pythia8; 

int main() {

  Pythia pythia;

  int ievent;
  bool n1,n2;
  bool flag;

  double pion_hist[nbins];					// The total histogram for pi-
  double pion_histS1[nbins];


  double jetpTminIn;
  double sigmaInelastic;
  double sigmaRatio;
  double sigmaRatiotot=0.;

  for(int ih=0; ih<nbins; ih++){		//Initialize the histograms to be zero everywhere:
    pion_histS1[ih] = 0.;
    pion_hist[ih] = 0.;
  }

  // Loop over number of subruns.
  for (int iSubrun = 0; iSubrun < nSubrun; ++iSubrun) {

    pythia.readFile("setup_main_pp_LHC_pi-_weighted.dat", iSubrun);
    jetpTminIn=pythia.parm("PhaseSpace:pTHatMin");
    pythia.setJetpTmin(jetpTminIn);
    pythia.init();
    pythia.settings.listChanged();

    sigmaInelastic = pythia.totalCrossSection() - pythia.elasticCrossSection();
    if(!iSubrun==0){
      sigmaRatio = pythia.jetCrossSection() / sigmaInelastic;
      sigmaRatiotot += sigmaRatio;
//      cout << "inelastic = " << scientific << setprecision(2) << sigmaInelastic << endl;
//      cout << "ratio_jet/inelastic = " << scientific << setprecision(2) << sigmaRatio << endl;
    }

    n1=false;						//if neutron, n1=true
    n2=false;

    //creating a file for events:
//    stringstream ss;
//    ss << iSubrun;
//    string subrun = ss.str();
//    string e_name = /raid/erebos/radius/shaitan/chanwook/MARTINI/main/EventOutput/ALICE_pp_" + subrun + ".dat"; 
//    fstream foutEvents( e_name.c_str() ,ios::out);   

    // The main loop for pythia event
    ievent=0;
    while(ievent<nEvent){
      flag=pythia.next(n1,n2);
      if(flag==true){
        ievent+=1;

    for (int i=0; i<pythia.event.size(); ++i){
      //Output the event
//  	  if(pythia.event[i].isFinal()){
//	      foutEvents << pythia.event[i].id() << " " << pythia.event[i].pT() << " " 
//		       << pythia.event[i].y() << " " << pythia.event[i].eta() << endl;
//	    }

      if (pythia.event[i].isFinal() && pythia.event[i].id() == -211 ){
        int pbin = (int)((double)nbins*pythia.event[i].pT()/ptmax);
        if(pbin >= 0 && pbin < nbins){
          if(!iSubrun==0)
            pion_hist[pbin] += sigmaRatio*(double)nbins/(double)nEvent;
          else
            pion_histS1[pbin] += (double)nbins/(double)nEvent;
        }
      }
    }
      }

    }
		//End of the main loop

//  	foutEvents << "EndOfEvent" << endl;

//    if(iSubrun==0)
//      foutEvents << "inelXSec = " << sigmaInelastic << endl;
//    else
//      foutEvents << "jetXSec[" << iSubrun << "] = " << pythia.jetCrossSection() << endl;
 
//    foutEvents.close();
  }
  //End of loop for Subrun

  cout << "sigmaRatiotot = " << scientific << setprecision(2) << sigmaRatiotot << endl;

  for(int it=0; it<nbins; it++){
    pion_hist[it] += pion_histS1[it]*(1.- sigmaRatiotot);
  }

  string weight_name = "pp_LHC_pi-.dat";
  fstream foutweight( weight_name.c_str() ,ios::out);   

  for(int ib=0; ib<nbins; ib++){
    foutweight << ptmax*(double)ib/(double)nbins << " " << pion_hist[ib]<< endl;
  }

  foutweight.close();

  return 0;
}  
