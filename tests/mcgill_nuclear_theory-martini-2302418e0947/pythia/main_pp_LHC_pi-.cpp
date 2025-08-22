#include "Pythia.h"
#define ET_min 0.
#define etaMax 1.0
#define nbins 500
#define nEvent 1
#define ptmax 100.
#define Subrun 0

using namespace Pythia8; 

int main() {

  Pythia pythia;
  pythia.readFile("setup_main_pp_LHC_pi-.dat");

  if (Subrun==0) {
    pythia.readString("PhaseSpace:pTHatMin = 10.");
    pythia.readString("PhaseSpace:pTHatMax = 30.");
  }
  else if (Subrun==1) {
    pythia.readString("PhaseSpace:pTHatMin = 30.");
    pythia.readString("PhaseSpace:pTHatMax = 50.");
  }
  else if (Subrun==2) {
    pythia.readString("PhaseSpace:pTHatMin = 50.");
    pythia.readString("PhaseSpace:pTHatMax = 70.");
  }
  else if (Subrun==3) {
    pythia.readString("PhaseSpace:pTHatMin = 70.");
    pythia.readString("PhaseSpace:pTHatMax = 90.");
  }
  else if (Subrun==4) {
    pythia.readString("PhaseSpace:pTHatMin = 90.");
    pythia.readString("PhaseSpace:pTHatMax = 0.");
  }

  double jetpTminIn;
  jetpTminIn=pythia.parm("PhaseSpace:pTHatMin");
  pythia.setJetpTmin(jetpTminIn);

  pythia.init();
  pythia.settings.listChanged();

  int ievent=0;
  bool n1,n2;
  bool flag;

  double sigmaInelastic;
  double sigmaRatio;

  double bin[nbins];
  double pion_hist[nbins];          // The histogram for pi-

//  if(!Subrun==0){
    sigmaInelastic = pythia.totalCrossSection() - pythia.elasticCrossSection();
    sigmaRatio = pythia.jetCrossSection() / sigmaInelastic;
    cout << "jet = " << scientific << setprecision(5)  << pythia.jetCrossSection() 
         << " inelastic = " << sigmaInelastic << endl;
    cout << "ratio_jet/inelastic = " << scientific << setprecision(5) << sigmaRatio << endl;
//  }

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    bin[ih] = ptmax*((double)ih+0.5)/(double)nbins;
    pion_hist[ih] = 0.;
  }

  n1=false;						//if neutron, n1=true
  n2=false;

  //Creating a file for event output

  // The main loop
  while(ievent<nEvent){
    flag=pythia.next(n1,n2);
    if(flag==true){
      ievent+=1;
  
  for (int i=0; i<pythia.event.size(); ++i){
    if (pythia.event[i].mT() > ET_min
				&&fabs(pythia.event[i].eta()) < etaMax
				&&pythia.event[i].isFinal() 
//				&& pythia.event[i].id() == -211){
	      && pythia.event[i].isCharged()){
      int pbin = (int)((double)nbins*pythia.event[i].pT()/ptmax);

      if(pbin >= 0 && pbin < nbins){
        pion_hist[pbin] += ((double)nbins/ptmax)/(double)nEvent;
      }
    }
  }
    }
  }
  stringstream ss;
  ss << Subrun;
  string subrun = ss.str();
/*
//  string pion_name = "pp_LHC_pi-_|eta|<1.0_" + subrun + ".dat";
  string pion_name = "pTRtest_" + subrun + ".dat";
  fstream foutpion( pion_name.c_str() ,ios::out);   

  for(int ib=0; ib<nbins; ib++){
    foutpion << bin[ib] << " " << pion_hist[ib]/(2.*M_PI*bin[ib]) << endl;
  }

  foutpion.close();
*/
  return 0;
}  
