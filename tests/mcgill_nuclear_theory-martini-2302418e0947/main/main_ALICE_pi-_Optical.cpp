#include "MARTINI.h"
#define ET_min 0.
#define etaMax 1.0
#define numEvents 5000
#define nbins 500
#define ptmax 100.

using namespace std;

int main(int argc, char* argv[]) 
{
  MARTINI martini;
  martini.readFile("setup_main_ALICE_pi-_Optical.dat");

  string val = argv[2];
  int Subrun = atoi(val.c_str());

  if (Subrun==0) {
    martini.readString("General:JetPTMin = 0.");
    martini.readString("General:JetPTMax = 20.");
  }
  else if (Subrun==1) {
    martini.readString("General:JetPTMin = 20.");
    martini.readString("General:JetPTMax = 40.");
  }
  else if (Subrun==2) {
    martini.readString("General:JetPTMin = 40.");
    martini.readString("General:JetPTMax = 60.");
  }
  else if (Subrun==3) {
    martini.readString("General:JetPTMin = 60.");
    martini.readString("General:JetPTMax = 80.");
  }
  else if (Subrun==4) {
    martini.readString("General:JetPTMin = 80.");
    martini.readString("General:JetPTMax = 0.");
  }

  martini.init(argc, argv);
  martini.settings.listChanged();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  Event fullEvent;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  double bin[nbins];
  double pion_hist[nbins];          // The histogram for particle
  double pion_SE[nbins];

  double jetXSec,totaljetXSec;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    bin[ih] = ptmax*((double)ih+0.5)/(double)nbins;
    pion_hist[ih] = 0.;
    pion_SE[ih] = 0.;
  }

  jetXSec = martini.pythia.jetCrossSection();
  totaljetXSec = martini.pythia.totaljetCrossSection();

  //Make a tag for the output:

  string n_string = argv[1];
  n_string = "PbPb_charged_"+n_string;

  //creating a file for events
//  string e_name = "/raid/erebos/chanwook/EventOutput/ALICE_" + n_string + "_" + subrun + ".dat";
//  fstream foutEvents( e_name.c_str() ,ios::out);

  while(event_counter<numEvents)      // The main loop
    {
      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      martini.generateEvent(plist[0]); 
	  // version that samples number of collisions with Glauber model
      
      int counter = 0;
      if (martini.returnEvolution() == 1)         // evolve in medium if settings want that
	{
	  for(int i=0; i<mt; i++) // loop over all time steps 
	    {
	      counter = martini.evolve(plist, counter, i);
	      counter+=1;
	    }
	}
      // fragmentation
      fullEvent.clear();
      martini.pythia.event.clear();
      int frag_success = martini.fragmentation( plist );
      if(frag_success){
	event_counter++;
	fullEvent = martini.pythia.event;

	//Finally, bin the minus pions:
	int eventSize = fullEvent.size();
	for(int ip=0; ip<eventSize; ip++){
	  if(fullEvent[ip].mT() > ET_min
	     && fabs(fullEvent[ip].eta()) < etaMax
	     && fullEvent[ip].isFinal() 
//	     && fullEvent[ip].id() == -211){
	     && fullEvent[ip].isCharged()){
	    int pbin = (int)((double)nbins*fullEvent[ip].pT()/ptmax);
	    if(pbin >= 0 && pbin < nbins){
	      pion_hist[pbin] += ((double)nbins/ptmax)/(double)numEvents;
        pion_SE[pbin] += 1./(double)numEvents;
	    }
	  }
	  if(fullEvent[ip].isFinal()){
//	    foutEvents << fullEvent[ip].m() << " " << fullEvent[ip].charge() << " " << fullEvent[ip].pT() << " " 
//		       << fullEvent[ip].phi() << " " << fullEvent[ip].eta() << endl;
	  }
	}
//	foutEvents << "EndOfEvent" << endl;
    }

  string pion_name = n_string+"_RS8_Edep1_"+argv[2]+".dat";
  fstream foutpion( pion_name.c_str() ,ios::out);

  for(int ib=0; ib<nbins; ib++){
    pion_hist[ib] /= (2.*M_PI*bin[ib]);
    if (pion_SE[ib] < 1.0) pion_SE[ib] = sqrt((pion_SE[ib]-pion_SE[ib]*pion_SE[ib])/(double)numEvents)/(2.*M_PI*bin[ib]*(ptmax/(double)nbins));
    else pion_SE[ib] = 0.;
    foutpion << bin[ib] << " " << pion_hist[ib] << " " << pion_hist[ib]-pion_SE[ib] << " " << pion_hist[ib]+pion_SE[ib] << endl;
  }
  foutpion.close();
//  foutEvents.close();

  cout << "\nnumEvents = " << numEvents << "\n" 
       << "subrun = " << Subrun << "\n"
       << scientific << setprecision(5) << "\n"
       << "totalsigmaJet = " << totaljetXSec << "\n"
       << "sigmaJet = " << jetXSec << endl;

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
