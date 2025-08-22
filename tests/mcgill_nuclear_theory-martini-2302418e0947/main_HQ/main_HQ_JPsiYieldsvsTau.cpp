#include "MARTINI.h"
#include <sstream>
#define numJets 2000
#define nybins 100
#define DNDY_MAX 5.
#define nptbins 50
#define Y_MAX 0.35
#define PT_MAX 5.

using namespace std;

int main(int argc, char* argv[]) 
{
  MARTINI martini;
  Event fullEvent;
  
  martini.readFile("setup_files/setup_main_HQ_JPsiYieldsvsTau.dat");
  martini.init(argc, argv);
  martini.settings.listAll();
  
  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;
  
  int numEvents;
  int event_counter = 0;            // Counter of events
  
  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  
  double NJPsi_diagonal[mt], NJPsi_recombinant[mt];
  for(int it=0; it<mt; it++){
    NJPsi_diagonal[it] = NJPsi_recombinant[it] = 0.;
  }
  
  while(event_counter<numJets)      // The main loop
    {
      event_counter++;
      cout << "event_counter = " << event_counter << endl;
      
      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list
      
      martini.generateEventHeavyQuarks(plist[0]); 
      cout << "OK after generateEvent_HQ..." << endl;
      int eventSize = plist[0]->size();
      cout << "eventSize = " << eventSize << endl;
      
      int counter = 0;
      if (martini.returnEvolution() == 1)
	{
	  for(int i=0; i<mt; i++)
	    {
	      counter = martini.evolveHeavyQuarks(plist, counter, i);
	      martini.hadronizeHeavyQuarks(plist);
	      eventSize = martini.pythia.event.size();
	      for(int ie=0; ie<eventSize; ie++){
		Vec4 itsP = martini.pythia.event[ie].p();
		double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
		double itsM = itsP.mCalc();
		double pz = itsP.pz();
		double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );
		int status = martini.pythia.event[ie].status();
		
		if(fabs(y)<Y_MAX && pt<PT_MAX){
		  if(martini.pythia.event[ie].id()==443){
		    if(status == 81) NJPsi_diagonal[i] += 1.;
		    else NJPsi_recombinant[i] += 1.;
		  }
		  //if(martini.pythia.event[ie].id()==443000){
		  //if(status == 81) NExcitedJPsi_diagonal += 1.;
		  //else NExcitedJPsi_recombined += 1.;
		  //}
		  //if(martini.pythia.event[ie].id()==553){
		  //if(status == 81) NUpsilon_diagonal += 1.;
		  //else NUpsilon_recombined += 1.;
		  //}
		  //if(martini.pythia.event[ie].id()==553000){
		  //if(status == 81) NExcitedUpsilon_diagonal += 1.;
		  //NExcitedUpsilon_recombined += 1.;
		  //}
		  //if(martini.pythia.event[ie].id()==543){
		  //if(status == 81) NBc_diagonal += 1.;
		  //else NBc_recombined += 1.;
		  //}
		  //if(martini.pythia.event[ie].id()==543000){
		  //if(status == 81) NExcitedBc_diagonal += 1.;
		  //else NExcitedBc_recombined += 1.;
		  //}
		}
	      }
	      counter+=1;
	    }
	}
    }
  
  //Output the properly normalized numbers of quarkonia:
  for(int it=0; it<mt; it++){
    NJPsi_diagonal[it] /= (double)numJets;
    NJPsi_recombinant[it] /= (double)numJets;
    cout << "At t = " << it*dtfm << " fm/c, NJPsi_diagonal = " << NJPsi_diagonal[it] << " and NJPsi_recombinant = " << NJPsi_recombinant[it] << endl;
  }
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
