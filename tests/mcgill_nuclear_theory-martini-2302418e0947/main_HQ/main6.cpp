#include "MARTINI.h"
#include <sstream>
#define ET_min 0.
#define etaMax 0.7
#define numJets 200000
#define nbins 500
#define ptmax 50.

using namespace std;

int main(int argc, char* argv[]) 
{
  MARTINI martini;

  martini.readFile("setup_main6.dat");
  martini.init(argc, argv);
  martini.settings.listAll();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  Event fullEvent;
  int numEvents;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  double charged_hist[nbins];          // The histogram for pi+

  //Make a tag for the output:
  string n_string = argv[1];
  n_string = "pp_higherpthat_"+n_string;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    charged_hist[ih] = 0.;
  }

  //Output the events:
  string e_name = "RAAEvents_"+n_string+".dat";
  fstream foutEvents( e_name.c_str() ,ios::out); 

  while(event_counter<numJets)      // The main loop
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
	//Finally, bin the charged pions:
	int eventSize = fullEvent.size();
	for(int ip=0; ip<eventSize; ip++){
	  if(fullEvent[ip].mT() > ET_min
	     && fabs(fullEvent[ip].eta()) < etaMax
	     && fullEvent[ip].isFinal() 
	     && fullEvent[ip].isCharged()){
	    int pbin = (int)((double)nbins*fullEvent[ip].pT()/ptmax);
	    if(pbin >= 0 && pbin < nbins){
	      charged_hist[pbin] += (double)nbins/(double)numJets;
	    }
	  }
	  if(fullEvent[ip].isFinal()){
	    //foutEvents << fullEvent[ip].m() << " " << fullEvent[ip].pT() << " " << fullEvent[ip].phi() << " " << fullEvent[ip].eta() << endl;
	  }
	}
	//foutEvents << "EndOfEvent" << endl;
      }
    }
  
  string charged_name = "charged_"+n_string+".dat";
  fstream foutcharged( charged_name.c_str() ,ios::out);   
  
  for(int ib=0; ib<nbins; ib++){
    foutcharged << ptmax*(double)ib/(double)nbins << " " << charged_hist[ib] << endl;
  }
  
  foutEvents.close();
  foutcharged.close();
  
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
