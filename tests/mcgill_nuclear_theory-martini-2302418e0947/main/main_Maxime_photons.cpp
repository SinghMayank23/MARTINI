#include "MARTINI.h"
#include <sstream>
#define ET_min 0.
#define etaMax 10.
#define numJets 1
#define nbins 50
#define ptmax 10.

using namespace std;

int main(int argc, char* argv[]) 
{
  MARTINI martini;
  Event fullEvent;

  martini.readFile("setup_main_Maxime_photons.dat");
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

  Vec4 itsP;                          // Its P

  //Make a tag for the output:
  string n_string = argv[1];
  string dir_string = argv[3];
  n_string = "test_ptmin1_"+n_string;
  ////Output the events:
  //string e_name = dir_string+"/improvedphotons_"+n_string+".dat";
  //string e_name_f = dir_string+"/improvedphotons_fragmented_"+n_string+".dat";
  string e_name = "improvedphotons_"+n_string+".dat";
  string e_name_f = "improvedphotons_fragmented_"+n_string+".dat";
  fstream foutEvents( e_name.c_str() ,ios::out); 
  fstream foutFEvents( e_name_f.c_str() ,ios::out); 

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      int totaljets = 0;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      totaljets = martini.generateEvent(plist[0]); 
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

      // Fragmentation is not needed or wanted here, operate directly on plist:
      int eventSize = plist[0]->size();
      for(int ip=0; ip<eventSize; ip++){
	itsP = plist[0]->at(ip).p();
	if(plist[0]->at(ip).id() == 22){
	  foutEvents << plist[0]->at(ip).status() << " " << plist[0]->at(ip).source() << " " << itsP.px() << " " << itsP.py() << " " << itsP.pz() << endl;
	}
      }
      foutEvents << "EndOfEvent" << endl;
    
      //Now fragment:
      fullEvent.clear();
      martini.pythia.event.clear();
      for(int ijets=0; ijets<totaljets; ijets++){
	int frag_success = martini.fragmentation( plist , ijets);
	if(frag_success){
	  fullEvent = martini.pythia.event;
	  //Finally, bin the charged pions:
	  eventSize = fullEvent.size();
	  for(int ip=0; ip<eventSize; ip++){
	    if(fullEvent[ip].id() == 22
	       && (fullEvent[ip].status() < 90 || fullEvent[ip].status() > 99) ){ //Exclude decay photons
	      //&& fullEvent[ip].isFinal() ){
	      foutFEvents << fullEvent[ip].status() << " " << fullEvent[ip].px() << " " << fullEvent[ip].py() << " " 
			  << fullEvent[ip].pz() << endl;
	    }
	  }
	}
      }
      foutFEvents << "EndOfEvent" << endl;
    }
  
  
  foutEvents.close();
  foutFEvents.close();
  
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
