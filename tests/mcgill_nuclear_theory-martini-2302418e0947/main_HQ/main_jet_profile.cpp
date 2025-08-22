#include "MARTINI.h"
#include <sstream>
#define ET_min 0.
#define yMax 0.3
#define nptbins 30
#define PT_MAX 15.
#define numJets 1000000

using namespace std;

int main(int argc, char* argv[]) 
{
  MARTINI martini;
  Event fullEvent;

  martini.readFile("setup_main_jet_profile.dat");
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
  n_string = "pp_276TeV_ptmin0_factMultFac10_renormMultFac10_"+n_string;
  ////Output the events:
  //string e_name = dir_string+"/improvedphotons_"+n_string+".dat";
  //string e_name_f = dir_string+"/improvedphotons_fragmented_"+n_string+".dat";
  //string e_name = "partons_"+n_string+".dat";
  //fstream foutEvents( e_name.c_str() ,ios::out); 
  //string f_name = "neutralpions_"+n_string+".dat";
  //fstream foutFEvents( f_name.c_str() ,ios::out); 

  //Start binning at runtime now:
  double dN2piptdptdy_photons[nptbins];
  for(int ipt=0; ipt<nptbins; ipt++){
    dN2piptdptdy_photons[ipt] = 0.;
  }

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      int totaljets = 0;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      totaljets = martini.generateEvent(plist[0]); 
	  // version that samples number of collisions with Glauber model
      
      // Pre-evolution jet profile:
      int eventSize = plist[0]->size();
      for(int ip=0; ip<eventSize; ip++){
	itsP = plist[0]->at(ip).p();
	double E = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py() +itsP.pz()*itsP.pz() );
	double y = 0.5*log(( E+itsP.pz() )/( E-itsP.pz() ) );
	if(fabs(y)<yMax && plist[0]->at(ip).id()==22){
	  //foutEvents << plist[0]->at(ip).id() << " " << plist[0]->at(ip).status() << " " << itsP.px() << " " << itsP.py() << " " << itsP.pz() << endl;
	  int ptbin = (int)(nptbins*sqrt(itsP.px()*itsP.px()+itsP.py()*itsP.py() )/PT_MAX);
	  if(ptbin < nptbins){
	    dN2piptdptdy_photons[ptbin] += (nptbins/PT_MAX)*(0.5/yMax)/(numJets*2.*M_PI*PT_MAX*((double)ptbin+0.5)/((double)nptbins));
	  }
	}
      }
      //foutEvents << "EndOfEvent" << endl;
    
      //int counter = 0;
      //if (martini.returnEvolution() == 1)         // evolve in medium if settings want that
      //{
      //  for(int i=0; i<mt; i++) // loop over all time steps 
      //    {
      //      counter = martini.evolve(plist, counter, i);
      //      counter+=1;
      //    }
      //}

      ////Now fragment:
      //fullEvent.clear();
      //martini.pythia.event.clear();
      //for(int ijets=0; ijets<totaljets; ijets++){
      //int frag_success = martini.fragmentation( plist , ijets);
      //if(frag_success){
      //  fullEvent = martini.pythia.event;
      //  //Finally, record the neutral pions:
	//  eventSize = fullEvent.size();
      //for(int ip=0; ip<eventSize; ip++){
      //    if(fullEvent[ip].id() == 111 ){
      //      //foutFEvents << fullEvent[ip].status() << " " << fullEvent[ip].px() << " " << fullEvent[ip].py() << " " 
	//      //	  << fullEvent[ip].pz() << endl;
      //  }
      //  }
      //}
      //}
      //foutFEvents << "EndOfEvent" << endl;
    }
  string photon_name = "dN2piptdptdy_photons_"+n_string+".dat";
  fstream foutPhotons( photon_name.c_str() ,ios::out); 
  for(int ipt=0; ipt<nptbins; ipt++){
    foutPhotons << ipt*PT_MAX/nptbins << " " << dN2piptdptdy_photons[ipt] << endl;
  }
  
  
  //foutEvents.close();
  //foutFEvents.close();
  
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
