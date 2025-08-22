#include "MARTINI.h"
#include <sstream>
#define numJets 1000
#define nybins 100
#define DNDY_MAX 5.
#define nptbins 50

#define R_MAX 15.
#define DP_MAX 15.
//#define Y_MAX 0.35
#define Y_MAX 2.
//#define Y_MAX 5.
#define PT_MAX 5.

using namespace std;

int main(int argc, char* argv[]) 
{
  double yields[10][12];
  for(int ib=0; ib<10; ib++){
    for(int iy=0; iy<12; iy++){
      yields[ib][iy] = 0.;
    }
  }


  for(int ib=2; ib<11; ib+=2){
    cout << "ib = " << ib << endl;
    string setup_tag = "setup_files/setup_main_HQ_JPsiYields_LHC_b";
    stringstream bstring;
    bstring << ib;
    string setup_name = setup_tag + bstring.str() + ".dat";
    
    MARTINI martini;
    Event fullEvent;
    
    martini.readFile(setup_name.c_str() );
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
    
    while(event_counter<numJets)      // The main loop
      {
	event_counter++;
	cout << "event_counter = " << event_counter << endl;
	
	plist[0]->clear();          // clear the parton list
	plist[1]->clear();          // clear the parton list
	
	martini.generateEventHeavyQuarks(plist[0]); 
	cout << "OK after generateEvent_HQ..." << endl;
	int eventSize = plist[0]->size();
	//cout << "eventSize = " << eventSize << endl;
	
	int counter = 0;
	//if (martini.returnEvolution() == 1)
	//{
	    for(int i=0; i<mt; i++)
	      {
		//counter = martini.evolveHeavyQuarks(plist, counter, i);
		counter = martini.evolveAndHadronizeHeavyQuarks(plist, counter, i);
		counter+=1;
	      }
	    //}
	    //cout << "OK after evolution and hadronization..." << endl;
	eventSize = plist[0]->size();
        cout << "eventSize = " << eventSize << endl ;
	for(int ip=0; ip<eventSize; ip++){
	  double dx = plist[0]->at(ip).x()-plist[0]->at(plist[0]->at(ip).antiI() ).x();
	  double dy = plist[0]->at(ip).y()-plist[0]->at(plist[0]->at(ip).antiI() ).y();
	  double dz = plist[0]->at(ip).z()-plist[0]->at(plist[0]->at(ip).antiI() ).z();
	  Vec4 P = plist[0]->at(ip).p();
	  Vec4 antiP = plist[0]->at(plist[0]->at(ip).antiI() ).p();
	  double r = sqrt(dx*dx+dy*dy+dz*dz );
	  double dp = sqrt( (P.px()-antiP.px() )*(P.px()-antiP.px() ) 
			   +(P.py()-antiP.py() )*(P.py()-antiP.py() ) 
			   +(P.pz()-antiP.pz() )*(P.pz()-antiP.pz() ) );
	}

	if (martini.returnEvolution() == 0)                                                                                                                                  
	  {
	    martini.hadronizeHeavyQuarks(plist);
	  }
	
	eventSize = martini.pythia.event.size();
	//cout << "eventSize = " << eventSize << endl;
	
	for(int ie=0; ie<eventSize; ie++){
	  Vec4 itsP = martini.pythia.event[ie].p();
	  double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	  double itsM = itsP.mCalc();
	  double pz = itsP.pz();
	  //double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );
	  //Pseudorapidity is more appropriate here:
	  double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	  int status = martini.pythia.event[ie].status();
	  
	  if(fabs(y)<Y_MAX && pt<PT_MAX){
	    if(martini.pythia.event[ie].id()==443){
	      if(status == 81) yields[ib-1][0] += 1.;
	      else yields[ib-1][1] += 1.;
	    }
	    if(martini.pythia.event[ie].id()==443000){
	      if(status == 81) yields[ib-1][2] += 1.;
	      else yields[ib-1][3] += 1.;
	    }
	    if(martini.pythia.event[ie].id()==553){
	      if(status == 81) yields[ib-1][4] += 1.;
	      else yields[ib-1][5] += 1.;
	    }
	    if(martini.pythia.event[ie].id()==553000){
	      if(status == 81) yields[ib-1][6] += 1.;
	      else yields[ib-1][7] += 1.;
	    }
	    if(martini.pythia.event[ie].id()==543){
	      if(status == 81) yields[ib-1][8] += 1.;
	      else yields[ib-1][9] += 1.;
	    }
	    if(martini.pythia.event[ie].id()==543000){
	      if(status == 81) yields[ib-1][10] += 1.;
	      else yields[ib-1][11] += 1.;
	    }
	  }
	}
      }
    
    //Output the properly normalized numbers of quarkonia:
    for(int iy=0;iy<12; iy++){
      yields[ib-1][iy] /= (double)numJets;
    }
    
    delete plist[0];
    delete plist[1];
    delete plist;  
    delete plistInitial;
  }
  stringstream s;
  s << argv[1];

  string directoryName = "HQTestingAndResults/";

  //string yieldsName = directoryName+"N_QbarQ_RHIC_evolutionOff_nuclearEffectsOff_moveBeforeTau02_CMFrame_dt00045_"+s.str()+".dat";
  string yieldsName = directoryName+"N_QbarQ_CMS_rescaledpp_"+s.str()+".dat";
  ofstream yieldsOut;
  yieldsOut.open(yieldsName.c_str() );

  for(int ib=2; ib<11; ib+=2){
    yieldsOut << (double)ib;
    for(int iy=0; iy<12; iy++){
      yieldsOut << " " << yields[ib-1][iy];
    }
    yieldsOut << endl;
  }

}
