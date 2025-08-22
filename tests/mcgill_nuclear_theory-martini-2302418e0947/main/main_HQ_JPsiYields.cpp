#include "MARTINI.h"
#include <sstream>
#define numJets 5000
#define nybins 100
#define DNDY_MAX 5.
#define nptbins 50

//#define Y_MAX 0.35
#define R_MAX 15.
#define DP_MAX 15.
#define Y_MAX 0.35
#define PT_MAX 5.

using namespace std;

int main(int argc, char* argv[]) 
{
  ofstream yieldsOut;
  yieldsOut.open("N_QbarQ_evolutionOn_nuclearEffectsOn_1S1Mass_notMovedBeforeTau0.dat");
  double dN4pir2dr[nptbins], dN4pip2dp[nptbins];
  for(int ipt=0; ipt<nptbins; ipt++){
    dN4pir2dr[ipt] = dN4pip2dp[ipt] = 0.;
  }
  for(int ib=1; ib<13; ib++){
    cout << "ib = " << ib << endl;
    string setup_tag = "setup_files/setup_main_HQ_JPsiYields_b";
    stringstream bstring;
    bstring << ib;
    string setup_name = setup_tag + bstring.str() + ".dat";
    
    //Counting the number of bound states:
    double NJPsi_diagonal, NExcitedJPsi_diagonal, NUpsilon_diagonal, NExcitedUpsilon_diagonal, NBc_diagonal, NExcitedBc_diagonal;
    double NJPsi_recombined, NExcitedJPsi_recombined, NUpsilon_recombined, NExcitedUpsilon_recombined, NBc_recombined, NExcitedBc_recombined;
    NJPsi_diagonal = NExcitedJPsi_diagonal = 0.;
    NUpsilon_diagonal = NExcitedUpsilon_diagonal = 0.;
    NBc_diagonal = NExcitedBc_diagonal = 0.;
    NJPsi_recombined = NExcitedJPsi_recombined = 0.;
    NUpsilon_recombined = NExcitedUpsilon_recombined = 0.;
    NBc_recombined = NExcitedBc_recombined = 0.;
    
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
	//cout << "event_counter = " << event_counter << endl;
	
	plist[0]->clear();          // clear the parton list
	plist[1]->clear();          // clear the parton list
	
	martini.generateEvent_HQ(plist[0]); 
	//cout << "OK after generateEvent_HQ..." << endl;
	int eventSize = plist[0]->size();
	//cout << "eventSize = " << eventSize << endl;
	
	int counter = 0;
	if (martini.returnEvolution() == 1)
	  {
	    for(int i=0; i<mt; i++)
	      {
		counter = martini.evolve_HQ(plist, counter, i);
		counter+=1;
	      }
	  }
	//cout << "OK after evolution..." << endl;
	eventSize = plist[0]->size();
        //cout << "eventSize = " << eventSize << endl ;
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
	  //cout << "Particle " << ip << "'s tFinal = " << plist[0]->at(ip).tFinal() << endl;
	  //cout << "Delta p = " << dp << endl;
	  //cout << "r = " << r << endl;
	  //int rbin = (int) (nptbins*r/RMAX);
	  //int dpbin = (int) (nptbins*dp/DPMAX);
	  //if(rbin < nptbins){
	  //dN4pir2dr[rbin] += ((double)nptbins/RMAX)/(4.*M_PI*RMAX*RMAX*((double)rbin+0.5)*((double)rbin+0.5));
	  //}
	  //if(dpbin < nptbins){
	  //dN4pip2dp[dpbin] += ((double)nptbins/DPMAX)/(4.*M_PI*DPMAX*DPMAX*((double)dpbin+0.5)*((double)dpbin+0.5));
	  //}
	}
	martini.fragmentation_HQ(plist);
	//cout << "OK after fragmentation_HQ..." << endl;
	
	eventSize = martini.pythia.event.size();
	//cout << "eventSize = " << eventSize << endl;
	
	for(int ie=0; ie<eventSize; ie++){
	  Vec4 itsP = martini.pythia.event[ie].p();
	  double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	  double itsM = itsP.mCalc();
	  //double itsM = martini.pythia.event[ie].mass();
	  //cout << "itsM: " << itsM << endl;
	  double pz = itsP.pz();
	  double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );
	  int status = martini.pythia.event[ie].status();
	  
	  if(fabs(y)<Y_MAX && pt<PT_MAX){
	    if(martini.pythia.event[ie].id()==443){
	      if(status == 81) NJPsi_diagonal += 1.;
	      else NJPsi_recombined += 1.;
	    }
	    if(martini.pythia.event[ie].id()==443000){
	      if(status == 81) NExcitedJPsi_diagonal += 1.;
	      else NExcitedJPsi_recombined += 1.;
	    }
	    if(martini.pythia.event[ie].id()==553){
	      if(status == 81) NUpsilon_diagonal += 1.;
	      else NUpsilon_recombined += 1.;
	    }
	    if(martini.pythia.event[ie].id()==553000){
	      if(status == 81) NExcitedUpsilon_diagonal += 1.;
	      NExcitedUpsilon_recombined += 1.;
	    }
	    if(martini.pythia.event[ie].id()==543){
	      if(status == 81) NBc_diagonal += 1.;
	      else NBc_recombined += 1.;
	    }
	    if(martini.pythia.event[ie].id()==543000){
	      if(status == 81) NExcitedBc_diagonal += 1.;
	      else NExcitedBc_recombined += 1.;
	    }
	  }
	}
	//cout << "OK after binning..." << endl;
      }
    
    //Output the properly normalized numbers of quarkonia:
    NJPsi_diagonal /= (double)numJets;
    NExcitedJPsi_diagonal /= (double)numJets;
    NUpsilon_diagonal /= (double)numJets;
    NExcitedUpsilon_diagonal /= (double)numJets;
    NBc_diagonal /= (double)numJets;
    NExcitedBc_diagonal /= (double)numJets;
    NJPsi_recombined /= (double)numJets;
    NExcitedJPsi_recombined /= (double)numJets;
    NUpsilon_recombined /= (double)numJets;
    NExcitedUpsilon_recombined /= (double)numJets;
    NBc_recombined /= (double)numJets;
    NExcitedBc_recombined /= (double)numJets;
    
    yieldsOut << (double)ib << " " << NJPsi_diagonal << " " << NJPsi_recombined << " " << NExcitedJPsi_diagonal << " " << NExcitedJPsi_recombined << " " 
	      << NUpsilon_diagonal << " " << NUpsilon_recombined  << " " << NExcitedUpsilon_diagonal << " " << NExcitedUpsilon_recombined  << " " 
	      << NBc_diagonal << " " << NBc_recombined  << " " << NExcitedBc_diagonal << " " << NExcitedBc_recombined << endl;

    delete plist[0];
    delete plist[1];
    delete plist;  
    delete plistInitial;
  }
}
