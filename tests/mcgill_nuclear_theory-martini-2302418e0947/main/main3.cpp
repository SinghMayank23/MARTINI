#include "MARTINI.h"
#include <sstream>
#define ET_min 2.7
//#define ET_min 2.
#define yMax 4.5
#define ETgt_min 100.
#define ETlt_min 25.
#define Rjet 0.4

// This main program runs MARTINI events and bins them appropiately for comparisons with the analysis
// of ATLAS. Instead of working with "jet cones", we simply split the azimuthal range into two halves, 
// where one half is centered around the leading hadron's azimuth. -CFY 12/13/2010

int main(int argc, char* argv[]) 
{
  // Create MARTINI object - it will automatically generate PYTHIA object and set its standard values
  // those values can be changed below before doing the explicit PYTHIA initialization:
  MARTINI martini;

  martini.readFile("setup.dat");
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
  int numJets = 500;               // Number of accepted jets required
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  double ETmax;                     // Maximum transverse energy of the fragmented nucleons
  int ETmax_i;                      // Number on the list of highest energy hadron
  double phiTriggerNucleon;         // The azimuth of the highest energy nucleon

  double ETmax_opp;                     // Maximum transverse energy of the "opposite-side" nucleon
  int ETmax_i_opp;                      // Number on the list of opposing nucleon
  double phiOppositeNucleon;         // The azimuth of the opposing nucleon

  double Px1, Py1, Px2, Py2;        //Used to determine \Delta \phi for the dijets
  double E1, E2, ET_gt, ET_lt;      // E_> and E_<
  double e1_calc, e2_calc, pz1, pz2; //For calculating the jet rapidities
  double phi1, phi2, delta_phi;

  int nbins = 30;                   // Number of bins
  double dNdA_hist[nbins];          // The histogram for dN/dA
  double dNdphi_hist[nbins];        // The histogram for dN/d\phi
  int fullEventSize;

  stringstream file_n;
  file_n << argc;

  string e_name = "Event_"+file_n.str()+".dat";
  fstream foutEvent(e_name.c_str(),ios::out); //For outputting the entire event

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
  }

  while(event_counter<numJets)     // The main loop
    {

      plist[0]->clear();           // clear the parton list
      plist[1]->clear();           // clear the parton list

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
      numEvents+=martini.fragmentation( plist );
      fullEvent = martini.pythia.event;

      //Now that the partons have undergone fragmentation, find the hadron of highest transverse energy:
      fullEventSize = fullEvent.size();
      ETmax = 0.; ETmax_i = 0;
      for(int ip=0; ip<fullEventSize; ip++){
	double yp = 0.5*log((fullEvent[ip].eCalc()+fullEvent[ip].pz())
			    /(fullEvent[ip].eCalc()-fullEvent[ip].pz()) );
	if( fullEvent[ip].mT() > ETmax
	    && fabs(yp) < yMax 
	    && fullEvent[ip].isHadron() && fullEvent[ip].isCharged()
	    && fullEvent[ip].isFinal() ){
	  ETmax = fullEvent[ip].mT();
	  ETmax_i = ip;
	}
      }
      phiTriggerNucleon = fullEvent[ETmax_i].phi();
      
      //In a naive way, determine E_> and E_< :
      Px1 = 0.; Py1 = 0.; Px2 = 0.; Py2 = 0.;
      E1 = 0.; E2 = 0.;
      for(int ip = 0; ip<fullEventSize; ip++){

	double yp = 0.5*log((fullEvent[ip].eCalc()+fullEvent[ip].pz())
			    /(fullEvent[ip].eCalc()-fullEvent[ip].pz()) );

	if(fullEvent[ip].mT() > ET_min && 
	   fabs(yp) < yMax &&
	   fullEvent[ip].isHadron()
	   && fullEvent[ip].isCharged() && fullEvent[ip].isFinal() ){

	  double dphi = max(fullEvent[ip].phi(), phiTriggerNucleon)
	    -min(fullEvent[ip].phi(), phiTriggerNucleon);
	  if(dphi > M_PI){ dphi = 2.*M_PI-dphi;}

	  if(dphi < 0.5*M_PI){
	    E1  += fullEvent[ip].mT();
	    Px1 += fullEvent[ip].pT()*cos(fullEvent[ip].phi());
	    Py1 += fullEvent[ip].pT()*sin(fullEvent[ip].phi());
	  }
	  else{
	    E2 += fullEvent[ip].mT();
	    Px2 += fullEvent[ip].pT()*cos(fullEvent[ip].phi());
	    Py2 += fullEvent[ip].pT()*sin(fullEvent[ip].phi());
	  }
	}
      }
      ET_gt = max(E1, E2); ET_lt = min(E1, E2);
      //If this energy triggers ATLAS, determine E_> and E_< :
      if(ET_gt > ETgt_min && ET_lt > ETlt_min){

	//Output all final charged particles:
	for(int ip = 0; ip<fullEventSize; ip++){
	  if(fullEvent[ip].isCharged() && fullEvent[ip].isFinal()){
	    foutEvent << fullEvent[ip].m() << " " << fullEvent[ip].pT() << " " << fullEvent[ip].y() << endl;
	  }
	}
	foutEvent << "EndofEvent" << endl;

	cout << ET_gt << " " << ET_lt << endl;

	if(Px1 > 0.){
	  phi1 = atan(Py1/Px1);
	}
	else{
	  phi1 = M_PI+atan(Py1/Px1);
	}
	if(phi1<0.){phi1 += 2.*M_PI;}
	if(Px2 > 0.){
	  phi2 = atan(Py2/Px2);
	}
	else{
	  phi2 = M_PI+atan(Py2/Px2);
	}
	if(phi2<0.){phi2 += 2.*M_PI;}

	//Put this event into the appropriate dN/dA bin:
	cout << " |E> - E<|/(E> + E<) = " << fabs(ET_gt-ET_lt)/(ET_gt+ET_lt) << endl;
	int A_bin = (int)( ((double)nbins)*(fabs(ET_gt-ET_lt)/(ET_gt+ET_lt)) );
	cout << "A_bin = " << A_bin << endl;
	dNdA_hist[A_bin] += (double)nbins/(double)numJets;

	//Next, put this event into the appropriate dN/d\phi bin:
	delta_phi = max(phi1, phi2)-min(phi1, phi2);
	if(delta_phi > M_PI){delta_phi = 2.*M_PI-delta_phi;}
	int phi_bin = (int)( (double)nbins*delta_phi/M_PI );
	dNdphi_hist[phi_bin] += (double)nbins/(double)numJets;
	
	//Finally, add to the counter:
	event_counter += 1;
      }
      //else{ cout << "Event did not trigger ATLAS" << endl;}
    }

  //Finally, output the histograms. We really need to think of more information to output here and how:
  string A_name = "dNdA_"+file_n.str()+".dat";
  fstream foutdNdA(A_name.c_str(),ios::out); 
  for(int ib=0; ib<nbins; ib++)
    {
      foutdNdA << (double)ib/(double)nbins << " " << dNdA_hist[ib] << endl;
    }
  foutdNdA.close();

  string phi_name = "dNdphi_"+file_n.str()+".dat";  
  fstream foutdNdphi(phi_name.c_str(),ios::out); 
  for(int ib=0; ib<nbins; ib++)
    {
      foutdNdphi << M_PI*(double)ib/(double)nbins << " " << dNdphi_hist[ib] << endl;
    }
  foutdNdA.close();
  
  foutEvent.close();
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}

