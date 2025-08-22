#include "MARTINI.h"
#define ET_min 2.7
#define yMax 4.5
//#define ETgt_min 2.
//#define ETlt_min 2.
#define ETgt_min 100.
#define ETlt_min 25.
#define Rjet 0.4

int main(int argc, char* argv[]) 
{
  // create MARTINI object - it will automatically generate PYTHIA object and set its standard values
  // those values can be changed below before doing the explicit PYTHIA init.
  cout << "argc = " << argc << endl;
  stringstream file_n;
  file_n << argc;
  cout << "argc = " << file_n.str() << endl;

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
  int numJets = 1;               // Number of accepted jets required
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  double ETmax;                      // Maximum transverse energy of the fragmented nucleons
  double ETmax_opposite;             // Maximum transverse energy of the "opposite-side" nucleons
  int ETmax_i;                       // Number on the list of the highest energy hadron
  int ETmax_i_opposite;              // Number on the list of the highest energy "opposite-side" nucleon
  double E1, E2, ET_gt, ET_lt;       // E_> and E_<
  double phiLeadingNucleon;         // The azimuth of the highest energy nucleon
  double phiOppositeNucleon;        // The azimuth of the highest energy "opposite side" nucleon
  double Px1, Py1, Px2, Py2;        //The transverse momenta of the legs of the dijet
  double phi1, phi2;

  int nbins = 20;                  // Number of bins for dN/dA
  double dNdA_hist[nbins];          // The histogram for dN/dA
  double dNdphi_hist[nbins];        // The histogram for dN/dphi
  int fullEventSize;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
  }

  //Output the events:
  string e_name = "Events_"+file_n.str()+".dat";
  fstream foutEvents(e_name.c_str(),ios::out); 

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
      numEvents+=martini.fragmentation( plist );
      fullEvent = martini.pythia.event;
      
      //Now that the partons have undergone fragmentation, find the hadron of highest transverse energy:
      fullEventSize = fullEvent.size();
      ETmax = 0.; ETmax_i = 0;
      for(int ip=0; ip<fullEventSize; ip++){
	//Also in this loop, output the event in a text-based way:
	if(fullEvent[ip].isCharged() && fullEvent[ip].isHadron() ){
	  foutEvents << fullEvent[ip].m() << " " << fullEvent[ip].pT() << " "
		     << fullEvent[ip].phi() << " " << fullEvent[ip].y() << endl;
	}
	double yp = 0.5*log((fullEvent[ip].eCalc()+fullEvent[ip].pz())
			    /(fullEvent[ip].eCalc()-fullEvent[ip].pz()) );
	if( fullEvent[ip].mT() > ETmax
	    && fabs(yp) < yMax 
	    && fullEvent[ip].isHadron() 
	    && fullEvent[ip].isCharged()
	    && fullEvent[ip].isFinal() ){
	  ETmax = fullEvent[ip].mT();
	  ETmax_i = ip;
	}
      }
      foutEvents << "EndOfEvent" << endl;
      phiLeadingNucleon = fullEvent[ETmax_i].phi();
      
      //Here, we will look for the highest momentum hadron going in the "opposite transverse direction"
      //of the particle ETmax_i, and then determine E1 and E2 from cutting out cones around these two 
      //particles. This is all very much ad-hoc, but as long as this procedure works for a majority of 
      //the events, we are in business, statistically speaking:
      ETmax_opposite = 0.; ETmax_i_opposite = 0;
      for(int ip = 0; ip<fullEventSize; ip++){
	double yp = 0.5*log((fullEvent[ip].eCalc()+fullEvent[ip].pz())
			    /(fullEvent[ip].eCalc()-fullEvent[ip].pz()) );
	if(fullEvent[ip].mT() > ETmax_opposite && 
	   fabs(yp) < yMax 
	   && fullEvent[ip].isHadron()
	   && fullEvent[ip].isCharged() && fullEvent[ip].isFinal() ){
	  double dphi = max(fullEvent[ip].phi(), phiLeadingNucleon)
	    -min(fullEvent[ip].phi(), phiLeadingNucleon);
	  if(dphi > M_PI){ dphi = 2.*M_PI-dphi;}
	  if(dphi > 0.5*M_PI){
	      ETmax_opposite = fullEvent[ip].mT();
	      ETmax_i_opposite = ip;
	  }
	}
      }
      phiOppositeNucleon = fullEvent[ETmax_i_opposite].phi();

      //In a naive way, determine E_> and E_< :
      E1 = 0.; E2 = 0.;
      Px1 = 0.; Py1 = 0.; Px2 = 0.; Py2 = 0.;
      for(int ip = 0; ip<fullEventSize; ip++){
	double yp = 0.5*log((fullEvent[ip].eCalc()+fullEvent[ip].pz())
			    /(fullEvent[ip].eCalc()-fullEvent[ip].pz()) );
	if(fullEvent[ip].mT() > ET_min && 
	   fabs(yp) < yMax 
	   && fullEvent[ip].isHadron()
	   && fullEvent[ip].isCharged() && fullEvent[ip].isFinal() ){
	  double dphileading = max(fullEvent[ip].phi(), phiLeadingNucleon)
	    -min(fullEvent[ip].phi(), phiLeadingNucleon);
	  if(dphileading > M_PI){ dphileading = 2.*M_PI-dphileading;}

	  double dphiopposite = max(fullEvent[ip].phi(), phiOppositeNucleon)
	    -min(fullEvent[ip].phi(), phiOppositeNucleon);
	  if(dphiopposite > M_PI){ dphiopposite = 2.*M_PI-dphiopposite;}

	  if(dphileading < Rjet){ 
	    E1 += fullEvent[ip].mT();
	    Px1 += fullEvent[ip].pT()*cos(fullEvent[ip].phi());
	    Py1 += fullEvent[ip].pT()*sin(fullEvent[ip].phi());
	  }

	  if(dphiopposite < Rjet){ 
	    E2 += fullEvent[ip].mT();
	    Px2 += fullEvent[ip].pT()*cos(fullEvent[ip].phi());
	    Py2 += fullEvent[ip].pT()*sin(fullEvent[ip].phi());
	  }
	}
      }

      ET_gt = max(E1, E2); ET_lt = min(E1, E2);
      //If this energy triggers ATLAS, determine E_> and E_< :
      if(ET_gt > ETgt_min && ET_lt > ETlt_min){
	cout << ET_gt << " " << ET_lt << endl;
	
	//Put this event into the appropriate dN/dA bin:
	cout << " |E> - E<|/(E> + E<) = " << fabs(ET_gt-ET_lt)/(ET_gt+ET_lt) << endl;
	int proper_bin = (int)( ((double)nbins)*(fabs(ET_gt-ET_lt)/(ET_gt+ET_lt)) );
	cout << "proper_bin = " << proper_bin << endl;
	dNdA_hist[proper_bin] += 1./(double)(nbins*numJets);

	//And into the appropriate dN/dphi bin:
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
	
	double Deltaphi = max(phi1, phi2) - min(phi1, phi2);
	if(Deltaphi > M_PI){Deltaphi = 2.*M_PI-Deltaphi;}
	int phi_bin = (int)( ((double)nbins)*Deltaphi/M_PI);
	dNdphi_hist[phi_bin] += 1./(double)(nbins*numJets);

	//Finally, add to the counter:
	event_counter += 1;
      }
    }
  string A_name = "dNdA_"+file_n.str()+".dat";
  fstream foutdNdA(A_name.c_str(),ios::out); 
  string phi_name = "dNdphi_"+file_n.str()+".dat";
  fstream foutdNdphi(phi_name.c_str(), ios::out);
  for(int ib=0; ib<nbins; ib++)
    {
      foutdNdA << (double)ib/(double)nbins << " " << dNdA_hist[ib] << endl;
      foutdNdphi << M_PI*(double)ib/(double)nbins << " " << dNdphi_hist[ib] << endl;
    }
  foutdNdA.close();
  foutdNdphi.close();
  foutEvents.close();

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}

