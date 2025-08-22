#include "MARTINI.h"
#include <sstream>
#define ET_min 0.7
#define etaMax 4.5
#define ETgt_min 100.
#define etaLeadingJetMax 2.8
#define ETlt_min 25.
#define Rjet 0.4
#define numJets 800
#define nbins 20

using namespace std;

//We are now using the anti-kt algorithm that most, if not all, of the 
//LHC collaborations are using:
int main(int argc, char* argv[]) 
{
  MARTINI martini;

  martini.readFile("setup_dijets_pp_celljet.dat");
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

  //Energies of the jets:
  double ET_gt, ET_lt;
  int i_greatest_ET, i_greatest_opposite_ET;
  double phiLeadingJet;
  double etaLeadingJet;

  double dNdA_hist[nbins];          // The histogram for dN/dA
  double dNdphi_hist[nbins];        // The histogram for dN/dphi
  int fullEventSize;

  //Make a tag for the output:
  string n_string = argv[1];
  n_string = "pp_partons_celljet_Qin_pthatmin25_"+n_string;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
  }

  //Output the events:
  string e_name = "Events_"+n_string+".dat";
  fstream foutEvents( e_name.c_str() ,ios::out); 

  CellJet celljet = CellJet(etaMax, 90, 63, 2);
  //CellJet celljet = CellJet(etaMax);
  //CellJet celljet = CellJet(etaMax, 900, 630, 2);

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

      cout << "Size of plist[0] = " << plist[0]->size() << endl;
      celljet.analyze(martini.pythia.event, 25., 0.4, 1.5, cout);
      
      //Next, find the jet with highest transverse energy:
      ET_gt = 0.; i_greatest_ET = 0;
      int jetsSize = celljet.size();
      cout << "jetsSize = " << jetsSize << endl;
      for(int ij=0; ij<jetsSize; ij++){
	if(celljet.eT(ij) > ET_gt){
	  ET_gt = celljet.eT(ij);
	  i_greatest_ET = ij;
	  phiLeadingJet = celljet.phiWeighted(ij);
	  etaLeadingJet = celljet.etaWeighted(ij);
	}
      }

      //Next, find the dijet of greatest energy on the opposite azimuth of the leading
      //jet:
      ET_lt = 0.; i_greatest_opposite_ET = 0;
      for(int ij=0; ij<jetsSize; ij++){
	double delta_phi = max(phiLeadingJet, celljet.phiWeighted(ij))
	  -min(phiLeadingJet, celljet.phiWeighted(ij) );
	if(delta_phi > M_PI){delta_phi = 2.*M_PI-delta_phi;}
	if(celljet.eT(ij) > ET_lt
	   && delta_phi > 0.5*M_PI){
	  ET_lt = celljet.eT(ij);
	  i_greatest_opposite_ET = ij;
	}
      }
      
      //cout << "ET_gt = " << ET_gt << ", ET_lt = " << ET_lt << endl;     

      //Now that the highest energy jet and the highest energy opposing jet has been determined, 
      //determine if they make it into the ATLAS analysis and bin them accordingly:

      if(ET_gt > ETgt_min && ET_lt > ETlt_min && fabs(etaLeadingJet) < etaLeadingJetMax){
	
	//Put this event into the appropriate dN/dA bin:
	cout << " |E> - E<|/(E> + E<) = " << fabs(ET_gt-ET_lt)/(ET_gt+ET_lt) << endl;
	int proper_bin = (int)( ((double)nbins)*(fabs(ET_gt-ET_lt)/(ET_gt+ET_lt)) );
	//cout << "proper_bin = " << proper_bin << endl;
	if(proper_bin>=0 && proper_bin<nbins){
	  dNdA_hist[proper_bin] += (double)nbins/(double)numJets;}
	
	double Deltaphi = max(phiLeadingJet, celljet.phiWeighted(i_greatest_opposite_ET) ) 
	  - min(phiLeadingJet, celljet.phiWeighted(i_greatest_opposite_ET) );
	if(Deltaphi > M_PI){Deltaphi = 2.*M_PI-Deltaphi;}
	int phi_bin = (int)( ((double)nbins)*Deltaphi/M_PI);
	if(phi_bin>=0 && phi_bin<nbins){
	  dNdphi_hist[phi_bin] += (double)nbins/(double)numJets;}
	
	//Finally, add to the counter:
	event_counter += 1;
      }
    }
  
  string A_name = "dNdA_"+n_string+".dat";
  fstream foutdNdA(A_name.c_str(),ios::out); 
  string phi_name = "dNdphi_"+n_string+".dat";
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
