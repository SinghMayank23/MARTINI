#include "MARTINI.h"
#include "fastjet/ClusterSequence.hh"
#include <sstream>
#include <math.h>
#define ET_min 0.1
#define etaMax 2.
#define ETgt_min 120.
#define etaLeadingJetMax 2.
#define ETlt_min 50.
#define Rjet 0.5
#define B_IN_TESLA 3.8
#define RD_IN_METERS 1.
#define numJets 5
#define nbins 20

using namespace fastjet;
using namespace std;

//Rewritten for the analysis at the CMS:

//We are now using the anti-kt algorithm that most, if not all, of the 
//LHC collaborations are using:
int main(int argc, char* argv[]) 
{
  // create MARTINI object - it will automatically generate PYTHIA object and set its standard values
  // those values can be changed below before doing the explicit PYTHIA init.
  MARTINI martini;

  martini.readFile("setup_main_CMS_dNdA.dat");
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
  n_string = "PbPb_CMS_alphas27"+n_string;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
  }

  //Output the events:
  string e_name = "Events_"+n_string+".dat";
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
      numEvents+=martini.fragmentation( plist );
      fullEvent = martini.pythia.event;
      
      //Now that the partons have undergone fragmentation, place all final, charged hadrons detected by 
      //ATLAS into a PseudoJet vector:
      fullEventSize = fullEvent.size();

      vector<PseudoJet> particles;
      for(int ip=0; ip<fullEventSize; ip++){
	//Charged particles have a minimum p_t, neutral particles do not:
	if(fabs(fullEvent[ip].eta()) < etaMax 
	   //&& ( (fullEvent[ip].isCharged() && fullEvent[ip].mT() > ET_min) 
	   //	|| !(fullEvent[ip].isCharged()) )
	   && fullEvent[ip].mT() > ET_min
	   && fullEvent[ip].isFinal() 
	   && fullEvent[ip].pT() > 0.149896*fabs(fullEvent[ip].charge())*B_IN_TESLA*RD_IN_METERS ){
	  double Ep = fullEvent[ip].eCalc();
	  double azimuth = fullEvent[ip].phi()+asin( 0.149896*fullEvent[ip].charge()*B_IN_TESLA*RD_IN_METERS/fullEvent[ip].pT() );
	  double pxp = fullEvent[ip].pT()*cos( azimuth );
	  double pyp = fullEvent[ip].pT()*sin( azimuth );
	  double pzp = fullEvent[ip].pz();
	  particles.push_back( PseudoJet(pxp, pyp, pzp, Ep) );
	}
      }
      
      //Using the anti-kt algorithm:
      JetDefinition jet_def(antikt_algorithm, Rjet);
      
      ClusterSequence cs(particles, jet_def); 
      vector<PseudoJet> jets = cs.inclusive_jets();
      
      //Next, find the jet with highest transverse energy within the rapidity range:
      ET_gt = 0.; i_greatest_ET = 0;
      int jetsSize = jets.size();
      cout << "jetsSize = " << jetsSize << endl;
      for(int ij=0; ij<jetsSize; ij++){
	if(jets[ij].mperp() > ET_gt
	   && fabs(jets[ij].pseudorapidity()) < etaLeadingJetMax){
	  ET_gt = jets[ij].mperp();
	  i_greatest_ET = ij;
	  phiLeadingJet = jets[i_greatest_ET].phi();
	  etaLeadingJet = jets[i_greatest_ET].pseudorapidity();
	}
      }

      //Next, find the dijet of greatest energy on the opposite azimuth of the leading
      //jet within the rapidity range:
      ET_lt = 0.; i_greatest_opposite_ET = 0;
      for(int ij=0; ij<jetsSize; ij++){
	double delta_phi = max(phiLeadingJet, jets[ij].phi())
	  -min(phiLeadingJet, jets[ij].phi() );
	if(delta_phi > M_PI){delta_phi = 2.*M_PI-delta_phi;}
	if(jets[ij].mperp() > ET_lt
	   && delta_phi > (2./3.)*M_PI
	   && fabs(jets[ij].pseudorapidity()) < etaLeadingJetMax){
	  ET_lt = jets[ij].mperp();
	  i_greatest_opposite_ET = ij;
	}
      }
      
      cout << "ET_gt = " << ET_gt << ", ET_lt = " << ET_lt << endl;     

      //Now that the highest energy jet and the highest energy opposing jet has been determined, 
      //determine if they make it into the ATLAS analysis and bin them accordingly:

      if(ET_gt > ETgt_min && ET_lt > ETlt_min){
	
	//First of all, if this dijet triggers ATLAS, output ALL final hadrons for 
	//later tweaks of the analysis:
	for(int ip=0; ip<fullEventSize; ip++){
	  if(fullEvent[ip].isFinal() ){
	    foutEvents << fullEvent[ip].m() << " " << fullEvent[ip].charge() << " " << fullEvent[ip].pT()*cos(fullEvent[ip].phi() ) << " "
	           << fullEvent[ip].pT()*sin(fullEvent[ip].phi() ) << " " << fullEvent[ip].pz() << endl;
	  }
	}
	foutEvents << "EndOfEvent" << endl;


	//Put this event into the appropriate dN/dA bin:
	cout << " |E> - E<|/(E> + E<) = " << fabs(ET_gt-ET_lt)/(ET_gt+ET_lt) << endl;
	int proper_bin = (int)( ((double)nbins)*(fabs(ET_gt-ET_lt)/(ET_gt+ET_lt)) );
	//cout << "proper_bin = " << proper_bin << endl;
	if(proper_bin >=0 && proper_bin<nbins){
	  dNdA_hist[proper_bin] += (double)nbins/(double)numJets;
	}

	double Deltaphi = max(phiLeadingJet, jets[i_greatest_opposite_ET].phi() ) 
	  - min(phiLeadingJet, jets[i_greatest_opposite_ET].phi() );
	if(Deltaphi > M_PI){Deltaphi = 2.*M_PI-Deltaphi;}
	int phi_bin = (int)( ((double)nbins)*Deltaphi/M_PI);
	if(phi_bin >=0 && phi_bin<nbins){
	  dNdphi_hist[phi_bin] += (double)nbins/(double)numJets;
	}
	
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
