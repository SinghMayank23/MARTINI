#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#define ET_min 0.1
#define etaMax 2.
#define ETgt_min 120.
#define etaLeadingJetMax 2.
#define ETlt_min 50.
#define Rjet 0.5
#define B_IN_TESLA 3.8
#define RD_IN_METERS 1.
#define nEvent 10000
#define nbins 20

using namespace fastjet;
using namespace Pythia8; 

int main() {

  Pythia pythia;
  pythia.readFile("setup_main_pp_LHC.dat");
  pythia.readString("PhaseSpace:pTHatMin = 80.");

  double jetpTminIn;
  jetpTminIn=pythia.parm("PhaseSpace:pTHatMin");
  pythia.setJetpTmin(jetpTminIn);

  pythia.init();

  int ievent=0;
  bool flag;

  //Energies of the jets:
  double ET_gt, ET_lt;
  int i_greatest_ET, i_greatest_opposite_ET;
  double phiLeadingJet;
  double etaLeadingJet;

  double dNdA_hist[nbins];          // The histogram for dN/dA
  double dNdphi_hist[nbins];        // The histogram for dN/dphi

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
  }

  //Make a tag for the output:
  string time_name = "elapse_dijet.dat";

  //creating a file for events
  fstream fouttime( time_name.c_str() ,ios::out | ios::app);
  fouttime << "ievent : " << ievent << "	" << asctime(timeinfo) << endl;

  // The main loop
  while(ievent<nEvent)
  {
    flag=pythia.next();
    if(flag==true){
  
      //Now that the partons have undergone fragmentation, place all final,
      //charged hadrons detected by ATLAS into a PseudoJet vector:

      vector<PseudoJet> particles;
      for (int i=0; i<pythia.event.size(); ++i){
    	//Charged particles have a minimum p_t, neutral particles do not:
  if (pythia.event[i].mT() > ET_min
	  	&&fabs(pythia.event[i].eta()) < etaMax
			&&pythia.event[i].isFinal()
      &&pythia.event[i].pT() > 0.149896*fabs(pythia.event[i].charge())*B_IN_TESLA*RD_IN_METERS ){
	  double Ep = pythia.event[i].eCalc();
	  double azimuth = pythia.event[i].phi()+asin( 0.149896*pythia.event[i].charge()*B_IN_TESLA*RD_IN_METERS/pythia.event[i].pT() );
	  double pxp = pythia.event[i].pT()*cos( azimuth );
	  double pyp = pythia.event[i].pT()*sin( azimuth );
	  double pzp = pythia.event[i].pz();
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
	double delta_phi = max(phiLeadingJet, jets[ij].phi()) - min(phiLeadingJet, jets[ij].phi() );
	if(delta_phi > M_PI){delta_phi = 2.*M_PI-delta_phi;}
	if(jets[ij].mperp() > ET_lt
	   && delta_phi > (2./3.)*M_PI
	   && fabs(jets[ij].pseudorapidity()) < etaLeadingJetMax){
	  ET_lt = jets[ij].mperp();
	  i_greatest_opposite_ET = ij;
	}
      }
 
      //Now that the highest energy jet and the highest energy opposing jet has been determined, 
      //determine if they make it into the ATLAS analysis and bin them accordingly:

      if(ET_gt > ETgt_min && ET_lt > ETlt_min){

	//Put this event into the appropriate dN/dA bin:
  cout << "ET_gt = " << ET_gt << ", ET_lt = " << ET_lt
	     << " |E> - E<|/(E> + E<) = " << fabs(ET_gt-ET_lt)/(ET_gt+ET_lt) << endl;
	int proper_bin = (int)( ((double)nbins)*(fabs(ET_gt-ET_lt)/(ET_gt+ET_lt)) );
	//cout << "proper_bin = " << proper_bin << endl;
	if(proper_bin >=0 && proper_bin<nbins){
	  dNdA_hist[proper_bin] += (double)nbins/(double)nEvent;
	}

	double Deltaphi = max(phiLeadingJet, jets[i_greatest_opposite_ET].phi() ) 
	  - min(phiLeadingJet, jets[i_greatest_opposite_ET].phi() );
	if(Deltaphi > M_PI){Deltaphi = 2.*M_PI-Deltaphi;}
	int phi_bin = (int)( ((double)nbins)*Deltaphi/M_PI);
	if(phi_bin >=0 && phi_bin<nbins){
	  dNdphi_hist[phi_bin] += (double)nbins/(double)nEvent;
	}
	
	//Finally, add to the counter:
	ievent += 1;
  cout << "ievent = " << ievent << endl;

  if (ievent % 100 == 0) {
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    fouttime << "ievent : " << ievent << "	" << asctime(timeinfo) << endl;
  }
      }  

    }
  }


  string pion_name = "dNdA_pp_CMS_dijet_80GeV.dat";
  fstream foutdNdA( pion_name.c_str() ,ios::out);   
  string phi_name = "dNdphi_pp_CMS_dijet_80GeV.dat";
  fstream foutdNdphi(phi_name.c_str(), ios::out);

  for(int ib=0; ib<nbins; ib++)
    {
      foutdNdA << (double)ib/(double)nbins << " " << dNdA_hist[ib] << endl;
      foutdNdphi << M_PI*(double)ib/(double)nbins << " " << dNdphi_hist[ib] << endl;
    }
  foutdNdA.close();
  foutdNdphi.close();
  fouttime.close();
  remove( time_name.c_str() );

  cout << "\nnEvent = " << nEvent << endl; 


  return 0;
}  
