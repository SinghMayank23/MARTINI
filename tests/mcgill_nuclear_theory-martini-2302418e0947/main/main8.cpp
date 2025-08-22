#include "MARTINI.h"
#include "fastjet/ClusterSequence.hh"
#include <sstream>
#include <math.h>
#define ET_min 0.7
#define etaMax 4.5
#define ETgt_min 100.0
#define etaLeadingJetMax 2.8
#define ETlt_min 25.
#define Rjet 0.4
#define pT_h_min 2.0
#define numJets 10
#define nbins 20
#define jTRange 5.0
#define zRange 1.0
#define AjOfHadrons true

using namespace fastjet;
using namespace std;

//We are now using the anti-kt algorithm that most, if not all, of the 
//LHC collaborations are using:
int main(int argc, char* argv[]) 
{
  // create MARTINI object - it will automatically generate PYTHIA object and set its standard values
  // those values can be changed below before doing the explicit PYTHIA init.
  MARTINI martini;

  martini.readFile("setup_files/setup_main8.dat");
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
  double phiLeadingJet, phiSubleadingJet;
  double etaLeadingJet, etaSubleadingJet;

  double dNdA_hist[nbins];          // The histogram for dN/dA
  double dNdphi_hist[nbins];        // The histogram for dN/dphi
  double dNdjT_leadingJet_hist[nbins];        // The histogram for dN/dj_T of the leading jet
  double dNdjT_subleadingJet_hist[nbins];        // The histogram for dN/dj_T of the sub-leading jet
  double dNdz_leadingJet_hist[nbins];        // The histogram for dN/dz of the leading jet
  double dNdz_subleadingJet_hist[nbins];        // The histogram for dN/dz of the sub-leading jet
  int fullEventSize;

  //Make a tag for the output:
  string n_string = argv[1];
  n_string = "PbPb_276TeV_alpha27_BFieldOn_pCut1p5_"+n_string;
  //n_string = "PbPb_CYbestattempt_pthatmin25_alpha27_maxTime20_"+n_string;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdA_hist[ih] = 0.;
    dNdphi_hist[ih] = 0.;
    dNdjT_leadingJet_hist[ih] = 0.;
    dNdjT_subleadingJet_hist[ih] = 0.;
    dNdz_leadingJet_hist[ih] = 0.;
    dNdz_subleadingJet_hist[ih] = 0.;
  }

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

      if( AjOfHadrons ){
	// fragmentation
	fullEvent.clear();
	martini.pythia.event.clear();
	numEvents+=martini.fragmentation( plist );
	fullEvent = martini.pythia.event;
      }

      //Now that the partons have undergone fragmentation, place all final, charged hadrons detected by 
      //ATLAS into a PseudoJet vector:
      if( AjOfHadrons){
	fullEventSize = fullEvent.size();
      }
      else{
	fullEventSize = plist[0]->size();
      }
      vector<PseudoJet> particles;
      for(int ip=0; ip<fullEventSize; ip++){
	//WARNING: THE B-FIELD IS OFF!
	//Charged particles have a minimum p_t, neutral particles do not:
	if( AjOfHadrons){
	  if(fabs(fullEvent[ip].eta()) < etaMax 
	     && ( (fullEvent[ip].isCharged() //&& fullEvent[ip].mT() > ET_min) 
		   && fullEvent[ip].pT() > 0.3597*fabs(fullEvent[ip].charge()) )
		  || !(fullEvent[ip].isCharged()) )
	     ////&& fullEvent[ip].mT() > ET_min
	     && fullEvent[ip].isFinal() ){
	    double Ep = fullEvent[ip].eCalc();
	    double azimuth = fullEvent[ip].phi()+asin( 0.3597*fullEvent[ip].charge()/fullEvent[ip].pT() );
	    double pxp = fullEvent[ip].pT()*cos( azimuth );
	    double pyp = fullEvent[ip].pT()*sin( azimuth );
	    double pzp = fullEvent[ip].pz();
	    particles.push_back( PseudoJet(pxp, pyp, pzp, Ep) );
	  }
	}
	else{
	  double pp = sqrt( plist[0]->at(ip).p().px()*plist[0]->at(ip).p().px()
			    + plist[0]->at(ip).p().py()*plist[0]->at(ip).p().py()
			    + plist[0]->at(ip).p().pz()*plist[0]->at(ip).p().pz() );
	  double etap = 0.5*log((pp + plist[0]->at(ip).p().pz() )/(pp - plist[0]->at(ip).p().pz() ) ) ;
	  if(fabs(etap) < etaMax ){
	    double Ep = sqrt( plist[0]->at(ip).mass()*plist[0]->at(ip).mass() + pp*pp );
	    particles.push_back( PseudoJet(plist[0]->at(ip).p().px(), 
					   plist[0]->at(ip).p().py(),
					   plist[0]->at(ip).p().pz(), Ep ) );
	  }
	}
      }
      
      //Using the anti-kt algorithm:
      JetDefinition jet_def(antikt_algorithm, Rjet);
      
      ClusterSequence cs(particles, jet_def); 
      vector<PseudoJet> jets = cs.inclusive_jets();
      
      //Next, find the jet with highest transverse energy:
      ET_gt = 0.; i_greatest_ET = 0;
      int jetsSize = jets.size();
      cout << "jetsSize = " << jetsSize << endl;
      for(int ij=0; ij<jetsSize; ij++){
	if(jets[ij].mperp() > ET_gt){
	  ET_gt = jets[ij].mperp();
	  i_greatest_ET = ij;
	  phiLeadingJet = jets[i_greatest_ET].phi();
	  etaLeadingJet = jets[i_greatest_ET].pseudorapidity();
	}
      }
      
      //Next, find the dijet of greatest energy on the opposite azimuth of the leading
      //jet:
      ET_lt = 0.; i_greatest_opposite_ET = 0;
      for(int ij=0; ij<jetsSize; ij++){

	double delta_phi = max(phiLeadingJet, jets[ij].phi())
	  -min(phiLeadingJet, jets[ij].phi() );
	if(delta_phi > M_PI){delta_phi = 2.*M_PI-delta_phi;}

	if(jets[ij].mperp() > ET_lt
	   && delta_phi > 0.5*M_PI){
	  ET_lt = jets[ij].mperp();
	  i_greatest_opposite_ET = ij;
	  phiSubleadingJet = jets[ij].phi();
	  etaSubleadingJet = jets[ij].pseudorapidity();

	}
      }
      
      cout << "ET_gt = " << ET_gt << ", ET_lt = " << ET_lt << endl;     
      
      //Now that the highest energy jet and the highest energy opposing jet has been determined, 
      //determine if they make it into the ATLAS analysis and bin them accordingly:

      if(ET_gt > ETgt_min && ET_lt > ETlt_min && fabs(etaLeadingJet) < etaLeadingJetMax){
	
	cout << "etaLeadingJet = " << etaLeadingJet << endl;
	cout << "etaSubleadingJet = " << etaSubleadingJet << endl;

	//First of all, if this dijet triggers ATLAS, output ALL final hadrons for 
	//later tweaks of the analysis:
	for(int ip=0; ip<fullEventSize; ip++){
	  if(fabs(fullEvent[ip].eta()) < etaMax
	     && fullEvent[ip].mT() > ET_min
	     && fullEvent[ip].isFinal() ){
	    //foutEvents << fullEvent[ip].m() << " " << fullEvent[ip].charge() << " " << fullEvent[ip].pT()*cos(fullEvent[ip].phi() ) << " "
	    //       << fullEvent[ip].pT()*sin(fullEvent[ip].phi() ) << " " << fullEvent[ip].pz() << endl;
	  }
	}
	//foutEvents << "EndOfEvent" << endl;
	
	
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

	//Bin the appropriate hadrons into the dN/dj_T histograms:
	if( AjOfHadrons ){
	  for(int ih=0; ih<fullEventSize; ih++){
	    if(fullEvent[ih].isFinal()
	       && fabs(fullEvent[ih].eta() ) < etaMax
	       && fullEvent[ih].pT() > pT_h_min){ //Brian Cole says that the hadrons have a minimum p_T for this binning. -8/16/2011                                         

	      double deltaEtaLeading = etaLeadingJet - fullEvent[ih].eta();
	      double deltaPhiLeading = max(phiLeadingJet, fullEvent[ih].phi() ) - min(phiLeadingJet, fullEvent[ih].phi() );
	      if(deltaPhiLeading > M_PI) deltaPhiLeading = 2.*M_PI-deltaPhiLeading ;
	      double deltaRLeading = sqrt(deltaEtaLeading*deltaEtaLeading+deltaPhiLeading*deltaPhiLeading );

	      if(deltaRLeading < Rjet){
		double jT = fullEvent[ih].pT() * sin(deltaRLeading);
		int jTbin = (int)( jT/jTRange * (double)nbins );
		double zh = fullEvent[ih].pT() * cos(deltaRLeading)/ET_gt;
		int zbin = (int)( zh/zRange *(double)nbins );
		if(jTbin < nbins){
		  dNdjT_leadingJet_hist[jTbin] += ((double)nbins/jTRange /(2.*M_PI*( (double)jTbin+0.5)*jTRange/((double)nbins) ) )/(double)numJets;
		}
		if(zbin < nbins){
		  dNdz_leadingJet_hist[zbin] += ((double)nbins/zRange )/(double)numJets;
		}
	      }

	      double deltaEtaSubleading = etaSubleadingJet - fullEvent[ih].eta();
	      double deltaPhiSubleading = max(phiSubleadingJet, fullEvent[ih].phi() ) - min(phiSubleadingJet, fullEvent[ih].phi() );
	      if(deltaPhiSubleading > M_PI) deltaPhiSubleading = 2.*M_PI-deltaPhiSubleading ;
	      double deltaRSubleading = sqrt(deltaEtaSubleading*deltaEtaSubleading+deltaPhiSubleading*deltaPhiSubleading );

	      if(deltaRSubleading < Rjet){
		double jT = fullEvent[ih].pT() * sin(deltaRSubleading);
		double zh = fullEvent[ih].pT() * cos(deltaRSubleading)/ET_lt;
		int jTbin = (int)( jT/jTRange * (double)nbins );
		int zbin = (int)(zh/zRange *(double)nbins );
		if(jTbin < nbins){
		  dNdjT_subleadingJet_hist[jTbin] += ((double)nbins/jTRange /(2.*M_PI*( (double)jTbin+0.5)*jTRange/((double)nbins) ) )/(double)numJets ;
		}
		if(zbin < nbins){
		  dNdz_subleadingJet_hist[zbin] += ((double)nbins/zRange )/(double)numJets ;
		}
	      }
	    }
	  }
	}
	
	//Finally, add to the counter:
	event_counter += 1;
      }
    }
  
  string A_name = "dNdA_"+n_string+".dat";
  fstream foutdNdA(A_name.c_str(),ios::out);
  string phi_name = "dNdphi_"+n_string+".dat";
  fstream foutdNdphi(phi_name.c_str(), ios::out);
  string jT_leadingJet_name = "dN2pijTdjT_leadingJet_"+n_string+".dat";
  fstream jT_leadingJetOut(jT_leadingJet_name.c_str(), ios::out);
  string jT_subleadingJet_name = "dN2pijTdjT_subleadingJet_"+n_string+".dat";
  fstream jT_subleadingJetOut(jT_subleadingJet_name.c_str(), ios::out);
  string z_leadingJet_name = "dNdz_leadingJet_"+n_string+".dat";
  fstream z_leadingJetOut(z_leadingJet_name.c_str(), ios::out);
  string z_subleadingJet_name = "dNdz_subleadingJet_"+n_string+".dat";
  fstream z_subleadingJetOut(z_subleadingJet_name.c_str(), ios::out);
  for(int ib=0; ib<nbins; ib++)
    {
      foutdNdA << (double)ib/(double)nbins << " " << dNdA_hist[ib] << endl;
      foutdNdphi << M_PI*(double)ib/(double)nbins << " " << dNdphi_hist[ib] << endl;
      jT_leadingJetOut << jTRange * (double)ib/(double)nbins << " " << dNdjT_leadingJet_hist[ib] << endl;
      jT_subleadingJetOut << jTRange * (double)ib/(double)nbins << " " << dNdjT_subleadingJet_hist[ib] << endl;
      z_leadingJetOut << zRange * (double)ib/(double)nbins << " " << dNdz_leadingJet_hist[ib] << endl;
      z_subleadingJetOut << zRange * (double)ib/(double)nbins << " " << dNdz_subleadingJet_hist[ib] << endl;
    }
  foutdNdA.close();
  foutdNdphi.close();
  jT_leadingJetOut.close();
  jT_subleadingJetOut.close();
  z_leadingJetOut.close();
  z_subleadingJetOut.close();
  
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
