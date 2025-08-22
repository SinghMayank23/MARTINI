#include "MARTINI.h"
#include "fastjet/ClusterSequence.hh"
#include <sstream>
#include <math.h>
#define etaMax 4.5
#define ET_min 75.
#define ET_max 100.
#define etaJetMax 2.8
#define Rjet 0.4
#define pT_h_min 2.0
#define numJets 100
#define nbins 20
#define jTRange 5.
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

//  martini.readFile("setup_main_singlejets.dat");
  martini.readFile("setup_main_CMS_dNdA.dat");
  martini.init(argc, argv);
  martini.settings.listChanged();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  Event fullEvent;
  int numEvents;
  int jetCounter = 0; // Number of jets that make it into the analysis

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  //Energies of the jets:
  double ETJet;
  double phiJet;
  double etaJet;

  double dNdz_hist[nbins];        // The histogram for dN/dz
  double dNdjT_hist[nbins];        // The histogram for dN/dj_T
  int fullEventSize;

  //Make a tag for the output:
  string n_string = argv[1];
  n_string = "pp__"+n_string;

  //Initialize the histograms to be zero everywhere:
  for(int ih=0; ih<nbins; ih++){
    dNdz_hist[ih] = 0.;
    dNdjT_hist[ih] = 0.;
  }

  while(jetCounter<numJets)      // The main loop
    {
      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      martini.generateEvent(plist[0]); 
     
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
	fullEventSize = fullEvent.size();
      }
      else{
	fullEventSize = plist[0]->size();
      }

      //Now that the partons have undergone fragmentation, place all final, charged hadrons detected by 
      //ATLAS into a PseudoJet vector:
      vector<PseudoJet> particles;
      for(int ip=0; ip<fullEventSize; ip++){
	//Charged particles have a minimum p_t, neutral particles do not:
	if( AjOfHadrons){
	  if(fullEvent[ip].mT() > ET_min
	     //&& ( (fullEvent[ip].isCharged() && fullEvent[ip].mT() > ET_min) 
	     //	|| !(fullEvent[ip].isCharged()) )
	     && fabs(fullEvent[ip].eta()) < etaMax 
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
      
      //Now that the jets have been reconstructed, bin the hadrons differentially in 
      //j_T and z_T:
      int jetsSize = jets.size();
      for(int ij=0; ij<jetsSize; ij++){
	ETJet = jets[ij].mperp();
	etaJet = jets[ij].pseudorapidity();
	if( ETJet > ET_min && ETJet < ET_max
	    && fabs(etaJet) < etaJetMax ){
	  jetCounter += 1;
	  cout << "jetCounter = " << jetCounter << endl;
	  phiJet = jets[ij].phi();

	  if( AjOfHadrons ){
	    for(int ih=0; ih<fullEventSize; ih++){
	      if(fullEvent[ih].isFinal()
		 && fabs(fullEvent[ih].eta() ) < etaMax
		 && fullEvent[ih].pT() > pT_h_min){ //Brian Cole says that the hadrons have a minimum p_T for this binning. -8/16/2011
		double deltaEta = etaJet - fullEvent[ih].eta();
		double deltaPhi = max(phiJet, fullEvent[ih].phi() ) - min(phiJet, fullEvent[ih].phi() );
		if(deltaPhi > M_PI) deltaPhi = 2.*M_PI-deltaPhi ;
		
		double deltaR = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi );
		if(deltaR < Rjet){
		  double jT = fullEvent[ih].pT() * sin(deltaR);
		  int jTbin = (int)( jT/jTRange * (double)nbins );
		  if(jTbin < nbins){
		    dNdjT_hist[jTbin] += ((double)nbins/jTRange /(2.*M_PI*( (double)jTbin+0.5)*jTRange/((double)nbins) ) )/(double)numJets;
		  }
		}
	      }
	    }
	  }
	}
      }  
    }

  string jT_name = "dN2pijTdjT_"+n_string+"_80GeV.dat";
  fstream jTOut(jT_name.c_str(),ios::out); 
  for(int ib=0; ib<nbins; ib++)
    {
      jTOut << jTRange*(double)ib/(double)nbins << " " << dNdjT_hist[ib] << endl;
    }
  jTOut.close();
  
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
