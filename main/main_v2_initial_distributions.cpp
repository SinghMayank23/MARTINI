#include "MARTINI.h"
#include <sstream>
#define numJets 25000
#define nybins 200
#define DNDY2_MAX 10.
#define DNDY_MAX 10.
#define nptbins 50
//#define Y_MAX 0.35
#define Y_MAX 2.
#define PT_MAX 15.
#define R_MAX 10.
#define outputting true
#define nrapbins 200
#define ntau 200
//#define largeptbins 8

using namespace std;

int main(int argc, char* argv[]) 
{
  //The first bit of testing for heavy quarks!
  MARTINI martini;
  Event fullEvent;

  martini.readFile("setup_main_JPsi_v2.dat");
  if (argc > 1)
  {
    int job_num = stoi(argv[1]);
    martini.init(job_num);
  } else {
    martini.init(0);
  }
  martini.settings.listAll();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  string directory_name = "results/";
  string HQ_name = "LHC";

  int numEvents;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  double zmax = maxTime;

  Vec4 itsP;                          // Its P

  double** dN2piptdptdy_ccbar_y_beforeEvolution;
  double** dN2piptdptdy_ccbar_momrap_beforeEvolution;
  dN2piptdptdy_ccbar_y_beforeEvolution = new double* [nybins];
  dN2piptdptdy_ccbar_momrap_beforeEvolution = new double* [nybins];
  for(int iy=0; iy<nybins; iy++){
    dN2piptdptdy_ccbar_y_beforeEvolution[iy] = new double [nptbins];
    dN2piptdptdy_ccbar_momrap_beforeEvolution[iy] = new double [nptbins];
    for(int ip=0; ip<nptbins; ip++){
      dN2piptdptdy_ccbar_y_beforeEvolution[iy][ip] = 0.;
      dN2piptdptdy_ccbar_momrap_beforeEvolution[iy][ip] = 0.;
    }
  }

  double dNdy_ccbar_beforeEvolution[nybins];
  double dNdmomrap_ccbar_beforeEvolution[nybins];

  for(int iy=0; iy<nybins; iy++){
    dNdy_ccbar_beforeEvolution[iy] = 0.;
    dNdmomrap_ccbar_beforeEvolution[iy] = 0.;
  }

  cout << numJets << endl;

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      cout << "event_counter  " << event_counter << endl;
      martini.generateEventHeavyQuarks_w_colllist_OR_IPGfile(plist[0]); 

      int eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//Pseudorapidity is more appropriate here:
	double y      = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double momrap = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdy_ccbar_beforeEvolution[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}
	int imomrap = (int)( (momrap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(imomrap>=0 && imomrap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdmomrap_ccbar_beforeEvolution[imomrap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	int ipt = (int)((pt/PT_MAX)*nptbins );

	if(iy>=0 && iy<nybins){
	if(abs(plist[0]->at(ie).id())==4){
	  dN2piptdptdy_ccbar_y_beforeEvolution[iy][ipt] += (PT_MAX/(double)nptbins)*(DNDY_MAX*2./(double)nybins)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	}
	}
	if(imomrap>=0 && imomrap<nybins){
	if(abs(plist[0]->at(ie).id())==4){
	  dN2piptdptdy_ccbar_momrap_beforeEvolution[imomrap][ipt] += (PT_MAX/(double)nptbins)*(DNDY_MAX*2./(double)nybins)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	}
	}

      }
    }


  stringstream s;
  s << argv[1];
  
  if(outputting){

    string ccbarBeforeEvolutionname = directory_name+"dN2piptdpt_ccbar_y_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOut;
    ccbarBeforeEvolutionOut.open(ccbarBeforeEvolutionname.c_str());

    string bbbarBeforeEvolutionname = directory_name+"dN2piptdpt_ccbar_momrap_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOut;
    bbbarBeforeEvolutionOut.open(bbbarBeforeEvolutionname.c_str());

//    string ccbarBeforeHydroname = directory_name+"dN2piptdpt_ccbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydroOut;
//    ccbarBeforeHydroOut.open(ccbarBeforeHydroname.c_str());
//
//    string bbbarBeforeHydroname = directory_name+"dN2piptdpt_bbbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream bbbarBeforeHydroOut;
//    bbbarBeforeHydroOut.open(bbbarBeforeHydroname.c_str());
//
//    string ccbarBeforeEvolutionLargeptname = directory_name+"dN2piptdpt_ccbar_beforeEvolution_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeEvolutionLargeptOut;
//    ccbarBeforeEvolutionLargeptOut.open(ccbarBeforeEvolutionLargeptname.c_str());

//    string ccbarBeforeHydroLargeptname = directory_name+"dN2piptdpt_ccbar_beforeHydro_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydroLargeptOut;
//    ccbarBeforeHydroLargeptOut.open(ccbarBeforeHydroLargeptname.c_str());
//
//    string ccbarLargeptname = directory_name+"dN2piptdpt_ccbar_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarLargeptOut;
//    ccbarLargeptOut.open(ccbarLargeptname.c_str());
    for(int iy=0; iy<nybins; iy++){
    for(int ip=0; ip<nptbins; ip++){
      ccbarBeforeEvolutionOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_y_beforeEvolution[iy][ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_momrap_beforeEvolution[iy][ip]/(2.*event_counter) << endl;
    }
      ccbarBeforeEvolutionOut << endl; 
      bbbarBeforeEvolutionOut << endl; 
    }
  
    string ccbarBeforeEvolutiondNdyname = directory_name+"dNdy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdyOut;
    ccbarBeforeEvolutiondNdyOut.open(ccbarBeforeEvolutiondNdyname.c_str());

    string ccbarBeforeEvolutiondNdmomrapname = directory_name+"dNdmomrap_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdmomrapOut;
    ccbarBeforeEvolutiondNdmomrapOut.open(ccbarBeforeEvolutiondNdmomrapname.c_str());

    for(int iy=0; iy<nybins; iy++){
      ccbarBeforeEvolutiondNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
      ccbarBeforeEvolutiondNdmomrapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdmomrap_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
    }


    ccbarBeforeEvolutionOut.close();

    bbbarBeforeEvolutionOut.close();

    ccbarBeforeEvolutiondNdyOut.close();
    ccbarBeforeEvolutiondNdmomrapOut.close();
  }

  for(int iy = 0; iy < nybins; iy++){
    delete dN2piptdptdy_ccbar_y_beforeEvolution[iy];
    delete dN2piptdptdy_ccbar_momrap_beforeEvolution[iy];
  }

  delete dN2piptdptdy_ccbar_y_beforeEvolution;
  delete dN2piptdptdy_ccbar_momrap_beforeEvolution;

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
