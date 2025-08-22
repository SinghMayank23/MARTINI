#include "MARTINI.h"
#include <sstream>
#define numJets 50000
#define nybins 100
#define DNDY_MAX 5.
#define nptbins 50
//#define Y_MAX 0.35
#define Y_MAX 2.
#define PT_MAX 5.
#define R_MAX 10.
#define outputting true

using namespace std;

int main(int argc, char* argv[]) 
{
  //The first bit of testing for heavy quarks!
  MARTINI martini;
  Event fullEvent;

  martini.readFile("setup_files/setup_main_JPsi_v2.dat");
  martini.init(argc, argv);
  martini.settings.listAll();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  string directory_name = "HQTestingAndResults/LHC/";
  string HQ_name = "LHC";

  int numEvents;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  Vec4 itsP;                          // Its P

  //Define the histograms:
  double NJPsi_t[mt/10];
  double NJPsi_t_rec[mt/10];
  cout << "mt/10 = " << mt/10 << endl;

  double dN2piptdptdy_ccbar_beforeEvolution[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution[nptbins];
  
  double dN2piptdptdy_ccbar[nptbins];
  double dN2piptdptdy_bbbar[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar[nptbins];

  double dN2piptdptdy_DDbar[nptbins];
  double dN2piptdptdy_DBDBbar[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar[nptbins];
  double cos2Phi_dN2piptdptdy_DBDBbar[nptbins];

  double dN2piptdptdy_JPsi[nptbins];
  double dN2piptdptdy_ExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_JPsi[nptbins];
  double cos2Phi_dN2piptdptdy_ExcitedJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];

  double dN2piptdptdy_Upsilon[nptbins];
  double dN2piptdptdy_ExcitedUpsilon[nptbins];

  double dN2piptdptdy_bcbar[nptbins];
  double dN2piptdptdy_Excitedbcbar[nptbins];

  //Counting the number of bound states:
  double NJPsi, NExcitedJPsi, NUpsilon, NExcitedUpsilon, NBc, NExcitedBc;
  NJPsi = NExcitedJPsi = 0.;
  NUpsilon = NExcitedUpsilon = 0.;
  NBc = NExcitedBc = 0.;

  for(int ip=0; ip<nptbins; ip++){
    dN2piptdptdy_ccbar_beforeEvolution[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution[ip] = 0.;

    dN2piptdptdy_ccbar[ip] = 0.;
    dN2piptdptdy_bbbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar[ip] = 0.;

    dN2piptdptdy_DDbar[ip] = 0.;
    dN2piptdptdy_DBDBbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;

    dN2piptdptdy_JPsi[ip] = 0.;
    dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_JPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;

    dN2piptdptdy_Upsilon[ip] = 0.;
    dN2piptdptdy_ExcitedUpsilon[ip] = 0.;

    dN2piptdptdy_bcbar[ip] = 0.;
    dN2piptdptdy_Excitedbcbar[ip] = 0.;
  }

  double dNdy_ccbar_beforeEvolution[nybins];
  double dNdy_ccbar[nybins];
  double dNdy_DDbar[nybins];
  double dNdy_JPsi[nybins];
  double dNdy_ExcitedJPsi[nybins];

  for(int iy=0; iy<nybins; iy++){
    dNdy_ccbar_beforeEvolution[iy] = 0.;
    dNdy_ccbar[iy] = 0.;
    dNdy_DDbar[iy] = 0.;
    dNdy_JPsi[iy] = 0.;
    dNdy_ExcitedJPsi[iy] = 0.;
  }

  //double dN4pir2dr_ccbar_beforeEvolution[nptbins];
  //double dN4pir2dr_ccbar[nptbins];
  //for(int ir=0; ir<nptbins; ir++){
  //dN4pir2dr_ccbar_beforeEvolution[ir] = 0.;
  //dN4pir2dr_ccbar[ir] = 0.;
  //}

  ////Now, we test the fragmentation of the heavy quarks with 
  ////generateTestEvent_HQ:
  //ofstream NboundOut;
  //NboundOut.open("Nbound.dat");
  //for(int il=0; il<50; il++){
  //double Nbound = 0.;
  //int NSAMP = 10000;
  //double L = 0.7+0.5*(double)il/50.;
  //cout << "L = " << L << endl;
  //for(int isamp=0; isamp<NSAMP; isamp++){
  //  plist[0]->clear();
  //  plist[1]->clear();
  //  martini.generateTestEvent_HQ(plist[0], L);
  //  martini.fragmentation_HQ(plist);
  //  int eventSize = martini.pythia.event.size();
  //  //cout << "Event size in the fragmentation testing = " << eventSize << endl;
    //  for(int ie=0; ie<eventSize; ie++){
  //
  ////cout << "martini.pythia.event[" << ie << "].id() = " << martini.pythia.event[ie].id() << endl;
  //if(martini.pythia.event[ie].id() == 443){
  //  Nbound += 1./(double)NSAMP;
  //}
  //if(martini.pythia.event[ie].id() == 443000){
  //  Nbound += 1./(double)NSAMP;
  //}
  //  }
  //}
  //NboundOut << L << " " << Nbound << endl;
  //}
  //NboundOut.close();

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      martini.generateEventHeavyQuarks(plist[0]); 

      int eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//Pseudorapidity is more appropriate here:
	double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(martini.pythia.event[ie].id())==4){
	    dNdy_ccbar_beforeEvolution[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(martini.pythia.event[ie].id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	}

      }
	  
      int counter = 0;
      if (martini.returnEvolution() == 1)
	{
	  for(int i=0; i<mt; i++)
	    {
	      //counter = martini.evolveHeavyQuarks(plist, counter, i);
	      counter = martini.evolveAndHadronizeHeavyQuarks(plist, counter, i);
	      counter+=1;
	    }
	}

      eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//Pseudorapidity is more appropriate here:
	double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double cos2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	}

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(martini.pythia.event[ie].id())==4){
	    dNdy_ccbar[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==4){
	    dN2piptdptdy_ccbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(martini.pythia.event[ie].id())==5){
	    dN2piptdptdy_bbbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	}
      }
	  
      eventSize = martini.pythia.event.size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = martini.pythia.event[ie].p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = itsP.mCalc();
	double pz = itsP.pz();
	//Pseudorapidity is more appropriate here:
	double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double cos2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	}

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy >= 0 && iy < nybins){
	  if(abs(martini.pythia.event[ie].id())==411){
	    dNdy_DDbar[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dNdy_JPsi[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dNdy_ExcitedJPsi[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(martini.pythia.event[ie].id())==511){
	    dN2piptdptdy_DBDBbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DBDBbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dN2piptdptdy_JPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_JPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    }
	    NJPsi += 1.;
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dN2piptdptdy_ExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedJPsi += 1.;
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    }
	  }
	  if(martini.pythia.event[ie].id()==553){
	    dN2piptdptdy_Upsilon[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NUpsilon += 1.;
	  }
	  if(martini.pythia.event[ie].id()==553000){
	    dN2piptdptdy_ExcitedUpsilon[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedUpsilon += 1.;
	  }
	  if(martini.pythia.event[ie].id()==543){
	    dN2piptdptdy_bcbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NBc += 1.;
	  }
	  if(martini.pythia.event[ie].id()==543000){
	    dN2piptdptdy_Excitedbcbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NBc += 1.;
	  }
	}
      }
    }

  //Output the properly normalized numbers of quarkonia:
  NJPsi /= (double)numJets;
  NExcitedJPsi /= (double)numJets;
  NUpsilon /= (double)numJets;
  NExcitedUpsilon /= (double)numJets;
  NBc /= (double)numJets;
  NExcitedBc /= (double)numJets;

  cout << "The average number of JPsi particles per collision = " << NJPsi << endl;
  cout << "The average number of excited charmonia per collision = " << NExcitedJPsi << endl;
  cout << "The average number of Upsilon particles per collision = " << NUpsilon << endl;
  cout << "The average number of excited bottomonia per collision = " << NExcitedUpsilon << endl;
  cout << "The average number of B_c mesons per collision = " << NBc << endl;
  cout << "The average number of excited B_c mesons per collision = " << NExcitedBc << endl;

  stringstream s;
  s << argv[1];
  
  if(outputting){

    string ccbarBeforeEvolutionname = directory_name+"dN2piptdpt_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOut;
    ccbarBeforeEvolutionOut.open(ccbarBeforeEvolutionname.c_str());

    string bbbarBeforeEvolutionname = directory_name+"dN2piptdpt_bbbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOut;
    bbbarBeforeEvolutionOut.open(bbbarBeforeEvolutionname.c_str());

    string ccbarname = directory_name+"dN2piptdpt_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOut;
    ccbarOut.open(ccbarname.c_str());
    string cos2Phi_ccbarname = directory_name+"cos2Phi_dN2piptdpt_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOut;
    cos2Phi_ccbarOut.open(cos2Phi_ccbarname.c_str());

    string bbbarname = directory_name+"dN2piptdpt_bbbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOut;
    bbbarOut.open(bbbarname.c_str());
    string cos2Phi_bbbarname = directory_name+"cos2Phi_dN2piptdpt_bbbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOut;
    cos2Phi_bbbarOut.open(cos2Phi_bbbarname.c_str());

    string DDbarname = directory_name+"dN2piptdpt_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOut;
    DDbarOut.open(DDbarname.c_str());
    string cos2Phi_DDbarname = directory_name+"cos2Phi_dN2piptdpt_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOut;
    cos2Phi_DDbarOut.open(cos2Phi_DDbarname.c_str());

    string DBDBbarname = directory_name+"dN2piptdpt_DBDBbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DBDBbarOut;
    DBDBbarOut.open(DBDBbarname.c_str());
    string cos2Phi_DBDBbarname = directory_name+"cos2Phi_dN2piptdpt_DBDBbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DBDBbarOut;
    cos2Phi_DBDBbarOut.open(cos2Phi_DBDBbarname.c_str());

    string JPsiname = directory_name+"dN2piptdpt_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream JPsiOut;
    JPsiOut.open(JPsiname.c_str());
    string cos2Phi_JPsiname = directory_name+"cos2Phi_dN2piptdpt_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_JPsiOut;
    cos2Phi_JPsiOut.open(cos2Phi_JPsiname.c_str());

    string OnlyDirectJPsiname = directory_name+"dN2piptdpt_OnlyDirectJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream OnlyDirectJPsiOut;
    OnlyDirectJPsiOut.open(OnlyDirectJPsiname.c_str());
    string cos2Phi_OnlyDirectJPsiname = directory_name+"cos2Phi_dN2piptdpt_OnlyDirectJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_OnlyDirectJPsiOut;
    cos2Phi_OnlyDirectJPsiOut.open(cos2Phi_OnlyDirectJPsiname.c_str());

    string ExcitedJPsiname = directory_name+"dN2piptdpt_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedJPsiOut;
    ExcitedJPsiOut.open(ExcitedJPsiname.c_str());
    string cos2Phi_ExcitedJPsiname = directory_name+"cos2Phi_dN2piptdpt_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ExcitedJPsiOut;
    cos2Phi_ExcitedJPsiOut.open(cos2Phi_ExcitedJPsiname.c_str());

    string OnlyDirectExcitedJPsiname = directory_name+"dN2piptdpt_OnlyDirectExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream OnlyDirectExcitedJPsiOut;
    OnlyDirectExcitedJPsiOut.open(OnlyDirectExcitedJPsiname.c_str());
    string cos2Phi_OnlyDirectExcitedJPsiname = directory_name+"cos2Phi_dN2piptdpt_OnlyDirectExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_OnlyDirectExcitedJPsiOut;
    cos2Phi_OnlyDirectExcitedJPsiOut.open(cos2Phi_OnlyDirectExcitedJPsiname.c_str());

    string Upsilonname = directory_name+"dN2piptdpt_Upsilon_"+HQ_name+"_"+s.str()+".dat";
    ofstream UpsilonOut;
    UpsilonOut.open(Upsilonname.c_str());

    string ExcitedUpsilonname = directory_name+"dN2piptdpt_ExcitedUpsilon_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedUpsilonOut;
    ExcitedUpsilonOut.open(ExcitedUpsilonname.c_str());

    string bcbarname = directory_name+"dN2piptdpt_bcbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream bcbarOut;
    bcbarOut.open(bcbarname.c_str());

    string Excitedbcbarname = directory_name+"dN2piptdpt_Excitedbcbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedbcbarOut;
    ExcitedbcbarOut.open(Excitedbcbarname.c_str());
    
    for(int ip=0; ip<nptbins; ip++){
      ccbarBeforeEvolutionOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution[ip]/(2.*event_counter) << endl;

      ccbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar[ip]/(2.*event_counter) << endl;
      bbbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar[ip]/(2.*event_counter) << endl;

      DDbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar[ip]/(2.*event_counter) << endl;
      DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar[ip]/(2.*event_counter) << endl;
      cos2Phi_DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;

      JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      cos2Phi_ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;

      UpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Upsilon[ip]/(double)event_counter << endl;
      ExcitedUpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedUpsilon[ip]/(double)event_counter << endl;

      bcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bcbar[ip]/(double)event_counter << endl;
      ExcitedbcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Excitedbcbar[ip]/(double)event_counter << endl;
    }
    
    string ccbarBeforeEvolutiondNdyname = directory_name+"dNdy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdyOut;
    ccbarBeforeEvolutiondNdyOut.open(ccbarBeforeEvolutiondNdyname.c_str());

    string ccbardNdyname = directory_name+"dNdy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdyOut;
    ccbardNdyOut.open(ccbardNdyname.c_str());

    string DDbardNdyname = directory_name+"dNdy_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbardNdyOut;
    DDbardNdyOut.open(DDbardNdyname.c_str());

    string JPsidNdyname = directory_name+"dNdy_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream JPsidNdyOut;
    JPsidNdyOut.open(JPsidNdyname.c_str());

    string ExcitedJPsidNdyname = directory_name+"dNdy_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedJPsidNdyOut;
    ExcitedJPsidNdyOut.open(ExcitedJPsidNdyname.c_str());
    
    for(int iy=0; iy<nybins; iy++){
      ccbarBeforeEvolutiondNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
      ccbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar[iy]/(2.*(double)event_counter) << endl;
      DDbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_DDbar[iy]/(2.*(double)event_counter) << endl;
      JPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_JPsi[iy]/(2.*(double)event_counter) << endl;
      ExcitedJPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ExcitedJPsi[iy]/(2.*(double)event_counter) << endl;
    }
    
    ccbarBeforeEvolutionOut.close();

    bbbarBeforeEvolutionOut.close();

    ccbarOut.close();
    cos2Phi_ccbarOut.close();

    bbbarOut.close();
    cos2Phi_bbbarOut.close();

    DDbarOut.close();
    cos2Phi_DDbarOut.close();

    DBDBbarOut.close();
    cos2Phi_DBDBbarOut.close();

    JPsiOut.close();
    cos2Phi_JPsiOut.close();
    OnlyDirectJPsiOut.close();
    cos2Phi_OnlyDirectJPsiOut.close();

    ExcitedJPsiOut.close();
    cos2Phi_ExcitedJPsiOut.close();
    OnlyDirectExcitedJPsiOut.close();
    cos2Phi_OnlyDirectExcitedJPsiOut.close();

    UpsilonOut.close();

    ExcitedUpsilonOut.close();

    bcbarOut.close();

    ExcitedbcbarOut.close();
    
    ccbarBeforeEvolutiondNdyOut.close();
    ccbardNdyOut.close();
    DDbardNdyOut.close();
    JPsidNdyOut.close();
    ExcitedJPsidNdyOut.close();

  }
    
  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
