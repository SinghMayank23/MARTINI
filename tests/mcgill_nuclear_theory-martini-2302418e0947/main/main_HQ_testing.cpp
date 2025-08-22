#include "MARTINI.h"
#include <sstream>
#define numJets 2000
#define nybins 100
#define DNDY_MAX 5.
#define nptbins 50
//#define Y_MAX 0.35
#define Y_MAX 5.
#define PT_MAX 5.
#define R_MAX 10.
#define outputting false

using namespace std;

int main(int argc, char* argv[]) 
{
  //The first bit of testing for heavy quarks!
  MARTINI martini;
  Event fullEvent;

  martini.readFile("setup_files/setup_main_HQ_testing.dat");
  martini.init(argc, argv);
  martini.settings.listAll();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  string directory_name = "";
  //string directory_name = "HQ_results_and_analysis/dN2piptdptdy_AuAu_0t20/";
  string HQ_name = "RHIC_AuAu_0t20_UPotential_";

  int numEvents;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps

  Vec4 itsP;                          // Its P

  //Define the histograms:
  double dN2piptdptdy_ccbar_beforeEvolution[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution[nptbins];
  double dN2piptdptdy_ccbar[nptbins];
  double dN2piptdptdy_bbbar[nptbins];
  double dN2piptdptdy_DDbar[nptbins];
  double dN2piptdptdy_DBDBbar[nptbins];
  double dN2piptdptdy_JPsi[nptbins];
  double dN2piptdptdy_ExcitedJPsi[nptbins];
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
    dN2piptdptdy_DDbar[ip] = 0.;
    dN2piptdptdy_DBDBbar[ip] = 0.;
    dN2piptdptdy_JPsi[ip] = 0.;
    dN2piptdptdy_ExcitedJPsi[ip] = 0.;
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

  double dN4pir2dr_ccbar_beforeEvolution[nptbins];
  double dN4pir2dr_ccbar[nptbins];
  for(int ir=0; ir<nptbins; ir++){
    dN4pir2dr_ccbar_beforeEvolution[ir] = 0.;
    dN4pir2dr_ccbar[ir] = 0.;
  }

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
      //cout << "event_counter = " << event_counter << endl;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      martini.generateEvent_HQ(plist[0]); 
      //cout << "OK after generateEvent_HQ..." << endl;

      int eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );

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

	double dx = plist[0]->at(ie).x() - plist[0]->at(plist[0]->at(ie).antiI() ).x();
	double dy = plist[0]->at(ie).y() - plist[0]->at(plist[0]->at(ie).antiI() ).y();
	double dz = plist[0]->at(ie).z() - plist[0]->at(plist[0]->at(ie).antiI() ).z();
	double r = sqrt(dx*dx + dy*dy + dz*dz);
	int ir = (int)((r/R_MAX)*nptbins);
	if(ir < nptbins){
	  dN4pir2dr_ccbar_beforeEvolution[ir] += (R_MAX/(double)nptbins)/(4.*M_PI*R_MAX*R_MAX*((double)ir+0.5)*((double)ir+0.5)/(double)(nptbins*nptbins) ) ;
	}
      }
	  
      int counter = 0;
      if (martini.returnEvolution() == 1)
	{
	  for(int i=0; i<mt; i++)
	    {
	      counter = martini.evolve_HQ(plist, counter, i);
	      counter+=1;
	    }
	}
      //cout << "OK after evolve_HQ..." << endl;

      eventSize = plist[0]->size();

      //cout << "eventSize = " << eventSize << endl;

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );

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
	  }
	  if(abs(martini.pythia.event[ie].id())==5){
	    dN2piptdptdy_bbbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	}

	double dx = plist[0]->at(ie).x() - plist[0]->at(plist[0]->at(ie).antiI() ).x();
	double dy = plist[0]->at(ie).y() - plist[0]->at(plist[0]->at(ie).antiI() ).y();
	double dz = plist[0]->at(ie).z() - plist[0]->at(plist[0]->at(ie).antiI() ).z();
	double r = sqrt(dx*dx + dy*dy + dz*dz);
	int ir = (int)((r/R_MAX)*nptbins);
	if(ir < nptbins){
	  dN4pir2dr_ccbar[ir] += (R_MAX/(double)nptbins)/(4.*M_PI*R_MAX*R_MAX*((double)ir+0.5)*((double)ir+0.5)/(double)(nptbins*nptbins) ) ;
	}
      }
	  
      martini.fragmentation_HQ(plist);
      //cout << "OK after fragmentation_HQ..." << endl;

      eventSize = martini.pythia.event.size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = martini.pythia.event[ie].p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = itsP.mCalc();
	//double itsM = martini.pythia.event[ie].mass();
	//cout << "itsM: " << itsM << endl;
	double pz = itsP.pz();
	double y = 0.5*log( (sqrt(itsM*itsM+pt*pt+pz*pz) + pz)/(sqrt(itsM*itsM+pt*pt+pz*pz) - pz) );

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
	  }
	  if(abs(martini.pythia.event[ie].id())==511){
	    dN2piptdptdy_DBDBbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dN2piptdptdy_JPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NJPsi += 1.;
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dN2piptdptdy_ExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedJPsi += 1.;
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
    string bbbarname = directory_name+"dN2piptdpt_bbbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOut;
    bbbarOut.open(bbbarname.c_str());
    string DDbarname = directory_name+"dN2piptdpt_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOut;
    DDbarOut.open(DDbarname.c_str());
    string DBDBbarname = directory_name+"dN2piptdpt_DBDBbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DBDBbarOut;
    DBDBbarOut.open(DBDBbarname.c_str());
    string JPsiname = directory_name+"dN2piptdpt_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream JPsiOut;
    JPsiOut.open(JPsiname.c_str());
    string ExcitedJPsiname = directory_name+"dN2piptdpt_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedJPsiOut;
    ExcitedJPsiOut.open(ExcitedJPsiname.c_str());
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
      DDbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar[ip]/(2.*event_counter) << endl;
      DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;
      JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
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
      ccbarBeforeEvolutiondNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeEvolution[iy]/2. << endl;
      ccbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar[iy]/2. << endl;
      DDbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_DDbar[iy]/2. << endl;
      JPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_JPsi[iy]/2. << endl;
      ExcitedJPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ExcitedJPsi[iy]/2. << endl;
    }
    
    string ccbarBeforeEvolutiondNdrname = directory_name+"dN4pir2dr_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdrOut;
    ccbarBeforeEvolutiondNdrOut.open(ccbarBeforeEvolutiondNdrname.c_str());
    string ccbardNdrname = directory_name+"dN4pir2dr_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdrOut;
    ccbardNdrOut.open(ccbardNdrname.c_str());

    for(int ir=0; ir<nptbins; ir++){
      ccbarBeforeEvolutiondNdrOut << ((double)ir+0.5)*R_MAX/(double)nptbins << " " << dN4pir2dr_ccbar_beforeEvolution[ir]/(2.*event_counter) << endl;
      ccbardNdrOut << ((double)ir+0.5)*R_MAX/(double)nptbins << " " << dN4pir2dr_ccbar[ir]/(2.*event_counter) << endl;
    }

    ccbarBeforeEvolutionOut.close();
    bbbarBeforeEvolutionOut.close();
    ccbarOut.close();
    bbbarOut.close();
    DDbarOut.close();
    DBDBbarOut.close();
    JPsiOut.close();
    ExcitedJPsiOut.close();
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
