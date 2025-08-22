#include "MARTINI.h"
#include <sstream>
#define numJets 5000000
#define nybins 200
#define DNDY2_MAX 10.
#define DNDY_MAX 10.
#define nptbins 50
#define Y_MAX 1.0
#define Y_MAX_ALICE_v2 0.8
#define Y_MAX_ALICE_RAA 0.5
#define Y_MAX_CMS 1.
#define PT_MAX 15.
#define R_MAX 10.
#define outputting true

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

  Vec4 itsP;                          // Its P

  //Define the histograms:
  double NJPsi_t[mt/10];
  double NJPsi_t_rec[mt/10];
  cout << "mt/10 = " << mt/10 << endl;

  double dN2piptdptdy_ccbar_beforeEvolution_CMS[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_CMS[nptbins];
  double dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[nptbins];
  double dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[nptbins];
  
  double dN2piptdptdy_ccbar_CMS[nptbins];
  double dN2piptdptdy_bbbar_CMS[nptbins];
  double dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double dN2piptdptdy_bbbar_ALICE_RAA[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[nptbins];

  double dN2piptdptdy_DDbar_CMS[nptbins];
  double dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double dN2piptdptdy_DBDBbar[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double cos2Phi_dN2piptdptdy_DBDBbar[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_DBDBbar[nptbins];

  double dN2piptdptdy_JPsi[nptbins];
  double dN2piptdptdy_ExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_JPsi[nptbins];
  double cos2Phi_dN2piptdptdy_ExcitedJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_JPsi[nptbins];
  double sin2Phi_dN2piptdptdy_ExcitedJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];

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
    dN2piptdptdy_ccbar_beforeEvolution_CMS[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_CMS[ip] = 0.;

    dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ip] = 0.;

    dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ip] = 0.;

    dN2piptdptdy_ccbar_CMS[ip] = 0.;
    dN2piptdptdy_bbbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_CMS[ip] = 0.;

    dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;

    dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;

    dN2piptdptdy_DDbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_CMS[ip] = 0.;

    dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;

    dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;

    dN2piptdptdy_DBDBbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;
    sin2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;

    dN2piptdptdy_JPsi[ip] = 0.;
    dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_JPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_JPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;

    dN2piptdptdy_Upsilon[ip] = 0.;
    dN2piptdptdy_ExcitedUpsilon[ip] = 0.;

    dN2piptdptdy_bcbar[ip] = 0.;
    dN2piptdptdy_Excitedbcbar[ip] = 0.;
  }

  double dNdy_ccbar_beforeEvolution[nybins];
  double dNdy_ccbar[nybins];
  double dNdpseudorap_ccbar_beforeEvolution[nybins];
  double dNdpseudorap_ccbar[nybins];
  double dNdstrap_ccbar[nybins];
  double dNdy_DDbar[nybins];
  double dNdy_JPsi[nybins];
  double dNdy_ExcitedJPsi[nybins];

  for(int iy=0; iy<nybins; iy++){
    dNdy_ccbar_beforeEvolution[iy] = 0.;
    dNdy_ccbar[iy] = 0.;
    dNdpseudorap_ccbar_beforeEvolution[iy] = 0.;
    dNdpseudorap_ccbar[iy] = 0.;
    dNdstrap_ccbar[iy] = 0.;
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

  cout << numJets << endl;

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

     if(event_counter%10000 == 0) cout << "event_counter  " << event_counter << endl;
//      martini.generateEventHeavyQuarks_w_colllist_OR_IPGfile(plist[0]); 
      martini.generateEventHeavyQuarks(plist[0]); 

      int eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	double pseudorap = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y         = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdy_ccbar_beforeEvolution[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}
	int ipseudorap = (int)( (pseudorap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(ipseudorap>=0 && ipseudorap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdpseudorap_ccbar_beforeEvolution[ipseudorap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	}
	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

      }
	  
      //int counter = 0;
      //if (martini.returnEvolution() == 1)
      //  {
      //    for(int i=0; i<mt; i++)
      //      {
      //        //counter = martini.evolveHeavyQuarks(plist, counter, i);
      //        counter = martini.evolveAndHadronizeHeavyQuarks(plist, counter, i);
      //        counter+=1;
      //      }
      //  } else if (martini.returnEvolution() == 0)
      //  {
            int temp = martini.hadronizeHeavyQuarks(plist);
      //  }


      eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	double pseudorap = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );
	double cos2Phi = 0.;
	double sin2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	  sin2Phi = 2.*(itsP.px()*itsP.py())/(pt*pt);
	}

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdy_ccbar[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}
	int ipseudorap = (int)( (pseudorap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(ipseudorap>=0 && ipseudorap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdpseudorap_ccbar[ipseudorap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	}
	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
      }
	  
      eventSize = martini.pythia.event.size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = martini.pythia.event[ie].p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = itsP.mCalc();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	//double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y = 0.5*log( (sqrt(pt*pt+pz*pz + itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz + itsM*itsM) - pz) );
	double cos2Phi = 0.;
	double sin2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	  sin2Phi = 2.*(itsP.px()*itsP.py())/(pt*pt);
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

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==511){
	    dN2piptdptdy_DBDBbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DBDBbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DBDBbar[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dN2piptdptdy_JPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_JPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_JPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    }
	    NJPsi += 1.;
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dN2piptdptdy_ExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ExcitedJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedJPsi += 1.;
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
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

    string ccbarBeforeEvolutionnameCMS = directory_name+"dN2piptdpt_ccbar_beforeEvolution_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutCMS;
    ccbarBeforeEvolutionOutCMS.open(ccbarBeforeEvolutionnameCMS.c_str());

    string bbbarBeforeEvolutionnameCMS = directory_name+"dN2piptdpt_bbbar_beforeEvolution_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutCMS;
    bbbarBeforeEvolutionOutCMS.open(bbbarBeforeEvolutionnameCMS.c_str());

    string ccbarBeforeEvolutionnameALICEv2 = directory_name+"dN2piptdpt_ccbar_beforeEvolution_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutALICEv2;
    ccbarBeforeEvolutionOutALICEv2.open(ccbarBeforeEvolutionnameALICEv2.c_str());

    string bbbarBeforeEvolutionnameALICEv2 = directory_name+"dN2piptdpt_bbbar_beforeEvolution_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutALICEv2;
    bbbarBeforeEvolutionOutALICEv2.open(bbbarBeforeEvolutionnameALICEv2.c_str());

    string ccbarBeforeEvolutionnameALICERAA = directory_name+"dN2piptdpt_ccbar_beforeEvolution_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutALICERAA;
    ccbarBeforeEvolutionOutALICERAA.open(ccbarBeforeEvolutionnameALICERAA.c_str());

    string bbbarBeforeEvolutionnameALICERAA = directory_name+"dN2piptdpt_bbbar_beforeEvolution_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutALICERAA;
    bbbarBeforeEvolutionOutALICERAA.open(bbbarBeforeEvolutionnameALICERAA.c_str());

    string ccbarnameCMS = directory_name+"dN2piptdpt_ccbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutCMS;
    ccbarOutCMS.open(ccbarnameCMS.c_str());
    string cos2Phi_ccbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_ccbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutCMS;
    cos2Phi_ccbarOutCMS.open(cos2Phi_ccbarnameCMS.c_str());

    string bbbarnameCMS = directory_name+"dN2piptdpt_bbbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutCMS;
    bbbarOutCMS.open(bbbarnameCMS.c_str());
    string cos2Phi_bbbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_bbbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutCMS;
    cos2Phi_bbbarOutCMS.open(cos2Phi_bbbarnameCMS.c_str());

    string DDbarnameCMS = directory_name+"dN2piptdpt_DDbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutCMS;
    DDbarOutCMS.open(DDbarnameCMS.c_str());
    string cos2Phi_DDbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_DDbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutCMS;
    cos2Phi_DDbarOutCMS.open(cos2Phi_DDbarnameCMS.c_str());

    string ccbarnameALICEv2 = directory_name+"dN2piptdpt_ccbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutALICEv2;
    ccbarOutALICEv2.open(ccbarnameALICEv2.c_str());
    string cos2Phi_ccbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_ccbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutALICEv2;
    cos2Phi_ccbarOutALICEv2.open(cos2Phi_ccbarnameALICEv2.c_str());

    string bbbarnameALICEv2 = directory_name+"dN2piptdpt_bbbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutALICEv2;
    bbbarOutALICEv2.open(bbbarnameALICEv2.c_str());
    string cos2Phi_bbbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_bbbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutALICEv2;
    cos2Phi_bbbarOutALICEv2.open(cos2Phi_bbbarnameALICEv2.c_str());

    string DDbarnameALICEv2 = directory_name+"dN2piptdpt_DDbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutALICEv2;
    DDbarOutALICEv2.open(DDbarnameALICEv2.c_str());
    string cos2Phi_DDbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_DDbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutALICEv2;
    cos2Phi_DDbarOutALICEv2.open(cos2Phi_DDbarnameALICEv2.c_str());

    string ccbarnameALICERAA = directory_name+"dN2piptdpt_ccbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutALICERAA;
    ccbarOutALICERAA.open(ccbarnameALICERAA.c_str());
    string cos2Phi_ccbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_ccbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutALICERAA;
    cos2Phi_ccbarOutALICERAA.open(cos2Phi_ccbarnameALICERAA.c_str());

    string bbbarnameALICERAA = directory_name+"dN2piptdpt_bbbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutALICERAA;
    bbbarOutALICERAA.open(bbbarnameALICERAA.c_str());
    string cos2Phi_bbbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_bbbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutALICERAA;
    cos2Phi_bbbarOutALICERAA.open(cos2Phi_bbbarnameALICERAA.c_str());

    string DDbarnameALICERAA = directory_name+"dN2piptdpt_DDbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutALICERAA;
    DDbarOutALICERAA.open(DDbarnameALICERAA.c_str());
    string cos2Phi_DDbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_DDbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutALICERAA;
    cos2Phi_DDbarOutALICERAA.open(cos2Phi_DDbarnameALICERAA.c_str());

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
      ccbarBeforeEvolutionOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_CMS[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_CMS[ip]/(2.*event_counter) << endl;

      ccbarBeforeEvolutionOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ip]/(2.*event_counter) << endl;

      ccbarBeforeEvolutionOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ip]/(2.*event_counter) << endl;

      ccbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << endl;
      bbbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << endl;

      ccbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      bbbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << endl;

      ccbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      bbbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << endl;

      DDbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << endl;

      DDbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << endl;

      DDbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << endl;

      DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;
      cos2Phi_DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;

      JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_JPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      cos2Phi_ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;

      UpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Upsilon[ip]/(double)event_counter << endl;
      ExcitedUpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedUpsilon[ip]/(double)event_counter << endl;

      bcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bcbar[ip]/(double)event_counter << endl;
      ExcitedbcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Excitedbcbar[ip]/(double)event_counter << endl;
    }
  
    string ccbarBeforeEvolutiondNdyname = directory_name+"dNdy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdyOut;
    ccbarBeforeEvolutiondNdyOut.open(ccbarBeforeEvolutiondNdyname.c_str());

    string ccbarBeforeEvolutiondNdpseudorapname = directory_name+"dNdpseudorap_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdpseudorapOut;
    ccbarBeforeEvolutiondNdpseudorapOut.open(ccbarBeforeEvolutiondNdpseudorapname.c_str());

    string ccbardNdyname = directory_name+"dNdy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdyOut;
    ccbardNdyOut.open(ccbardNdyname.c_str());

    string ccbardNdpseudorapname = directory_name+"dNdpseudorap_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdpseudorapOut;
    ccbardNdpseudorapOut.open(ccbardNdpseudorapname.c_str());

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
      ccbarBeforeEvolutiondNdpseudorapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdpseudorap_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
      ccbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar[iy]/(2.*(double)event_counter) << endl;
      ccbardNdpseudorapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdpseudorap_ccbar[iy]/(2.*(double)event_counter) << endl;
      DDbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_DDbar[iy]/(2.*(double)event_counter) << endl;
      JPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_JPsi[iy]/(2.*(double)event_counter) << endl;
      ExcitedJPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ExcitedJPsi[iy]/(2.*(double)event_counter) << endl;
    }
    

    ccbarBeforeEvolutionOutCMS.close();
    bbbarBeforeEvolutionOutCMS.close();

    ccbarBeforeEvolutionOutALICEv2.close();
    bbbarBeforeEvolutionOutALICEv2.close();

    ccbarBeforeEvolutionOutALICERAA.close();
    bbbarBeforeEvolutionOutALICERAA.close();

    ccbarOutCMS.close();
    bbbarOutCMS.close();

    ccbarOutALICEv2.close();
    bbbarOutALICEv2.close();

    ccbarOutALICERAA.close();
    bbbarOutALICERAA.close();

    cos2Phi_ccbarOutCMS.close();
    cos2Phi_bbbarOutCMS.close();

    cos2Phi_ccbarOutALICEv2.close();
    cos2Phi_bbbarOutALICEv2.close();

    cos2Phi_ccbarOutALICERAA.close();
    cos2Phi_bbbarOutALICERAA.close();

    DDbarOutCMS.close();
    cos2Phi_DDbarOutCMS.close();

    DDbarOutALICEv2.close();
    cos2Phi_DDbarOutALICEv2.close();

    DDbarOutALICERAA.close();
    cos2Phi_DDbarOutALICERAA.close();

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
    ccbarBeforeEvolutiondNdpseudorapOut.close();
    ccbardNdyOut.close();
    ccbardNdpseudorapOut.close();
    DDbardNdyOut.close();
    JPsidNdyOut.close();
    ExcitedJPsidNdyOut.close();

  }

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
