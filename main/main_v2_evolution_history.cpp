#include "MARTINI.h"
#include <sstream>
#define numJets 250
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

  //Define the histograms:
  double NJPsi_t[mt/10];
  double NJPsi_t_rec[mt/10];
  cout << "mt/10 = " << mt/10 << endl;

  double hydro_T[1000];
  double prehydro_T[100];
  int hydro_counter[1000];
  int prehydro_counter[100];
  
  double **dNdmomrap_ccbar_evolution;
  dNdmomrap_ccbar_evolution = new double*[ntau];
  for (int it = 0; it < ntau; it++){
    dNdmomrap_ccbar_evolution[it] = new double[nrapbins];
    for(int iy=0; iy<nrapbins; iy++){
      dNdmomrap_ccbar_evolution[it][iy] = 0.;
    }
  }
  
  //double **dNdz_ccbar_evolution;
  //dNdz_ccbar_evolution = new double*[mt];
  //for (int it = 0; it < mt; it++){
  //  dNdz_ccbar_evolution[it] = new double[nzbins];
  //  for(int iy=0; iy<nzbins; iy++){
  //    dNdz_ccbar_evolution[it][iy] = 0.;
  //  }
  //}

  double **dNdsprap_ccbar_evolution;
  dNdsprap_ccbar_evolution = new double*[ntau];
  for (int it = 0; it < ntau; it++){
    dNdsprap_ccbar_evolution[it] = new double[nrapbins];
    for(int iy=0; iy<nrapbins; iy++){
      dNdsprap_ccbar_evolution[it][iy] = 0.;
    }
  }

  double ***dNdsprapmy_ccbar_evolution;
  dNdsprapmy_ccbar_evolution = new double**[ntau];
  for (int it = 0; it < ntau; it++){
    dNdsprapmy_ccbar_evolution[it] = new double*[nrapbins];
    for(int iy=0; iy<nrapbins; iy++){
      dNdsprapmy_ccbar_evolution[it][iy] = new double[2];
      dNdsprapmy_ccbar_evolution[it][iy][0] = 0.;
      dNdsprapmy_ccbar_evolution[it][iy][1] = 0.;
    }
  }

  double ***dNdsprapme_ccbar_evolution;
  dNdsprapme_ccbar_evolution = new double**[ntau];
  for (int it = 0; it < ntau; it++){
    dNdsprapme_ccbar_evolution[it] = new double*[nrapbins];
    for(int iy=0; iy<nrapbins; iy++){
      dNdsprapme_ccbar_evolution[it][iy] = new double[2];
      dNdsprapme_ccbar_evolution[it][iy][0] = 0.;
      dNdsprapme_ccbar_evolution[it][iy][1] = 0.;
    }
  }

  for (int ind = 0; ind < 100; ind++){
    prehydro_T[ind] = 0.;
    prehydro_counter[ind] = 0;
  }

  for (int ind = 0; ind < 1000; ind++){
    hydro_T[ind] = 0.;
    hydro_counter[ind] = 0;
  }

//  double delpt[largeptbins];
//  double midpt[largeptbins];
//  if (PT_MAX > 100.) {
//	  cout << "Error:: Cannot track largeptbins for diffusion stability such a large ptmax" << endl;
//	  exit(1);
//  }
//  delpt[0] = 100. - PT_MAX;
//  delpt[1] = 900.;
//  delpt[2] = 9000.;
//  delpt[3] = 90000.;
//  delpt[4] = 900000.;
//  delpt[5] = 9000000.;
//  delpt[6] = 90000000.;
//  delpt[7] = 900000000.;
//  midpt[0] = (100. + PT_MAX)/2.;
//  midpt[1] = 550.;
//  midpt[2] = 5500.;
//  midpt[3] = 55000.;
//  midpt[4] = 550000.;
//  midpt[5] = 5500000.;
//  midpt[6] = 55000000.;
//  midpt[7] = 550000000.;
//  if (largeptbins > 8) {
//	  cout << "Error:: Cannot track largeptbins for diffusion stability such a large largeptbins" << endl;
//	  exit(1);
//  }
//
  double dN2piptdptdy_ccbar_beforeEvolution[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution[nptbins];
//  double dN2piptdptdy_ccbar_beforeEvolution_largept[largeptbins];
  
  double dN2piptdptdy_ccbar[nptbins];
  double dN2piptdptdy_bbbar[nptbins];
//  double dN2piptdptdy_ccbar_largept[largeptbins];
  double cos2Phi_dN2piptdptdy_ccbar[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar[nptbins];

  double dN2piptdptdy_DDbar[nptbins];
  double dN2piptdptdy_DBDBbar[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar[nptbins];
  double cos2Phi_dN2piptdptdy_DBDBbar[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar[nptbins];
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

//  for(int ip=0; ip<largeptbins; ip++){
//    dN2piptdptdy_ccbar_beforeEvolution_largept[ip] = 0.;
//    dN2piptdptdy_ccbar_beforeHydro_largept[ip] = 0.;
//    dN2piptdptdy_ccbar_largept[ip] = 0.;
//  }

  for(int ip=0; ip<nptbins; ip++){
    dN2piptdptdy_ccbar_beforeEvolution[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution[ip] = 0.;

//    dN2piptdptdy_ccbar_beforeHydro[ip] = 0.;
//    dN2piptdptdy_bbbar_beforeHydro[ip] = 0.;

    dN2piptdptdy_ccbar[ip] = 0.;
    dN2piptdptdy_bbbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar[ip] = 0.;

    dN2piptdptdy_DDbar[ip] = 0.;
    dN2piptdptdy_DBDBbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar[ip] = 0.;
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
  double dNdmomrap_ccbar_beforeEvolution[nybins];
  double dNdmomrap_ccbar[nybins];
  double dNdstrap_ccbar[nybins];
  double dNdy_DDbar[nybins];
  double dNdy_JPsi[nybins];
  double dNdy_ExcitedJPsi[nybins];

  for(int iy=0; iy<nybins; iy++){
    dNdy_ccbar_beforeEvolution[iy] = 0.;
    dNdy_ccbar[iy] = 0.;
    dNdmomrap_ccbar_beforeEvolution[iy] = 0.;
    dNdmomrap_ccbar[iy] = 0.;
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

      cout << "event_counter  " << event_counter << endl;
      martini.generateEventHeavyQuarks_w_colllist_OR_IPGfile(plist[0]); 
//      martini.generateEventHeavyQuarks(plist[0]); 

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

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	//} else if (fabs(y)<Y_MAX) {
	//  int ipt = (int)(log(pt)/log(10.)) - 1;
	//  if(abs(martini.pythia.event[ie].id())==4 && ipt < largeptbins){
	//    dN2piptdptdy_ccbar_beforeEvolution_largept[ipt] += delpt[ipt]*(0.5/Y_MAX)/(2.*M_PI*midpt[ipt] );
	//  }
        }

      }

      double hydrotau0 = martini.returnHydroTau0();
      bool wrote_prehydro_list = false;

      int counter = 0;
      if (martini.returnEvolution() == 1)
	{
	  for(int i=0; i<mt; i++)
	    {
	      //counter = martini.evolveHeavyQuarks(plist, counter, i);
	      counter = martini.evolveAndHadronizeHeavyQuarks(plist, hydro_T, prehydro_T, hydro_counter, prehydro_counter, counter, i);
	      counter+=1;
	      double time = ((double)i+1.)*dtfm;
	      int eventsize = plist[0]->size();
	      for (int ie = 0; ie < eventsize; ie++){
	        double z = plist[0]->at(ie).z();
	        double tau = sqrt(time*time - z*z);
		for (int itau = 0; itau < ntau; itau++){
		  double ptau = ((double)itau + 1.)*0.1;
		  double prevtau = plist[0]->at(ie).prevtau();
		  if (tau > ptau && prevtau < ptau){
	             Vec4 itsP = plist[0]->at(ie).p();
	             double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	             double itsM = plist[0]->at(ie).mass();
	             double pz = itsP.pz();
	             double momrap = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );
	             double pseudorap = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	             int imomrap = (int)( (momrap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	             if(imomrap>=0 && imomrap<nrapbins){
	               if(abs(plist[0]->at(ie).id())==4 && plist[0]->at(ie).frozen()==0){
	                 dNdmomrap_ccbar_evolution[itau][imomrap] += (double)nrapbins/(2.*DNDY_MAX);
		       }
	//	       break;
		     }
		     double sprap = atanh(z/time);
		     int isprap = (int)( (sprap + DNDY_MAX)*(double)nrapbins/(2.*DNDY_MAX));
		     if (isprap >= 0 && isprap < nrapbins){
		       if (abs(plist[0]->at(ie).id())==4 && plist[0]->at(ie).frozen()==0){
		         dNdsprap_ccbar_evolution[itau][isprap] += (double)nrapbins/(2.*DNDY_MAX);
		         dNdsprapmy_ccbar_evolution[itau][isprap][0] += sprap - momrap;
		         dNdsprapmy_ccbar_evolution[itau][isprap][1] += 1.;
		         dNdsprapme_ccbar_evolution[itau][isprap][0] += sprap - pseudorap;
		         dNdsprapme_ccbar_evolution[itau][isprap][1] += 1.;
		       }
	//	       break;
		     }
	          }
	        }

		plist[0]->at(ie).prevtau(tau);
		//int iz = (int)( (z+zmax)*(double)nzbins/(2.*zmax) );
		//if(iz >= 0 && iz < nzbins){
	        //  if(abs(plist[0]->at(ie).id())==4){
	        //    dNdz_ccbar_evolution[i][iz] += (double)nzbins/(2.*zmax);
	        //  }
	        //}
	      }
	//      if (!wrote_prehydro_list && time >= hydrotau0)
	//	{
	//	  wrote_prehydro_list = true;
        //          int eventSize = plist[0]->size();

        //          for(int ie=0; ie<eventSize; ie++){
        //            Vec4 itsP = plist[0]->at(ie).p();
        //            double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
        //            double itsM = plist[0]->at(ie).mass();
        //            double pz = itsP.pz();
        //            //Pseudorapidity is more appropriate here:
        //            double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );

        //            //The dN/dy histograms:
        //            int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
        //            if(iy>=0 && iy<nybins){
        //              if(abs(martini.pythia.event[ie].id())==4){
        //                dNdy_ccbar_beforeHydro[iy] += (double)nybins/(2.*DNDY_MAX);
        //              }
        //            }

        //            if(fabs(y)<Y_MAX && pt<PT_MAX){
        //              int ipt = (int)((pt/PT_MAX)*nptbins );

        //              if(abs(martini.pythia.event[ie].id())==4){
        //                dN2piptdptdy_ccbar_beforeHydro[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
        //              }
        //              if(abs(martini.pythia.event[ie].id())==5){
        //                dN2piptdptdy_bbbar_beforeHydro[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
        //              }
	//            } else if (fabs(y)<Y_MAX) {
	//              int ipt = (int)(log(pt)/log(10.)) - 1;
	//              if(abs(martini.pythia.event[ie].id())==4 && ipt < largeptbins){
	//                dN2piptdptdy_ccbar_beforeHydro_largept[ipt] += delpt[ipt]*(0.5/Y_MAX)/(2.*M_PI*midpt[ipt] );
	//              }
        //            }

        //          }

	//	}
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
	double momrap = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );
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
	int imomrap = (int)( (momrap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(imomrap>=0 && imomrap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdmomrap_ccbar[imomrap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	//} else if (fabs(y)<Y_MAX) {
	//  int ipt = (int)(log(pt)/log(10.)) - 1;
	//  if(abs(martini.pythia.event[ie].id())==4 && ipt < largeptbins){
	//    dN2piptdptdy_ccbar_largept[ipt] += delpt[ipt]*(0.5/Y_MAX)/(2.*M_PI*midpt[ipt] );
	//  }
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

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
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

    string ccbarBeforeEvolutionname = directory_name+"dN2piptdpt_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOut;
    ccbarBeforeEvolutionOut.open(ccbarBeforeEvolutionname.c_str());

    string bbbarBeforeEvolutionname = directory_name+"dN2piptdpt_bbbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
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
    
    string PrehydroTname = directory_name+"PrehydroT_"+HQ_name+"_"+s.str()+".dat";
    ofstream PrehydroTOut;
    PrehydroTOut.open(PrehydroTname.c_str());
    
    string hydroTname = directory_name+"hydroT_"+HQ_name+"_"+s.str()+".dat";
    ofstream hydroTOut;
    hydroTOut.open(hydroTname.c_str());
    
    //for(int ip=0; ip<largeptbins; ip++){
    //  ccbarBeforeEvolutionLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_beforeEvolution_largept[ip]/(2.*event_counter) << endl;
    //  ccbarBeforeHydroLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_beforeHydro_largept[ip]/(2.*event_counter) << endl;
    //  ccbarLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_largept[ip]/(2.*event_counter) << endl;
    //}

    for(int ip=0; ip<nptbins; ip++){
      ccbarBeforeEvolutionOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution[ip]/(2.*event_counter) << endl;

//      ccbarBeforeHydroOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeHydro[ip]/(2.*event_counter) << endl;
//      bbbarBeforeHydroOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeHydro[ip]/(2.*event_counter) << endl;

      ccbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar[ip]/(2.*event_counter) << endl;
      bbbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar[ip]/(2.*event_counter) << endl;

      DDbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar[ip]/(2.*event_counter) << endl;
      DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar[ip]/(2.*event_counter) << endl;
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
  
    for (int ind = 0; ind < 100; ind++){
      //if (prehydro_counter[ind] == 0) break;
      PrehydroTOut << ind << " " << prehydro_T[ind] << " " << prehydro_counter[ind] << endl;
    }

    for (int ind = 0; ind < 1000; ind++){
      if (hydro_counter[ind] == 0) break;
      hydroTOut << ind << " " << hydro_T[ind] << " " << hydro_counter[ind] << endl;
    }

    string ccbarBeforeEvolutiondNdyname = directory_name+"dNdy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdyOut;
    ccbarBeforeEvolutiondNdyOut.open(ccbarBeforeEvolutiondNdyname.c_str());

    string ccbarBeforeEvolutiondNdmomrapname = directory_name+"dNdmomrap_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdmomrapOut;
    ccbarBeforeEvolutiondNdmomrapOut.open(ccbarBeforeEvolutiondNdmomrapname.c_str());

//    string ccbarBeforeHydrodNdyname = directory_name+"dNdy_ccbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydrodNdyOut;
//    ccbarBeforeHydrodNdyOut.open(ccbarBeforeHydrodNdyname.c_str());

    string ccbardNdsprapevolname = directory_name+"dNdsprap_ccbar_evolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdsprapevolOut;
    ccbardNdsprapevolOut.open(ccbardNdsprapevolname.c_str());

    string ccbardNdmomrapevolname = directory_name+"dNdmomrap_ccbar_evolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdmomrapevolOut;
    ccbardNdmomrapevolOut.open(ccbardNdmomrapevolname.c_str());

    string ccbardNdyname = directory_name+"dNdy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdyOut;
    ccbardNdyOut.open(ccbardNdyname.c_str());

    string ccbardNdmomrapname = directory_name+"dNdmomrap_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdmomrapOut;
    ccbardNdmomrapOut.open(ccbardNdmomrapname.c_str());

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
      ccbarBeforeEvolutiondNdmomrapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdmomrap_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
    //  ccbarBeforeHydrodNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeHydro[iy]/(2.*(double)event_counter) << endl;
      ccbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar[iy]/(2.*(double)event_counter) << endl;
      ccbardNdmomrapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdmomrap_ccbar[iy]/(2.*(double)event_counter) << endl;
      DDbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_DDbar[iy]/(2.*(double)event_counter) << endl;
      JPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_JPsi[iy]/(2.*(double)event_counter) << endl;
      ExcitedJPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ExcitedJPsi[iy]/(2.*(double)event_counter) << endl;
    }

    for(int itau = 0; itau < ntau; itau++){
      double tau = 0.1*(double)itau + 0.1;
      //double rapmax1 = acosh(maxTime/tau);
      //double rapmax2 = asinh(zmax/tau);
      //double rapmax = min(rapmax1, rapmax2);
      for(int irap = 0; irap < nrapbins; irap++){
        double rap = ((double)irap+0.5)*2.*DNDY_MAX/(double)nrapbins - DNDY_MAX;
	//double timet = tau*cosh(rap);
	//double spacz = tau*sinh(rap);
        //double izd = (spacz + zmax)*(double)nzbins/(2.*zmax);
	//double itd = timet/dtfm;
	//if (itd < 0 || izd < 0) continue;

	//int it = static_cast<int>(itd);
	//int iz = static_cast<int>(izd);

	//if (it >= (mt - 1) || iz >= (nzbins - 1)) continue;

	//double rest = itd - (double)it;
	//double resz = izd - (double)iz;

	//double value = dNdz_ccbar_evolution[it][iz]*(1.-rest)*(1.-resz) + dNdz_ccbar_evolution[it+1][iz]*rest*(1.-resz)
	//             + dNdz_ccbar_evolution[it][iz+1]*(1.-rest)*resz + dNdz_ccbar_evolution[it+1][iz+1]*rest*resz;
	//value *= tau;
        //ccbardNdsprapevolOut << tau << " " << rap << " " << value/(2.*(double)event_counter) << endl; 
        ccbardNdsprapevolOut << tau << " " << rap << " " << dNdsprap_ccbar_evolution[itau][irap]/(2.*(double)event_counter) << " " << dNdsprapmy_ccbar_evolution[itau][irap][0]/dNdsprapmy_ccbar_evolution[itau][irap][1] << " " << dNdsprapme_ccbar_evolution[itau][irap][0]/dNdsprapme_ccbar_evolution[itau][irap][1] << endl; 
	
      }
      ccbardNdsprapevolOut << endl;
    }

    for(int itau = 0; itau < ntau; itau++){
      double tau = 0.1*(double)itau + 0.1;
      //int itt = itau*10;
      //if (itt >= mt) break;
      //double time = (double)itt*dtfm;
      for(int iy = 0; iy < nrapbins; iy++){
        ccbardNdmomrapevolOut << tau << " " << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdmomrap_ccbar_evolution[itau][iy]/(2.*(double)event_counter) << endl;
      }
      ccbardNdmomrapevolOut << endl;
    }

    PrehydroTOut.close();    
    hydroTOut.close();   

    ccbardNdmomrapevolOut.close();
    ccbardNdsprapevolOut.close();

    ccbarBeforeEvolutionOut.close();

    bbbarBeforeEvolutionOut.close();

//    ccbarBeforeHydroOut.close();

//    bbbarBeforeHydroOut.close();

//    ccbarBeforeEvolutionLargeptOut.close();
//    ccbarBeforeHydroLargeptOut.close();
//    ccbarLargeptOut.close();

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
    ccbarBeforeEvolutiondNdmomrapOut.close();
    //ccbarBeforeHydrodNdyOut.close();
    ccbardNdyOut.close();
    ccbardNdmomrapOut.close();
    DDbardNdyOut.close();
    JPsidNdyOut.close();
    ExcitedJPsidNdyOut.close();

  }

  for (int it = 0; it < mt; it++){
    delete dNdmomrap_ccbar_evolution[it];
    //delete dNdz_ccbar_evolution[it];
  }
  for (int it = 0; it < ntau; it++){
    delete dNdsprap_ccbar_evolution[it];
    for (int iy = 0; iy < nybins; iy++){
      delete dNdsprapmy_ccbar_evolution[it][iy];
      delete dNdsprapme_ccbar_evolution[it][iy];
    }
    delete dNdsprapmy_ccbar_evolution[it];
    delete dNdsprapme_ccbar_evolution[it];
  }

  delete dNdmomrap_ccbar_evolution;
  //delete dNdz_ccbar_evolution;
  delete dNdsprap_ccbar_evolution;
  delete dNdsprapmy_ccbar_evolution;
  delete dNdsprapme_ccbar_evolution;

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
