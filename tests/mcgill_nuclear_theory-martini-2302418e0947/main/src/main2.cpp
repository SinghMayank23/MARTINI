#include "MARTINI.h"

int main(int argc, char* argv[]) 
{
  // the following three commands are essential
  MARTINI martini;
  martini.readFile("setup.dat");
  martini.readString("General:Evolution = off");
  martini.init(argc, argv);
  // this lists MARTINI's current settings
  martini.settings.listAll();
  


  Parton jp1;                       // used for fixed energy runs
  Parton jp2;

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plist[0] = new vector<Parton>;    // plist[0] is the main list that holds the high momentum partons that are evolved

  double flucMeasure = 0.;
  Vec4 vecp, vecp2;                 // pythia four-vector object to store four-momenta
  int totalNNs;
  double p;                         // momentum |p| of the parton 
  int counter = 0;                  // counter will provide unique color ID for produced partons
  int numEvents = 0;                // will store number of actually accepted events (some are rejected by fragmentation)
  int mt;                           // maximal time steps
  int posy;                         // bin number 
  int p_imax;                       // maximum number of particles in list
  int id;                           // parton ID
  int pisum = 0;                    
  int piPlusSum = 0, piMinusSum=0;;                    
  int KPlusSum = 0, KMinusSum=0;                    
  const int bins = 80;              // number of bins
  const int phiBins = 6;            // number of bins
  double scale = 20.;               // maximum p_t in GeV in the binning (bin size = scale/bins)
  double etascale = 4.;             // maximum eta
  double tscale = PI;               // maximum theta
  double rscale = 10.;              // maximum r
  double xscale = 20.;              // maximum x
  int binning[bins];                // array with all the bins
  int binning_piplus[bins];         // array with all the bins
  int binning_Kplus[bins];          // array with all the bins
  int binning_piminus[bins];         // array with all the bins
  int binning_Kminus[bins];          // array with all the bins
  int binning_tmp[bins];            // array with all the bins
  int binning_sq[bins];             // array with all the bins for squares for error
  int binning_hm[bins];             // array with all the bins
  int binning_g[bins];              // array with all the bins for gluons
  int binning_q[bins];              // array with all the bins for quarks
  int binning_gamma[bins];          // array with all the bins for photons
  int binning_gamma_sq[bins];       // array with all the bins for photons for squares for error
  int binning_gamma_tmp[bins];      // array with all the bins for photons
  int binning_r[bins];              // array with all the bins for the radius in the transverse plane
  double binning_x[bins][bins];        // array with all the bins for the xy positions in the transverse plane
  int binning_x1_x2[bins][bins];    // array with bins in dalitz x1 and x2
  double binning_eta_phi[bins][bins];  // array with bins in eta and phi for lego plot
  int binning_phi_pt[phiBins][bins];// array with bins in phi and pt
  int binning_phi_tmp[phiBins][bins];// array with bins in phi and pt
  int binning_phi_sq[phiBins][bins];// array with bins in phi and pt for squares for error
  int binning_phi_pt_g[phiBins][bins];// array with bins in phi and pt
  int binning_phi_tmp_g[phiBins][bins];// array with bins in phi and pt
  int binning_phi_sq_g[phiBins][bins];// array with bins in phi and pt for squares for error
  int binning_phi_pt_q[phiBins][bins];// array with bins in phi and pt
  int binning_phi_tmp_q[phiBins][bins];// array with bins in phi and pt
  int binning_phi_sq_q[phiBins][bins];// array with bins in phi and pt for squares for error
  double binning_mom[bins][bins];   
  int binning_theta[bins];          // array with all the bins for the angle theta
  int totalSum = 0;
  int otherPartons = 0;                
  int hmsum = 0;
  double ymax = 1; //0.5 for one unit of rapidity around y=0.
  double pt2=0.;
  double pt2r=0.;
  double Etot=0;
  double NColTot=0.;
  double qtTotAll=0.;
  double Erun;
  double Ejet=10.;
  double r, theta;
  double EJetTotal;
  int posr, postheta, posxi, posyi, posphi;
  int countAll=0, countLarge=0;
  int totalJets = 0;
  int twoJetEvents = 0;
  int fakeJetEvents = 0;
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  int runs=martini.returnRuns();
  
  // make them both quarks with x GeV for testing:
  
  if ( martini.returnFixedEnergy() == 1 )
    {
      jp1.id(1); jp2.id(1);
      jp1.p(Ejet,0.,0.); jp2.p(Ejet,0.,0.);
      jp1.col(101); jp1.acol(102); jp2.col(103); jp2.acol(104);
      jp1.x(4.); jp1.y(0.); jp1.z(0.);
      jp2.x(4.); jp2.y(0.); jp2.z(0.);
      //jp1.splits(0); jp2.splits(0);
    }
  
  // init binning
  for(int iy=0; iy<bins; iy++)
    {
      binning_r[iy]=0;
      binning_theta[iy]=0;
      binning[iy]=0;
      binning_piplus[iy]=0;
      binning_Kplus[iy]=0;
      binning_piminus[iy]=0;
      binning_Kminus[iy]=0;
      binning_tmp[iy]=0;
      binning_sq[iy]=0;
      binning_hm[iy]=0;
      binning_q[iy]=0;
      binning_g[iy]=0;
      binning_gamma[iy]=0;
      binning_gamma_sq[iy]=0;
      binning_gamma_tmp[iy]=0;
    }

  for(int iy=0; iy<bins; iy++)
    for(int ix=0; ix<bins; ix++)
      {
	binning_x[iy][ix]=0.;
	binning_mom[iy][ix]=0.;
      }

  for(int iy=0; iy<bins; iy++)
    for(int ix=0; ix<phiBins; ix++)
      {
	binning_phi_sq[ix][iy]=0.;
	binning_phi_pt[ix][iy]=0.;
	binning_phi_tmp[ix][iy]=0.;
	binning_phi_sq_q[ix][iy]=0.;
	binning_phi_pt_q[ix][iy]=0.;
	binning_phi_tmp_q[ix][iy]=0.;
	binning_phi_sq_g[ix][iy]=0.;
	binning_phi_pt_g[ix][iy]=0.;
	binning_phi_tmp_g[ix][iy]=0.;
      }

  for(int iy=0; iy<bins; iy++)
    for(int ix=0; ix<bins; ix++)
      {
	binning_x1_x2[ix][iy]=0;
	binning_eta_phi[ix][iy]=0.;
      }

  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  
  ofstream foutt("test.dat",ios::app); 

  // ****************************************************************************************** //

  for (int j=0; j<runs; j++)      // loop over all events
    {
      martini.setQtTot(0.);
      martini.setNCol(0);
      pt2=0.;
      Erun=0.;

      plist[0]=martini.next();


      // BINNING ------------------------------------------------------------------------------ //

      // bin partons:
      p_imax = plist[0]->size();
      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	{
	  double pl,pt,theta, En, y, xini, yini;
	  id = plist[0]->at(p_i).id();
	  
	  pt = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
	  pl = plist[0]->at(p_i).p().pz(); //p_long
	  En = sqrt(pt*pt+pl*pl);
	  theta = atan(pt/pl);
	  //if (theta<0) cout << "theta=" << theta << endl;
	  if (theta<0) theta = PI+theta;
	  postheta = floor(theta*(bins/tscale));

	  r = sqrt(pow(plist[0]->at(p_i).xini(),2.)+pow(plist[0]->at(p_i).yini(),2.));
	  xini = plist[0]->at(p_i).xini();
	  yini = plist[0]->at(p_i).yini();
	  //cout << "x=" << xini << endl;
	  //cout << "y=" << yini << endl;
	  posr = floor(r*(bins/rscale));
	  posxi = floor((xini+10)*(bins/xscale));
	  posyi = floor((yini+10)*(bins/xscale));
	  //cout << "posxi=" << posxi << endl;
	  //cout << "posyi=" << posyi << endl;

	  //countSplits+=plist[0]->at(p_i).splits();
	  //avSplits+=plist[0]->at(p_i).splits()/plist[0]->size();
	  //if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[0]->at(p_i).splits();
	  // if (p>500.) countHardGluons++;
	  
	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
	  
	  if ( pt>0.) // bin initial positions of all partons that are within the intersting y range 
	                              //and above a certain p_t //// abs(y)<=ymax && 
	    {
	      //if( (id > 0 && id < 4) || id == 21 )
	      //{
		  if(posxi>0 && posxi<bins && posyi>0 && posyi<bins) 
		    {
		      binning_x[posxi][posyi]+=1./static_cast<double>(p_imax);
		      binning_mom[posxi][posyi]+=plist[0]->at(p_i).pini().px();
		      //cout << "adding one at " << posxi << ", " << posyi << endl;
		    }
		  if(posr>0 && posr<bins) binning_r[posr]+=1;
		  if(postheta>0 && postheta<bins) binning_theta[postheta]+=1;
		  //}
	    }
	}

      p_imax = plist[0]->size();
      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	{
	  double pl,pt,phi;
	  double En;
	  double y; // rapidity
	  id = plist[0]->at(p_i).id();
	  p = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
	  pt = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
	  Erun += sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
	  pt2 += (pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
	  posy = floor(pt*(bins/scale));
	  //cout << "pt=" << pt << ", posy=" << posy << endl;
	  pl = plist[0]->at(p_i).p().pz(); //p_long
	  En = sqrt(pt*pt+pl*pl);
	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
	  phi = asin(sqrt(pow(plist[0]->at(p_i).p().py(),2.))/pt);   // azimuthal angle
	  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi

	  //cout << "p=" << p << endl;
	  //countSplits+=plist[0]->at(p_i).splits();
	  //avSplits+=plist[0]->at(p_i).splits()/plist[0]->size();
	  //if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[0]->at(p_i).splits();
	  // if (p>500.) countHardGluons++;
	  
	  if (martini.returnFixedEnergy()) // don't use a rapidity cut for the fixed energy calculation
	    {
	      if( id > 0 && id < 4 ) if(posy>=0 && posy<bins) binning_q[posy]+=1;
	      if( id == 21 ) if(posy>=0 && posy<bins) binning_g[posy]+=1;
	    }
	  else if ( abs(y)<=ymax )
	    {
	      if( id > 0 && id < 4 )
		{
		  if(posy>=0 && posy<bins) 
		    {
		      binning_q[posy]+=1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_pt_q[posphi][posy] += 1;
		    }
		}
	      else if( id == 21 )
		{
		  if(posy>=0 && posy<bins) 
		    {
		      binning_g[posy]+=1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_pt_g[posphi][posy] += 1;
		    }
		}
	    }
	}
      
      for(int iy=0; iy<bins; iy++)
	{
	  for(int ix=0; ix<phiBins; ix++)
	    {
	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy]-binning_phi_tmp_q[ix][iy];
	      binning_phi_sq_q[ix][iy]+=binning_phi_tmp_q[ix][iy]*binning_phi_tmp_q[ix][iy];
	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy];
	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy]-binning_phi_tmp_g[ix][iy];
	      binning_phi_sq_g[ix][iy]+=binning_phi_tmp_g[ix][iy]*binning_phi_tmp_g[ix][iy];
	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy];
	    }
	}

      pt2r+=pt2/p_imax;
      Etot+=Erun/p_imax;
      NColTot+=static_cast<double>(martini.returnNCol())/static_cast<double>(p_imax);
      qtTotAll+=martini.returnQtTot()/p_imax;

      // fragmentation
      if (martini.returnFragmentationSwitch() == 1)
	{
	  if ( martini.returnFullEvent() == 0 ) 
	    {
	      numEvents+=martini.fragmentation( plist );
	      Vec4 boostVec;
	      Vec4 restP;
	      double eRest1, eRest2, eRest3;
	      double eRest1a, eRest2a, eRest3a;
	      double pxRest1, pyRest1, pzRest1, pRest1;
	      double pxRest2, pyRest2, pzRest2, pRest2;
	      double pxRest3, pyRest3, pzRest3, pRest3;
	      double px, py, pz;
	      double beta, cosPhi1, cosPhi2, cosPhi3, gamma;
	      double vx, vy, vz;
	      double EJetTotala;
	      double pl;
	      double En;
	      double y; // rapidity
	      double eta; //pseudo-rapidity
	      double phi; // angle with respect to the reaction plane
	      p_imax = martini.pythia.event.size();
	      
	      // jet finder
	      CellJet cellJet(3.); //3. is etaMax
	      double dalitzX[3];
	      double dalitzXa[3];
	      
	      cellJet.analyze2( martini.pythia.event, 50., 0.39, 2., 2.);
	      // jets will be automatically sorted with decreasing eT
	      if (cellJet.size()==3)
		{
		  // cout << "________________________" << endl;
		  
		  boostVec=cellJet.pMassive(0)+cellJet.pMassive(1)+cellJet.pMassive(2);
		  //cout << "boostVec=" << boostVec << endl;
		  boostVec/=boostVec.e();
		  vx = boostVec.px();
		  vy = boostVec.py();
		  vz = boostVec.pz();
		  
		  beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of boost velocity in units of c
		  //cout << "beta=" << sqrt(vx*vx+vy*vy+vz*vz) << endl;
		  
		  gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
		  
		  //cout << "beta=" << beta << endl;
		  //cout << "gamma=" << gamma << endl;
		  
		  //pRest = p * gamma * (1.-beta*cosPhi);                     
		  
		  eRest1 = gamma*cellJet.pMassive(0).e()-vx*gamma*cellJet.pMassive(0).px()
		    -vy*gamma*cellJet.pMassive(0).py()-vz*gamma*cellJet.pMassive(0).pz();
		  eRest2 = gamma*cellJet.pMassive(1).e()-vx*gamma*cellJet.pMassive(1).px()
		    -vy*gamma*cellJet.pMassive(1).py()-vz*gamma*cellJet.pMassive(1).pz();
		  eRest3 = gamma*cellJet.pMassive(2).e()-vx*gamma*cellJet.pMassive(2).px()
		    -vy*gamma*cellJet.pMassive(2).py()-vz*gamma*cellJet.pMassive(2).pz();
		  
		  EJetTotal = eRest1+eRest2+eRest3;
		  if (eRest1>eRest2 && eRest1>eRest3)
		    {
		      dalitzX[0] = 2.*eRest1/EJetTotal;
		      if ( eRest2 > eRest3 )
			{
			  dalitzX[1] = 2.*eRest2/EJetTotal;
			  dalitzX[2] = 2.*eRest3/EJetTotal;
			}
		      else 
			{
			  dalitzX[1] = 2.*eRest3/EJetTotal;
			  dalitzX[2] = 2.*eRest2/EJetTotal;
			}
		    }
		  else if (eRest2>eRest1 && eRest2>eRest3)
		    {
		      dalitzX[0] = 2.*eRest2/EJetTotal;
		      if ( eRest1 > eRest3 )
			{
			  dalitzX[1] = 2.*eRest1/EJetTotal;
			  dalitzX[2] = 2.*eRest3/EJetTotal;
			}
		      else 
			{
			  dalitzX[1] = 2.*eRest3/EJetTotal;
			  dalitzX[2] = 2.*eRest1/EJetTotal;
			}
		    }
		  else if (eRest3>eRest1 && eRest3>eRest2)
		    {
		      dalitzX[0] = 2.*eRest3/EJetTotal;
		      if ( eRest1 > eRest2 )
			{
			  dalitzX[1] = 2.*eRest1/EJetTotal;
			  dalitzX[2] = 2.*eRest2/EJetTotal;
			}
		      else 
			{
			  dalitzX[1] = 2.*eRest2/EJetTotal;
			  dalitzX[2] = 2.*eRest1/EJetTotal;
			}
		    }
		  
		  flucMeasure += (sqrt((eRest1-eRest2)*(eRest1-eRest2))+sqrt((eRest1-eRest3)*(eRest1-eRest3))
				  +sqrt((eRest2-eRest3)*(eRest2-eRest3)))/(eRest1+eRest2+eRest3);
		  		  
		  /*		  
		  pxRest1 = -vx*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(0).px() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(0).py()
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(0).pz();
		  pyRest1 = -vy*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(0).py() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(0).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(0).pz();
		  pzRest1 = -vz*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(0).pz() 
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(0).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(0).py();
 		  
		  pxRest2 = -vx*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(1).px() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(1).py()
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(1).pz();
		  pyRest2 = -vy*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(1).py() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(1).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(1).pz();
		  pzRest2 = -vz*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(1).pz() 
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(1).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(1).py();
 		 		  
		  pxRest3 = -vx*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(2).px() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(2).py()
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(2).pz();
		  pyRest3 = -vy*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(2).py() 
		    + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(2).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(2).pz();
		  pzRest3 = -vz*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(2).pz() 
		    + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(2).px()
		    + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(2).py();
 		  eRest1a = sqrt(pxRest1*pxRest1+pyRest1*pyRest1);
		  eRest2a = sqrt(pxRest2*pxRest2+pyRest2*pyRest2);
		  eRest3a = sqrt(pxRest3*pxRest3+pyRest3*pyRest3);
		 
		  
		  EJetTotala = eRest1a+eRest2a+eRest3a;
		  if (eRest1a>eRest2a && eRest1a>eRest3a)
		    {
		      dalitzXa[0] = 2.*eRest1a/EJetTotala;
		      if ( eRest2a > eRest3a )
			{
			  dalitzXa[1] = 2.*eRest2a/EJetTotala;
			  dalitzXa[2] = 2.*eRest3a/EJetTotala;
			}
		      else 
			{
			  dalitzXa[1] = 2.*eRest3a/EJetTotala;
			  dalitzXa[2] = 2.*eRest2a/EJetTotala;
			}
		    }
		  else if (eRest2a>eRest1a && eRest2a>eRest3a)
		    {
		      dalitzXa[0] = 2.*eRest2a/EJetTotala;
		      if ( eRest1a > eRest3a )
			{
			  dalitzXa[1] = 2.*eRest1a/EJetTotala;
			  dalitzXa[2] = 2.*eRest3a/EJetTotala;
			}
		      else 
			{
			  dalitzXa[1] = 2.*eRest3a/EJetTotala;
			  dalitzXa[2] = 2.*eRest1a/EJetTotala;
			}
		    }
		  else if (eRest3a>eRest1a && eRest3a>eRest2a)
		    {
		      dalitzXa[0] = 2.*eRest3a/EJetTotala;
		      if ( eRest1a > eRest2a )
			{
			  dalitzXa[1] = 2.*eRest1a/EJetTotala;
			  dalitzXa[2] = 2.*eRest2a/EJetTotala;
			}
		      else 
			{
			  dalitzXa[1] = 2.*eRest2a/EJetTotala;
			  dalitzXa[2] = 2.*eRest1a/EJetTotala;
			}
		    }
		  		  		  
		  
		    cout << "pxRest1=" << pxRest1 << endl;
		    cout << "pyRest1=" << pyRest1 << endl;
		    cout << "pzRest1=" << pzRest1 << endl;
		  
		    cout << "sum pxRest=" << pxRest1+pxRest2+pxRest3 << endl;
		    cout << "sum pyRest=" << pyRest1+pyRest2+pyRest3 << endl;
		    cout << "sum pzRest=" << pzRest1+pzRest2+pzRest3 << endl;
		    cout << "eRest1=" << eRest1 << endl;
		    cout << "eRest2=" << eRest2 << endl;
		    cout << "eRest3=" << eRest3 << endl;
		    cout << "eRest1a=" << eRest1a << endl;
		    cout << "eRest2a=" << eRest2a << endl;
		    cout << "eRest3a=" << eRest3a << endl;
		    cout << "e1=" << cellJet.pMassive(0).e() << endl;
		    cout << "e2=" << cellJet.pMassive(1).e() << endl;
		    cout << "e3=" << cellJet.pMassive(2).e() << endl;
		   
		    cout << "eT1=" << cellJet.eT(0) << endl;
		    cout << "eT2=" << cellJet.eT(1) << endl;
		    cout << "eT3=" << cellJet.eT(2) << endl;
		  */
		  
		  posxi = floor(dalitzX[0]*(bins));
		  posyi = floor(dalitzX[1]*(bins));
		  //cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
		  if (abs(cellJet.etaCenter(0))<2. && abs(cellJet.etaCenter(1))<2. && abs(cellJet.etaCenter(2)<2.))
		    {
		      cout << "dalitzx1=" << dalitzX[0] << endl;
		      //if (dalitzX[0]>1.)
		      //{
			  cout << "dalitzx2=" << dalitzX[1] << endl;
			  cout << "dalitzx3=" << dalitzX[2] << endl;
			  //cout << "dalitzx2a=" << dalitzXa[1] << endl;
			  cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
			  //}
		      countAll++;
		      if (dalitzX[0]>0.9 && dalitzX[1]>0.8) countLarge++;
		      //cout << "dalitzx3=" << dalitzX[2] << endl;
		      binning_x1_x2[posxi][posyi]+=1;
		    }
		}

	      p_imax = martini.pythia.event.size();
	      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
		{
		  if (martini.pythia.event[p_i].isFinal()) totalSum++;
		  id = martini.pythia.event[p_i].id();
		  // count pi_0s (111) or pi+ (211)
		  if ( id == 111 ) // pions (pi0)
		    {		
		      p = martini.pythia.event[p_i].pT();                                // p_trans
		      posy = floor(p*(bins/scale));                              // bin number
		      phi = asin(sqrt(pow(martini.pythia.event[p_i].py(),2.))/p);        // azimuthal angle
		      posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
		      //cout << "p_T=" << p << ", px=" << pythia.event[p_i].px() << ", py=" << pythia.event[p_i].py() 
		      //     << ", p_T=" << sqrt(pow(pythia.event[p_i].px(),2.)+pow(pythia.event[p_i].py(),2.)) 
		      //     << ", phi=" << phi << ", phi_deg=" << phi/PI*180. 
		      //     << ", posphi=" << posphi << endl;
		      pl = martini.pythia.event[p_i].pz();                               // p_long
		      En = martini.pythia.event[p_i].e();                                // energy
		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		      if(posy>=0 && posy<bins && abs(y)<=ymax)
			{
			  binning[posy] += 1;
			  pisum += 1;
			  if(posphi>=0 && posphi<phiBins)
			    binning_phi_pt[posphi][posy] += 1;
			}
		    }
		}
	    }
	  else
	    {
	      Event fullEvent;
	      fullEvent.clear();
	      int threeJetsInThisEvent = 0;
	      for (int i=0; i<totalNNs; i++)
		{
		  p_imax = martini.pythia.event.size();
		  for ( int p_i=0; p_i<p_imax; p_i++) 
		    {
		      if (martini.pythia.event[p_i].isFinal()) 
			{
			  fullEvent.append(martini.pythia.event[p_i]);
			}
		    }
		  numEvents+=martini.fragmentation( plist, i );
		  Vec4 boostVec;
		  Vec4 restP;
		  double eRest1, eRest2, eRest3;
		  double pxRest1, pyRest1, pzRest1, pRest1;
		  double pxRest2, pyRest2, pzRest2, pRest2;
		  double pxRest3, pyRest3, pzRest3, pRest3;
		  double px, py, pz;
		  double beta, cosPhi1, cosPhi2, cosPhi3, gamma;
		  double vx, vy, vz;
		  double pl;
		  double En;
		  double y; // rapidity
		  double eta; //pseudo-rapidity
		  double phi; // angle with respect to the reaction plane
		  p_imax = martini.pythia.event.size();
		  
		  // jet finder
		  CellJet cellJet;
		  double dalitzX[3];
		  
		  cellJet.analyze2( martini.pythia.event, 120., 0.38, 2., 5.);
		  // jets will be automatically sorted with decreasing eT
		  if (cellJet.size()>=1)
		    {
		      cout << "found " << cellJet.size() << " jets in NN #" << i << endl;
		      twoJetEvents+=cellJet.size();
		    }
		  if (cellJet.size()==3)
		    {
		      // cout << "________________________" << endl;
		      
		      boostVec=cellJet.pMassive(0)+cellJet.pMassive(1)+cellJet.pMassive(2);
		      //cout << "boostVec=" << boostVec << endl;
		      boostVec/=boostVec.e();
		      vx = boostVec.px();
		      vy = boostVec.py();
		      vz = boostVec.pz();
		      
		      beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of boost velocity in units of c
		      //cout << "beta=" << sqrt(vx*vx+vy*vy+vz*vz) << endl;
		      
		      gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
		      
		      //cout << "beta=" << beta << endl;
		      //cout << "gamma=" << gamma << endl;
		      
		      //pRest = p * gamma * (1.-beta*cosPhi);                     
		      
		      eRest1 = gamma*cellJet.pMassive(0).e()-vx*gamma*cellJet.pMassive(0).px()
			-vy*gamma*cellJet.pMassive(0).py()-vz*gamma*cellJet.pMassive(0).pz();
		      eRest2 = gamma*cellJet.pMassive(1).e()-vx*gamma*cellJet.pMassive(1).px()
			-vy*gamma*cellJet.pMassive(1).py()-vz*gamma*cellJet.pMassive(1).pz();
		      eRest3 = gamma*cellJet.pMassive(2).e()-vx*gamma*cellJet.pMassive(2).px()
			-vy*gamma*cellJet.pMassive(2).py()-vz*gamma*cellJet.pMassive(2).pz();
		      
		      EJetTotal = eRest1+eRest2+eRest3;
		      if (eRest1>eRest2 && eRest1>eRest3)
			{
			  dalitzX[0] = 2.*eRest1/EJetTotal;
			  if ( eRest2 > eRest3 )
			    {
			      dalitzX[1] = 2.*eRest2/EJetTotal;
			      dalitzX[2] = 2.*eRest3/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest3/EJetTotal;
			      dalitzX[2] = 2.*eRest2/EJetTotal;
			    }
			}
		      else if (eRest2>eRest1 && eRest2>eRest3)
			{
			  dalitzX[0] = 2.*eRest2/EJetTotal;
			  if ( eRest1 > eRest3 )
			    {
			      dalitzX[1] = 2.*eRest1/EJetTotal;
			      dalitzX[2] = 2.*eRest3/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest3/EJetTotal;
			      dalitzX[2] = 2.*eRest1/EJetTotal;
			    }
			}
		      else if (eRest3>eRest1 && eRest3>eRest2)
			{
			  dalitzX[0] = 2.*eRest3/EJetTotal;
			  if ( eRest1 > eRest2 )
			    {
			      dalitzX[1] = 2.*eRest1/EJetTotal;
			      dalitzX[2] = 2.*eRest2/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest2/EJetTotal;
			      dalitzX[2] = 2.*eRest1/EJetTotal;
			    }
			}
		      
		      flucMeasure += (sqrt((eRest1-eRest2)*(eRest1-eRest2))+sqrt((eRest1-eRest3)*(eRest1-eRest3))
				      +sqrt((eRest2-eRest3)*(eRest2-eRest3)))/(eRest1+eRest2+eRest3);
		      
		      /*		  
					 pxRest1 = -vx*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(0).px() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(0).py()
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(0).pz();
					 pyRest1 = -vy*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(0).py() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(0).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(0).pz();
					 pzRest1 = -vz*gamma*cellJet.pMassive(0).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(0).pz() 
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(0).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(0).py();
					 
					 pxRest2 = -vx*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(1).px() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(1).py()
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(1).pz();
					 pyRest2 = -vy*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(1).py() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(1).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(1).pz();
					 pzRest2 = -vz*gamma*cellJet.pMassive(1).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(1).pz() 
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(1).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(1).py();
					 
					 pxRest3 = -vx*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vx*vx/(beta*beta))*cellJet.pMassive(2).px() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(2).py()
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(2).pz();
					 pyRest3 = -vy*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vy*vy/(beta*beta))*cellJet.pMassive(2).py() 
					 + (gamma-1.)*vx*vy/(beta*beta)*cellJet.pMassive(2).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(2).pz();
					 pzRest3 = -vz*gamma*cellJet.pMassive(2).e() + (1.+(gamma-1.)*vz*vz/(beta*beta))*cellJet.pMassive(2).pz() 
					 + (gamma-1.)*vx*vz/(beta*beta)*cellJet.pMassive(2).px()
					 + (gamma-1.)*vy*vz/(beta*beta)*cellJet.pMassive(2).py();
					 eRest1a = sqrt(pxRest1*pxRest1+pyRest1*pyRest1);
					 eRest2a = sqrt(pxRest2*pxRest2+pyRest2*pyRest2);
					 eRest3a = sqrt(pxRest3*pxRest3+pyRest3*pyRest3);
					 
					 
					 EJetTotala = eRest1a+eRest2a+eRest3a;
					 if (eRest1a>eRest2a && eRest1a>eRest3a)
					 {
					 dalitzXa[0] = 2.*eRest1a/EJetTotala;
					 if ( eRest2a > eRest3a )
					 {
					 dalitzXa[1] = 2.*eRest2a/EJetTotala;
					 dalitzXa[2] = 2.*eRest3a/EJetTotala;
					 }
					 else 
					 {
					 dalitzXa[1] = 2.*eRest3a/EJetTotala;
					 dalitzXa[2] = 2.*eRest2a/EJetTotala;
					 }
					 }
					 else if (eRest2a>eRest1a && eRest2a>eRest3a)
					 {
					 dalitzXa[0] = 2.*eRest2a/EJetTotala;
					 if ( eRest1a > eRest3a )
					 {
					 dalitzXa[1] = 2.*eRest1a/EJetTotala;
					 dalitzXa[2] = 2.*eRest3a/EJetTotala;
					 }
					 else 
					 {
					 dalitzXa[1] = 2.*eRest3a/EJetTotala;
					 dalitzXa[2] = 2.*eRest1a/EJetTotala;
					 }
					 }
					 else if (eRest3a>eRest1a && eRest3a>eRest2a)
					 {
					 dalitzXa[0] = 2.*eRest3a/EJetTotala;
					 if ( eRest1a > eRest2a )
					 {
					 dalitzXa[1] = 2.*eRest1a/EJetTotala;
					 dalitzXa[2] = 2.*eRest2a/EJetTotala;
					 }
					 else 
					 {
					 dalitzXa[1] = 2.*eRest2a/EJetTotala;
					 dalitzXa[2] = 2.*eRest1a/EJetTotala;
					 }
					 }
					 
					 
					 cout << "pxRest1=" << pxRest1 << endl;
					 cout << "pyRest1=" << pyRest1 << endl;
					 cout << "pzRest1=" << pzRest1 << endl;
					 
					 cout << "sum pxRest=" << pxRest1+pxRest2+pxRest3 << endl;
					 cout << "sum pyRest=" << pyRest1+pyRest2+pyRest3 << endl;
					 cout << "sum pzRest=" << pzRest1+pzRest2+pzRest3 << endl;
					 cout << "eRest1=" << eRest1 << endl;
					 cout << "eRest2=" << eRest2 << endl;
					 cout << "eRest3=" << eRest3 << endl;
					 cout << "eRest1a=" << eRest1a << endl;
					 cout << "eRest2a=" << eRest2a << endl;
					 cout << "eRest3a=" << eRest3a << endl;
					 cout << "e1=" << cellJet.pMassive(0).e() << endl;
					 cout << "e2=" << cellJet.pMassive(1).e() << endl;
					 cout << "e3=" << cellJet.pMassive(2).e() << endl;
					 
					 cout << "eT1=" << cellJet.eT(0) << endl;
					 cout << "eT2=" << cellJet.eT(1) << endl;
					 cout << "eT3=" << cellJet.eT(2) << endl;
		      */
		      

		      posxi = floor(dalitzX[0]*(bins));
		      posyi = floor(dalitzX[1]*(bins));
		      //cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
		      if (abs(cellJet.etaCenter(0))<2. && abs(cellJet.etaCenter(1))<2. && abs(cellJet.etaCenter(2)<2.))
			{
			  threeJetsInThisEvent++;
			  cout << "dalitzx1=" << dalitzX[0] << endl;
			  cout << "eT1=" << cellJet.eT(0) << endl;
			  cout << "eT2=" << cellJet.eT(1) << endl;
			  cout << "eT3=" << cellJet.eT(2) << endl;
			  if (dalitzX[0]>1.)
			    {
			      cout << "dalitzx2=" << dalitzX[1] << endl;
			      cout << "dalitzx3=" << dalitzX[2] << endl;
			      //cout << "dalitzx2a=" << dalitzXa[1] << endl;
			      cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
			    }
			  countAll++;
			  if (dalitzX[0]>0.9 && dalitzX[1]>0.8) countLarge++;
			  //cout << "dalitzx3=" << dalitzX[2] << endl;
			  binning_x1_x2[posxi][posyi]+=1;
			}
		    }
		  
		  for ( int p_i=0; p_i<p_imax; p_i++) 
		    {
		      if (martini.pythia.event[p_i].isFinal()) 
			{
			  totalSum++;
			}
		      id = martini.pythia.event[p_i].id();
		      // count pi_0s (111) or pi+ (211)
		      if ( id == 111 ) // pions (pi0)
			{		
			  p = martini.pythia.event[p_i].pT();                                // p_trans
			  posy = floor(p*(bins/scale));                              // bin number
			  phi = asin(sqrt(pow(martini.pythia.event[p_i].py(),2.))/p);        // azimuthal angle
			  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
			  //cout << "p_T=" << p << ", px=" << pythia.event[p_i].px() << ", py=" << pythia.event[p_i].py() 
			  //     << ", p_T=" << sqrt(pow(pythia.event[p_i].px(),2.)+pow(pythia.event[p_i].py(),2.)) 
			  //     << ", phi=" << phi << ", phi_deg=" << phi/PI*180. 
			  //     << ", posphi=" << posphi << endl;
			  pl = martini.pythia.event[p_i].pz();                               // p_long
			  En = martini.pythia.event[p_i].e();                                // energy
			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
			  if(posy>=0 && posy<bins && abs(y)<=ymax)
			    {
			      binning[posy] += 1;
			      pisum += 1;
			      if(posphi>=0 && posphi<phiBins)
				binning_phi_pt[posphi][posy] += 1;
			    }
			}
		    }
		}
	      cout << "three jets in this event: " << threeJetsInThisEvent << endl;
	       Vec4 boostVec;
		  Vec4 restP;
		  double eRest1, eRest2, eRest3;
		  double px, py, pz;
		  double beta, cosPhi1, cosPhi2, cosPhi3, gamma;
		  double vx, vy, vz;
		  double pl;
		  double En;
		  double y; // rapidity
		  double eta; //pseudo-rapidity
		  double phi; // angle with respect to the reaction plane
		  p_imax = fullEvent.size();
		  
		  // jet finder
		  CellJet cellJet;
		  double dalitzX[3];
		  int relevantJets=0;
		  cellJet.analyze2( fullEvent, 150., 0.4, 2., 10.); //second to last parameter is pTmin for each particle to be accepted
		  //as part of the jet
		  // jets will be automatically sorted with decreasing eT
		  //cout << "number of jets found=" << cellJet.size() << endl;
		  fakeJetEvents+=cellJet.size();
		  for (int i=0; i<cellJet.size(); i++)
		    {
		      if (abs(cellJet.etaCenter(i))<2.)
			relevantJets++;
		      cout << "jet " << i << " eta=" << cellJet.etaCenter(i) << endl;
		      cout << "jet " << i << " phi=" << cellJet.phiCenter(i) << endl;
		      cout << "jet " << i << " eT=" << cellJet.eT(i) << endl;
		    }
		  cout << "number of jets within specified eta range=" << relevantJets << endl;
		  
		  if (cellJet.size()>=3)
		    {
		      // cout << "________________________" << endl;
		      
		      boostVec=cellJet.pMassive(0)+cellJet.pMassive(1)+cellJet.pMassive(2);
		      //cout << "boostVec=" << boostVec << endl;
		      boostVec/=boostVec.e();
		      vx = boostVec.px();
		      vy = boostVec.py();
		      vz = boostVec.pz();
		      
		      beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of boost velocity in units of c
		      //cout << "beta=" << sqrt(vx*vx+vy*vy+vz*vz) << endl;
		      
		      gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
		      
		      //cout << "beta=" << beta << endl;
		      //cout << "gamma=" << gamma << endl;
		      
		      //pRest = p * gamma * (1.-beta*cosPhi);                     
		      
		      eRest1 = gamma*cellJet.pMassive(0).e()-vx*gamma*cellJet.pMassive(0).px()
			-vy*gamma*cellJet.pMassive(0).py()-vz*gamma*cellJet.pMassive(0).pz();
		      eRest2 = gamma*cellJet.pMassive(1).e()-vx*gamma*cellJet.pMassive(1).px()
			-vy*gamma*cellJet.pMassive(1).py()-vz*gamma*cellJet.pMassive(1).pz();
		      eRest3 = gamma*cellJet.pMassive(2).e()-vx*gamma*cellJet.pMassive(2).px()
			-vy*gamma*cellJet.pMassive(2).py()-vz*gamma*cellJet.pMassive(2).pz();
		      
		      EJetTotal = eRest1+eRest2+eRest3;
		      if (eRest1>eRest2 && eRest1>eRest3)
			{
			  dalitzX[0] = 2.*eRest1/EJetTotal;
			  if ( eRest2 > eRest3 )
			    {
			      dalitzX[1] = 2.*eRest2/EJetTotal;
			      dalitzX[2] = 2.*eRest3/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest3/EJetTotal;
			      dalitzX[2] = 2.*eRest2/EJetTotal;
			    }
			}
		      else if (eRest2>eRest1 && eRest2>eRest3)
			{
			  dalitzX[0] = 2.*eRest2/EJetTotal;
			  if ( eRest1 > eRest3 )
			    {
			      dalitzX[1] = 2.*eRest1/EJetTotal;
			      dalitzX[2] = 2.*eRest3/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest3/EJetTotal;
			      dalitzX[2] = 2.*eRest1/EJetTotal;
			    }
			}
		      else if (eRest3>eRest1 && eRest3>eRest2)
			{
			  dalitzX[0] = 2.*eRest3/EJetTotal;
			  if ( eRest1 > eRest2 )
			    {
			      dalitzX[1] = 2.*eRest1/EJetTotal;
			      dalitzX[2] = 2.*eRest2/EJetTotal;
			    }
			  else 
			    {
			      dalitzX[1] = 2.*eRest2/EJetTotal;
			      dalitzX[2] = 2.*eRest1/EJetTotal;
			    }
			}
		      
		      flucMeasure += (sqrt((eRest1-eRest2)*(eRest1-eRest2))+sqrt((eRest1-eRest3)*(eRest1-eRest3))
				      +sqrt((eRest2-eRest3)*(eRest2-eRest3)))/(eRest1+eRest2+eRest3);
		      
		      if (abs(cellJet.etaCenter(0))<2. && abs(cellJet.etaCenter(1))<2. && abs(cellJet.etaCenter(2)<2.))
			{
			  totalJets ++;
			  cout << "in full event: dalitzx1=" << dalitzX[0] << endl;
			  cout << "dalitzx2=" << dalitzX[1] << endl;
			  cout << "dalitzx3=" << dalitzX[2] << endl;
			  cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
			  cout << "ife eT1=" << cellJet.eT(0) << endl;
			  cout << "ife eT2=" << cellJet.eT(1) << endl;
			  cout << "ife eT3=" << cellJet.eT(2) << endl;
			  if (dalitzX[0]>1.)
			    {
			      cout << "dalitzx2=" << dalitzX[1] << endl;
			      cout << "dalitzx3=" << dalitzX[2] << endl;
			      //cout << "dalitzx2a=" << dalitzXa[1] << endl;
			      cout << "dalitz sum=" << dalitzX[0]+dalitzX[1]+dalitzX[2] << endl;
			    }
			}
		    }
		  if (cellJet.size()>2)
		    {
		      for ( int p_i=0; p_i<p_imax; p_i++) 
			{
			  if (fullEvent[p_i].isFinal()) 
			    {	
			      phi = fullEvent[p_i].phi();        // azimuthal angle
			      posphi = floor((phi+PI)*(bins/(2.*PI)));                // bin number in phi
			      En = fullEvent[p_i].pT();                                // energy
			      eta =fullEvent[p_i].eta();          // pseudo-rapidity
			      posyi = floor((eta+2.)*(bins/etascale));
			      //cout << "phi=" << phi << endl;
			      //cout << "eta=" << eta << endl;
			      //cout << "posphi=" << posphi << endl;
			      //cout << "poseta=" << posyi << endl;
			      if (abs(eta)<2.) binning_eta_phi[posyi][posphi]+=En;
			    }
			}
		      break;
		    }

	    }
	}
      
      for(int iy=0; iy<bins; iy++)
	{
	  binning_tmp[iy]=binning[iy]-binning_tmp[iy];
       	  binning_sq[iy]+=binning_tmp[iy]*binning_tmp[iy];
	  binning_tmp[iy]=binning[iy];
	  binning_gamma_tmp[iy]=binning_gamma[iy]-binning_gamma_tmp[iy];
       	  binning_gamma_sq[iy]+=binning_gamma_tmp[iy]*binning_gamma_tmp[iy];
	  binning_gamma_tmp[iy]=binning_gamma[iy];
	  for(int ix=0; ix<phiBins; ix++)
	    {
	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy]-binning_phi_tmp[ix][iy];
	      binning_phi_sq[ix][iy]+=binning_phi_tmp[ix][iy]*binning_phi_tmp[ix][iy];
	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy];
	    }
	}

      if(j%1000==0)
	{
	  cout << "#" << j  << endl; // write every 1000th time step
	  foutt << "#" << j << endl;
	}
    } // end loop over all events


  // ****************************************************************************************** //

  flucMeasure/=static_cast<double>(countAll);
  cout << " ---------------------------------------------- " << endl;
  cout << "total number of real+fake 3-jet events=" << totalJets << endl;
  cout << "number of real 3-jet events=" << countAll << endl;
  cout << "total # of jets=" << fakeJetEvents << endl;
  cout << "actual # of jets=" << twoJetEvents << endl;
  cout << "flucMeasure=" << flucMeasure << endl;
  foutt.close();
      
  // output parton positions:
  fstream foutth("theta.dat",ios::out); 
  foutth.precision(12);  
  foutth << "> Number_of_events = " << runs << endl; 
  foutth << "> quarks " << endl;
  for(int iy=0; iy<bins; iy++)
    {
      foutth << (iy)/(bins/tscale) << " " << bins*static_cast<double>(binning_theta[iy])/tscale << endl; //devide by runs later
    }
  foutth.close();

  fstream foutd("dalitz.dat",ios::out); 
  foutd.precision(12);  
  foutd << "> Dalitz plane = " << runs << endl; 
  foutd << "> x1,x2 " << endl;
  for(int ix=0; ix<bins; ix++)
    for(int iy=0; iy<bins; iy++)
      {
	foutd << (static_cast<double>(ix))/(bins)+1./(2.*bins) << " " << (static_cast<double>(iy))/(bins)+1./(2.*bins)
	      << " " << binning_x1_x2[ix][iy] << endl; //devide by runs later
	if (iy == bins-1) foutd << endl;
      }
  foutd.close();

  fstream foutle("lego.dat",ios::out); 
  foutle.precision(12);  
  foutle << "> lego plot = " << runs << endl; 
  foutle << "> eta, phi " << endl;
  for(int ix=0; ix<bins; ix++)
    for(int iy=0; iy<bins; iy++)
      {
	foutle << (static_cast<double>(ix))/(bins)*etascale+1./(2.*bins/etascale)-2. << " " << (static_cast<double>(iy))/(bins/(2.*PI))
	  +1./(2.*bins/(2.*PI))-PI << " " << binning_eta_phi[ix][iy] << endl; //devide by runs later
	if (iy == bins-1) foutle << endl;
      }
  foutle.close();


  fstream foutx("xy.dat",ios::out); 
  foutx.precision(12);  
  foutx << "> Number_of_events = " << runs << endl; 
  foutx << "> y,x " << endl;
  for(int iy=0; iy<bins; iy++)
    for(int ix=0; ix<bins; ix++)
      {
	foutx << (iy)/(bins/xscale) - 10 << " " << (ix)/(bins/xscale) -10 << " " << bins*bins*binning_x[iy][ix]/xscale/xscale
	      << " " << bins*bins*static_cast<double>(binning_mom[iy][ix])/xscale/xscale
	      << endl; //devide by runs later
	if (ix == bins-1) foutx << endl;
      }
  foutx.close();

  fstream foutr("r.dat",ios::out); 
  foutr.precision(12);  
  foutr << "> Number_of_events = " << runs << endl; 
  foutr << "> r " << endl;
  for(int iy=0; iy<bins; iy++)
    {
      foutr << (iy)/(bins/rscale) << " " << bins*static_cast<double>(binning_r[iy])/rscale  << endl; //divide by runs later
    }
  foutr.close();

  fstream foutl("x.dat",ios::out); 
  foutl.precision(12);  
  foutl << "> Number_of_events = " << runs << endl; 
  foutl << "> x " << endl;
  for(int ix=0; ix<bins; ix++)
    {
      foutl << (ix)/(bins/xscale) - 10 << " " << bins*static_cast<double>(binning_x[ix][20])/xscale  << endl; //divide by runs later
    }
  foutl.close();

  fstream fouty("y.dat",ios::out); 
  fouty.precision(12);  
  fouty << "> Number_of_events = " << runs << endl; 
  fouty << "> y " << endl;
  for(int iy=0; iy<bins; iy++)
    {
      fouty << (iy)/(bins/xscale) - 10 << " " << bins*static_cast<double>(binning_x[20][iy])/xscale  << endl; //divide by runs later
    }
  fouty.close();

  
  pt2r/=(runs*maxTime);
  Etot/=(runs);
  qtTotAll/=runs;
  NColTot/=runs;



  // output partons:
  fstream foutp("partons.dat",ios::out); 
  foutp.precision(12);  
  foutp << "> Number_of_events = " << runs << endl; 
  foutp << "> quarks gluons" << endl;
  for(int iy=0; iy<bins; iy++)
    {
      foutp <<  ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " << bins*static_cast<double>(binning_q[iy])/scale/2. //devide by runs later
	    << " " << bins*static_cast<double>(binning_g[iy])/scale/2. << endl; //devide by runs later
      //divide by 2 since we have 2 initial particles.
    }
  foutp.close();
  //cout << "total number of splittings per run=" << countSplits/runs 
  //	   << ", splittings per parton=" << avSplits/runs << endl;
  //cout << "total number of splittings of the initial partons per run=" << countSplitsOfInitialPartons/runs << endl;
  //cout << "hard gluons with p>500 GeV = " << countHardGluons << endl;
  
  fstream foutphiq("quarkphi.dat",ios::out); 
  foutphiq << "> quarks. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
  foutphiq.precision(12);  
  foutphiq << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<bins; iy++)
    {
      foutphiq << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		<< bins*static_cast<double>(binning_phi_pt_q[0][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[1][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[2][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[3][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[4][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[5][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphiq.close();

  fstream foutphig("gluonphi.dat",ios::out); 
  foutphig << "> gluons. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
  foutphig.precision(12);  
  foutphig << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<bins; iy++)
    {
      foutphig << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		<< bins*static_cast<double>(binning_phi_pt_g[0][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[1][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[2][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[3][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[4][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[5][iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphig.close();

  cout << endl;
  if (martini.returnFragmentationSwitch() == 1)
    {
      fstream fout("pi0.dat",ios::out); 
      fout << "> pi0:" << endl;
      fout.precision(12);  
      fout << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fout << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning[iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
	       << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      fout.close();

      fstream foutphi("pi0phi.dat",ios::out); 
      foutphi << "> pi0. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
      foutphi.precision(12);  
      foutphi << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutphi << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_phi_pt[0][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[1][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[2][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[3][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[4][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[5][iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
	       << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutphi.close();
      
    }  
  cout.precision(6);

  //cout << "K+/pi+=" << static_cast<double>(KPlusSum)/static_cast<double>(piPlusSum) <<endl;
  //cout << "K-/pi-=" << static_cast<double>(KMinusSum)/static_cast<double>(piMinusSum) <<endl;
  cout << "sigmaGen=" << martini.pythia.info.sigmaGen() << endl;
  cout << "pisum=" << pisum << endl;
  cout << "sum of negatively charged hadrons=" << hmsum << endl;
  cout << "used events=" << numEvents << endl;
  cout << "time steps=" << mt << endl;
  cout << "total number of particles=" << totalSum << endl; 
  cout << "total number of partons != u,d,s,g (mainly diquarks with p_t < 1 GeV):" << otherPartons << endl; 
  cout << "qhat =" << pt2r << " GeV^2/fm" << endl;
  cout << "<q_T^2>/mfp =" << qtTotAll/maxTime << " GeV^2/fm" << endl;
  cout << "dE/dx =" << (Ejet-Etot)/maxTime << " GeV/fm" << endl;
  cout << "number of collisions =" << NColTot << endl;
  cout << "partons: " << p_imax << endl;
  cout << "mfp=" << maxTime/NColTot << endl;
  delete plist[0];
  delete plist;

}
