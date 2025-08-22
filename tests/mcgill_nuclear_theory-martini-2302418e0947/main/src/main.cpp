#include "MARTINI.h"

int main(int argc, char* argv[]) 
{
  // create MARTINI object - it will automatically generate PYTHIA object and set its standard values
  // those values can be changed below before doing the explicit PYTHIA init.
  MARTINI martini;

  // how to read in a single variable: (it's just like in PYTHIA)
  // e.g. martini.readString("General:maxTime = 3");
  // all PYTHIA parameters can still be changed by the USER using martini.pythia.readString
  // see the PYTHIA manual for details

  // read in variables from a file - the file can also contain PYTHIA specific changes 
  // MARTINI::readFile will know the difference
  martini.readFile("setup.dat");

  // init MARTINI - sets up PYTHIA too - apply any changes to PYTHIA values AFTER this and BEFORE martini.initPythia
  // as for the moment MARTINI::init() may overwrite your PYTHIA settings if called after them
 
  martini.init(argc, argv);

  /// output all the MARTINI settings
  martini.settings.listAll();

  /// --- do things -------------------------------------------------//

  //output rate using many samples: generateSpectrum(jet energy, process)
  //processes: 1=q->qg, 2=g->qq, 3=g->gg, 4=q->qgamma
  //martini.generateSpectrum(2007.9947553, 4);
  
  //calculate distribution of quark going through brick
  //martini.quarkBrick(10);
  
  // sample thermal distribution as test
  //martini.thermalPlot();

  // test elastic collision part
  //martini.sampleElasticRate(50., 3);
  //martini.sampleElasticRateOmegaQ(500, 300.,1);
  martini.plotHydroDataXYZ();

  /// --- main thing :) -------------------------------------------------//

  // do full brick simulation with events from pythia and possibly fragmentation in the end
  //martini.pythiaEvents();


 // this is my main program - anything can go here
  /*
  Parton jp1;                       // used for fixed energy runs
  Parton jp2;
  
  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plist[0] = new vector<Parton>;    // plist[0] is the main list that holds the high momentum partons that are evolved
  plist[1] = new vector<Parton>;    // plist[1] is auxiliary list holding low mom. partons to be matched when fragmenting
                                    // plist[1] is not used currently 02/24/2009 
  Vec4 vecp, vecp2;                 // pythia four-vector object to store four-momenta
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
  const int bins = 20;              // number of bins
  const int phiBins = 6;            // number of bins
  double scale = 20.;               // maximum p_t in GeV in the binning (bin size = scale/bins)
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
  double mfpTot=0.;
  double qtTotAll=0.;
  double Erun;
  double Ejet=100.;
  double r, theta;
  int posr, postheta, posxi, posyi, posphi;

  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  int runs=martini.returnRuns();
  
  // make them both quarks with x GeV for testing:
  
  if ( martini.returnFixedEnergy() == 1 )
    {
      jp1.id(1); jp2.id(1);
      jp1.p(Ejet,0.,0.); 
      jp2.p(0.,Ejet,0.);
      jp1.col(101); jp1.acol(102); jp2.col(103); jp2.acol(104);
      jp1.x(0.); 
      jp1.y(0.); 
      jp1.z(0.);
      jp1.xini(jp1.x());
      jp1.yini(jp1.y());
      jp1.zini(jp1.z());
      jp2.x(0.); jp2.y(0.); jp2.z(0.);
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


  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  
  ofstream foutt("test.dat",ios::app); 

  for (int j=0; j<runs; j++)      // loop over all events
    {
      martini.setQtTot(0.);
      martini.setNCol(0);
      pt2=0.;
      Erun=0.;
      counter = 0;                // reset counter in the beginning of every event
     

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear auxiliary list (obsolete at the moment)
     
      if( martini.returnFixedEnergy() == 0 )
	{
	  martini.generateEvent(plist[0]); 
	  // version that samples number of collisions with Glauber model
	}
      else
	{
	  plist[0]->push_back( jp1 );
	  //plist[0]->push_back( jp2 );
	}
            
      if (martini.returnEvolution() == 1)         // evolve in medium if settings allow it

	{
	  for(int i=0; i<mt; i++) // loop over all time steps 
	    {
	      counter = martini.evolve(plist, counter, i);
	      counter+=1;
	    }
	}

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
	  
	  if ( pt>5.) // bin initial positions of all partons that are within the intersting y range 
	                              //and above a certain p_t //// abs(y)<=ymax && 
	    {
	      if( (abs(id) > 0 && abs(id) < 4) || id == 21 )
		{
		  if(posxi>0 && posxi<bins && posyi>0 && posyi<bins && abs(y)<=ymax && abs(xini)>0.00000001) 
		    {
		      binning_x[posxi][posyi]+=1./static_cast<double>(p_imax);
		      binning_mom[posxi][posyi]+=plist[0]->at(p_i).pini().px();
		      //cout << "adding one at " << posxi << ", " << posyi << endl;
		    }
		  if(posr>0 && posr<bins) binning_r[posr]+=1;
		  if(postheta>0 && postheta<bins) binning_theta[postheta]+=1;
		}
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
      mfpTot+=martini.returnMfp();
      qtTotAll+=martini.returnQtTot()/p_imax;

      // fragmentation
      if (martini.returnFragmentationSwitch() == 1)
	{

	  numEvents+=martini.fragmentation( plist );
	  
	  double pl;
	  double En;
	  double y; // rapidity
	  double eta; //pseudo-rapidity
	  double phi; // angle with respect to the reaction plane
	  p_imax = martini.pythia.event.size();
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	    {
	      //cout << martini.pythia.event[p_i].id() << " " 
	      //   << martini.pythia.event[p_i].p() << " ";
		
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
	      if ( id == 211 ) // pions (pi+)
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  phi = asin(sqrt(pow(martini.pythia.event[p_i].py(),2.))/p);        // azimuthal angle
		  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_piplus[posy] += 1;
		      piPlusSum += 1;
		    }
		}
	      if ( id == -211 ) // pions (pi-)
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_piminus[posy] += 1;
		      piMinusSum += 1;
		    }
		}
	      if ( id == 321 ) // Kaons (K+)
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_Kplus[posy] += 1;
		      KPlusSum += 1;
		    }
		}
	      if ( id == -321 ) // Kaons (K-)
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_Kminus[posy] += 1;
		      KMinusSum += 1;
		    }
		}
	      if ( id == 22 && ( martini.pythia.event[p_i].status()<90 || martini.pythia.event[p_i].status()>99 ) ) // photons
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_gamma[posy] += 1;
		    }
		}
	      if ( id < 0 && martini.pythia.event[p_i].isHadron() ) // hadrons
		{		
		  p = martini.pythia.event[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = martini.pythia.event[p_i].pz();                               // p_long
		  En = martini.pythia.event[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_hm[posy] += 1;
		      hmsum += 1;
		    }
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
      
      fstream fout2("hminus.dat",ios::out); 
      fout2 << endl << "> h-:" << endl;
      fout2.precision(12);  
      fout2 << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fout2 << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_hm[iy])/scale/(2.*ymax) << endl; 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	}
      fout2.close();

      fstream foutgamma("photons.dat",ios::out); 
      foutgamma << "> gamma:" << endl;
      foutgamma.precision(12);  
      foutgamma << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutgamma << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		    << bins*static_cast<double>(binning_gamma[iy])/scale/(2.*ymax) << " " 
		    << bins*bins*static_cast<double>(binning_gamma_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		    << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutgamma.close();

      fstream foutpiplus("pi+.dat",ios::out); 
      foutpiplus << "> pi+:" << endl;
      foutpiplus.precision(12);  
      foutpiplus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutpiplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_piplus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutpiplus.close();
      fstream foutpiminus("pi-.dat",ios::out); 
      foutpiminus << "> pi-:" << endl;
      foutpiminus.precision(12);  
      foutpiminus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutpiminus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_piminus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutpiminus.close();

      fstream foutKplus("K+.dat",ios::out); 
      foutKplus << "> K+:" << endl;
      foutKplus.precision(12);  
      foutKplus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutKplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_Kplus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutKplus.close();
      fstream foutKminus("K-.dat",ios::out); 
      foutKminus << "> K-:" << endl;
      foutKminus.precision(12);  
      foutKminus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutKminus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_Kminus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutKminus.close();
    }  
  cout.precision(6);

  cout << "K+/pi+=" << static_cast<double>(KPlusSum)/static_cast<double>(piPlusSum) <<endl;
  cout << "K-/pi-=" << static_cast<double>(KMinusSum)/static_cast<double>(piMinusSum) <<endl;
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
  cout << "mfp=" << mfpTot/runs  << " fm" << endl;
  cout << "average distance travelled=" << mfpTot/runs*NColTot  << " fm" << endl;
  delete plist[0];
  delete plist[1];
  delete plist;
  */  
}

