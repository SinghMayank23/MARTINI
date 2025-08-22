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
  //martini.plotHydroDataXYZ();

  /// --- main thing :) -------------------------------------------------//

  // do full brick simulation with events from pythia and possibly fragmentation in the end
  //martini.pythiaEvents();


 // this is my main program - anything can go here

  Parton jp1;                       // used for fixed energy runs
  Parton jp2;
  
  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;    // plist[0] is the main list that holds the high momentum partons that are evolved
  plist[1] = new vector<Parton>;    // plist[1] is not used currently 02/24/2009 
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
  int K0Sum = 0, K0SSum = 0, K0LSum=0;                    
  int hmsum = 0, hpsum=0, pSum=0, pbarSum=0;
  int rhoSum = 0, omegaSum=0, etaSum=0, phiSum=0;
  int eSum = 0, eplusSum=0, JPsiSum=0;
  int muSum = 0, muplusSum = 0;
  const int bins = 40;              // number of bins
  const int xbins = 20;              // number of bins
  const int ecbins = 20;            // number of bins
  const int phiBins = 6;            // number of bins
  const int Lbins = 30;             // number of bins in L binning
  const int corrPhiBins = 48;       // phi bins for correlation measurement
  const int ztBins = 20;            // zT bins for gamma-hadron correlation measurement
  const int pthBins = 20;           // pT^h bins for gamma-hadron correlation measurement

  double scale = 20.;               // maximum p_t in GeV in the binning (bin size = scale/bins)
  double tscale = PI;               // maximum theta
  double rscale = 10.;              // maximum r
  double xscale = 20.;              // maximum x
  double Lscale = 15;               // maximum distance traveled
  double ecscale = 20;
  double ztscale = 1.;              // maximum momentum fraction for gamma-hadron correlations
  double pthscale = 10.;            // maximum momentum of the associated hadron in gamma-hadron correlations


  int binning[bins];                // array with all the bins
  int binning_tmp[bins];            // array with all the bins
  int binning_sq[bins];             // array with all the bins for squares for error

  int binning_forward[bins];                // array with all the bins
  int binning_forward_tmp[bins];            // array with all the bins
  int binning_forward_sq[bins];             // array with all the bins for squares for error

  int binning_ec[bins];             // array with all the bins in number of elastic collisions
  int binning_ec_tmp[bins];         // array with all the bins
  int binning_ec_sq[bins];          // array with all the bins for squares for error

  int binning_mu[bins];             // array with all the bins
  int binning_mu_tmp[bins];         // array with all the bins
  int binning_mu_sq[bins];          // array with all the bins for squares for error

  int binning_muplus[bins];         // array with all the bins
  int binning_muplus_tmp[bins];     // array with all the bins
  int binning_muplus_sq[bins];      // array with all the bins for squares for error

  int binning_JPsi[bins];           // array with all the bins
  int binning_JPsi_tmp[bins];       // array with all the bins
  int binning_JPsi_sq[bins];        // array with all the bins for squares for error

  int binning_e[bins];              // array with all the bins
  int binning_e_tmp[bins];          // array with all the bins
  int binning_e_sq[bins];           // array with all the bins for squares for error

  int binning_eplus[bins];          // array with all the bins
  int binning_eplus_tmp[bins];      // array with all the bins
  int binning_eplus_sq[bins];       // array with all the bins for squares for error

  int binning_phiMeson[bins];       // array with all the bins
  int binning_phiMeson_tmp[bins];   // array with all the bins
  int binning_phiMeson_sq[bins];    // array with all the bins for squares for error

  int binning_eta[bins];            // array with all the bins
  int binning_eta_tmp[bins];        // array with all the bins
  int binning_eta_sq[bins];         // array with all the bins for squares for error

  int binning_omega[bins];          // array with all the bins
  int binning_omega_tmp[bins];      // array with all the bins
  int binning_omega_sq[bins];       // array with all the bins for squares for error

  int binning_rho[bins];            // array with all the bins
  int binning_rho_tmp[bins];        // array with all the bins
  int binning_rho_sq[bins];         // array with all the bins for squares for error
  
  int binning_p[bins];              // array with all the bins
  int binning_p_tmp[bins];          // array with all the bins
  int binning_p_sq[bins];           // array with all the bins for squares for error

  int binning_pbar[bins];           // array with all the bins
  int binning_pbar_tmp[bins];       // array with all the bins
  int binning_pbar_sq[bins];        // array with all the bins for squares for error

  int binning_piplus[bins];         // array with all the bins
  int binning_piplus_tmp[bins];     // array with all the bins for squares for error
  int binning_piplus_sq[bins];      // array with all the bins for squares for error

  int binning_piminus[bins];        // array with all the bins
  int binning_piminus_tmp[bins];    // array with all the bins
  int binning_piminus_sq[bins];     // array with all the bins

  int binning_Kplus[bins];          // array with all the bins
  int binning_Kplus_tmp[bins];      // array with all the bins
  int binning_Kplus_sq[bins];       // array with all the bins

  int binning_Kminus[bins];         // array with all the bins
  int binning_Kminus_tmp[bins];     // array with all the bins
  int binning_Kminus_sq[bins];      // array with all the bins

  int binning_K0S[bins];            // array with all the bins
  int binning_K0S_tmp[bins];        // array with all the bins
  int binning_K0S_sq[bins];         // array with all the bins

  int binning_K0[bins];             // array with all the bins
  int binning_K0_tmp[bins];         // array with all the bins
  int binning_K0_sq[bins];          // array with all the bins

  int binning_K0L[bins];            // array with all the bins
  int binning_K0L_tmp[bins];        // array with all the bins
  int binning_K0L_sq[bins];         // array with all the bins

  int binning_hm[bins];             // array with all the bins
  int binning_hm_tmp[bins];         // array with all the bins
  int binning_hm_sq[bins];          // array with all the bins

  int binning_hp[bins];             // array with all the bins
  int binning_hp_tmp[bins];         // array with all the bins
  int binning_hp_sq[bins];          // array with all the bins

  int binning_g[bins];              // array with all the bins for gluons
  int binning_q[bins];              // array with all the bins for quarks
  int binning_gamma[bins];          // array with all the bins for photons
  int binning_gamma_sq[bins];       // array with all the bins for photons for squares for error
  int binning_gamma_tmp[bins];      // array with all the bins for photons
  int binning_r[bins];              // array with all the bins for the radius in the transverse plane
  int binning_Eini1[bins];          // array with all the bins for the initial energy for all radii in the transverse plane
  double binning_Eini[bins];        // array with all the bins for the initial energy at an initial radius in the transverse plane
  double binning_x[xbins][xbins];     // array with all the bins for the xy positions in the transverse plane

  int binning_corr_phi[corrPhiBins];   // array with bins in phi for correlations
  int binning_corr_phi_sq[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi_tmp[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi2[corrPhiBins];   // array with bins in phi for correlations
  int binning_corr_phi2_sq[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi2_tmp[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi3[corrPhiBins];   // array with bins in phi for correlations
  int binning_corr_phi3_sq[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi3_tmp[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi4[corrPhiBins];   // array with bins in phi for correlations
  int binning_corr_phi4_sq[corrPhiBins];// array with bins in phi for square of correlations
  int binning_corr_phi4_tmp[corrPhiBins];// array with bins in phi for square of correlations

  int binning_corr_zT[ztBins];         // array with bins in zT for correlations
  int binning_corr_zT_sq[ztBins];      // array with bins in zT for square of correlations
  int binning_corr_zT_tmp[ztBins];     // array with bins in zT for square of correlations

  int binning_corr_pTh[ztBins];         // array with bins in pTh for correlations
  int binning_corr_pTh_sq[ztBins];      // array with bins in pTh for square of correlations
  int binning_corr_pTh_tmp[ztBins];     // array with bins in pTh for square of correlations

  int binning_corr_zT2[ztBins];         // array with bins in zT for correlations
  int binning_corr_zT2_sq[ztBins];      // array with bins in zT for square of correlations
  int binning_corr_zT2_tmp[ztBins];     // array with bins in zT for square of correlations

  int binning_corr_pTh2[ztBins];         // array with bins in pTh for correlations
  int binning_corr_pTh2_sq[ztBins];      // array with bins in pTh for square of correlations
  int binning_corr_pTh2_tmp[ztBins];     // array with bins in pTh for square of correlations

  int binning_corr_zT3[ztBins];         // array with bins in zT for correlations
  int binning_corr_zT3_sq[ztBins];      // array with bins in zT for square of correlations
  int binning_corr_zT3_tmp[ztBins];     // array with bins in zT for square of correlations

  int binning_corr_pTh3[ztBins];         // array with bins in pTh for correlations
  int binning_corr_pTh3_sq[ztBins];      // array with bins in pTh for square of correlations
  int binning_corr_pTh3_tmp[ztBins];     // array with bins in pTh for square of correlations

  int binning_phi_pt[phiBins][bins];   // array with bins in phi and pt
  int binning_phi_tmp[phiBins][bins];  // array with bins in phi and pt
  int binning_phi_sq[phiBins][bins];   // array with bins in phi and pt for squares for error

  int binning_phi_forward_pt[phiBins][bins];   // array with bins in phi and pt
  int binning_phi_forward_tmp[phiBins][bins];  // array with bins in phi and pt
  int binning_phi_forward_sq[phiBins][bins];   // array with bins in phi and pt for squares for error

  int binning_phi_pt_g[phiBins][bins]; // array with bins in phi and pt
  int binning_phi_tmp_g[phiBins][bins];// array with bins in phi and pt
  int binning_phi_sq_g[phiBins][bins]; // array with bins in phi and pt for squares for error
  int binning_phi_pt_q[phiBins][bins]; // array with bins in phi and pt
  int binning_phi_tmp_q[phiBins][bins];// array with bins in phi and pt
  int binning_phi_sq_q[phiBins][bins]; // array with bins in phi and pt for squares for error


  double binning_mom[xbins][xbins];   
  int binning_distance_q[Lbins];        // array with bins for the distance a quark has traveled
  int binning_distance_g[Lbins];        // array with bins for the distance a gluon has traveled
  int binning_distance_q_sq[Lbins];    // array with bins for the square of the distance a quark has traveled
  int binning_distance_g_sq[Lbins];    // array with bins for the square of the distance a gluon has traveled
  int binning_distance_q_tmp[Lbins];    // auxiliary array with bins for the distance a quark has traveled
  int binning_distance_g_tmp[Lbins];    // auxiliary array with bins for the distance a gluon has traveled
  int binning_distance_q_forward[Lbins];        // array with bins for the distance a quark has traveled
  int binning_distance_g_forward[Lbins];        // array with bins for the distance a gluon has traveled
  int binning_distance_q_sq_forward[Lbins];    // array with bins for the square of the distance a quark has traveled
  int binning_distance_g_sq_forward[Lbins];    // array with bins for the square of the distance a gluon has traveled
  int binning_distance_q_tmp_forward[Lbins];    // auxiliary array with bins for the distance a quark has traveled
  int binning_distance_g_tmp_forward[Lbins];    // auxiliary array with bins for the distance a gluon has traveled
  int binning_phi_distance_g[phiBins][Lbins]; // array with bins in phi and L
  int binning_phi_distance_tmp_g[phiBins][Lbins];// array with bins in phi and L
  int binning_phi_distance_sq_g[phiBins][Lbins]; // array with bins in phi and L for squares for error
  int binning_phi_distance_q[phiBins][Lbins]; // array with bins in phi and L
  int binning_phi_distance_tmp_q[phiBins][Lbins];// array with bins in phi and L
  int binning_phi_distance_sq_q[phiBins][Lbins]; // array with bins in phi and L for squares for error
  int binning_phi_distance_g_forward[phiBins][Lbins]; // array with bins in phi and L
  int binning_phi_distance_tmp_g_forward[phiBins][Lbins];// array with bins in phi and L
  int binning_phi_distance_sq_g_forward[phiBins][Lbins]; // array with bins in phi and L for squares for error
  int binning_phi_distance_q_forward[phiBins][Lbins]; // array with bins in phi and L
  int binning_phi_distance_tmp_q_forward[phiBins][Lbins];// array with bins in phi and L
  int binning_phi_distance_sq_q_forward[phiBins][Lbins]; // array with bins in phi and L for squares for error
  int binning_theta[bins];             // array with all the bins for the angle theta

  int totalSum = 0;
  int otherPartons = 0;                
  double ymaxCorr = 0.7; 
  double ymaxCorr2 = 1.0;
  double ymaxCorr3 = 0.35;
  double ymax = 0.5; //0.5 for one unit of rapidity around y=0.
  double ymaxPartons = 0.35; 
  double ymaxGammaCorr = 0.35; 
  double ymaxGammaCorr3 = 1.; 
  double pt2=0.;
  double pt2r=0.;
  double Etot=0;
  double NColTot=0.;
  double mfpTot=0.;
  double qtTotAll=0.;
  double Erun;
  double r, theta;
  double projector[3];
  double plong[3], ptrans[3];
  double pTransToInitial, pLongToInitial;
  int totalNNs;
  int posr, postheta, posxi, posyi, posphi, posphi2, posL, posZT, pospTh;
  int totalNumber_q;
  int totalNumber_g;
  int totalNumber_phi_q[phiBins];
  int totalNumber_phi_g[phiBins];
  int totalNumber_q_forward;
  int totalNumber_g_forward;
  int totalNumber_phi_q_forward[phiBins];
  int totalNumber_phi_g_forward[phiBins];
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  int runs=martini.returnRuns();
  double ptsq, dE, qhat, length, pini, ptini;
  int partonsInQhat;
  double qhatTotal = 0;
  double dETotal = 0;
  double ptsqTotal = 0;
  double dEdx, dEdxTotal = 0;
  int totalElasticCollisions = 0;
  int partonsInElasticCollisions = 0;
  int Ntrigger=0; // number of triggers used
  int NPhotontrigger=0; // number of gamma triggers used 
  int NPhotontrigger2=0; // number of gamma triggers used
  int NPhotontrigger3=0; // number of gamma triggers used
  int Ntrigger2=0; // number of gamma triggers used
  int Ntrigger3=0; // number of gamma triggers used
  int Ntrigger4=0; // number of gamma triggers used
  Event fullEvent;
  int initGluons;
  int initQuarks;
  double gluonFraction; // the fraction of gluons for my source term
  gluonFraction = 0.25;
  Random *random;
  random  = new Random();
  double rn;
  int Sid;
  double binning_history[20][static_cast<int>(20./dtfm)]; // bin the values of position and energy loss in bins of time
  int post;
  int elasticCollisions;
  int posElasticCollisions;
  double pl,pt, En, y, xini, yini;
  double phi;
  double L;
  double xi, yi, zi, xf, yf, zf;
  double py,px;
  double pxSq, pySq, pySqTrig;
  double pxAll, pyAll, pxAllTot, pyAllTot, pxAllTrig, pyAllTrig, pxAllTotTrig, pyAllTotTrig, pyTot, pyTotTrig;
  double pxSqSum, pxSqSumTrig, pxTot, pxTotTrig;
  int entriesForPy=0;
  int entriesForPyTrig=0;
  int ctrig=0;
  int ctrigTrig=0;
  int trig, trigTrig;
  pySq=0.;
  pxSqSum =0.;
  pySqTrig=0.;
  pxSqSumTrig =0.;
  // make them both quarks with x GeV for testing:
  
  double pxjet=martini.returnInitialPXjet();
  double pyjet=martini.returnInitialPYjet();
 
  if ( martini.returnFixedEnergy() == 1 )
    {
      jp1.id(1); // check below if this is changed in every run
      jp1.p(pxjet,pyjet,0.); 
      jp1.tini(0.);                    // set initial time
      jp1.pini(jp1.p());                    // set initial momentum
      jp1.col(101); jp1.acol(102);
      jp1.x(martini.returnInitialXjet()); 
      jp1.y(martini.returnInitialYjet()); 
      jp1.z(0.);
      jp1.xini(jp1.x());
      jp1.yini(jp1.y());
      jp1.zini(jp1.z());
      //jp1.splits(0); jp2.splits(0);
      cout << "x_0=" << jp1.x() << ", y_0=" << jp1.y() << ", px_0=" << pxjet << ", py_0=" << pyjet << endl;
    }

  // init binning

  initGluons = 0;
  initQuarks =0;

  for(int iy=0; iy<bins; iy++)
    {
      binning_r[iy]=0;
      binning_Eini[iy]=0;
      binning_Eini1[iy]=0;
      binning_theta[iy]=0;

      binning[iy]=0;
      binning_tmp[iy]=0;
      binning_sq[iy]=0;

      binning_forward[iy]=0;
      binning_forward_tmp[iy]=0;
      binning_forward_sq[iy]=0;

      binning_ec[iy]=0;
      binning_ec_tmp[iy]=0;
      binning_ec_sq[iy]=0;

      binning_mu[iy]=0;
      binning_mu_tmp[iy]=0;
      binning_mu_sq[iy]=0;

      binning_muplus[iy]=0;
      binning_muplus_tmp[iy]=0;
      binning_muplus_sq[iy]=0;

      binning_JPsi[iy]=0;
      binning_JPsi_tmp[iy]=0;
      binning_JPsi_sq[iy]=0;

      binning_e[iy]=0;
      binning_e_tmp[iy]=0;
      binning_e_sq[iy]=0;

      binning_eplus[iy]=0;
      binning_eplus_tmp[iy]=0;
      binning_eplus_sq[iy]=0;

      binning_phiMeson[iy]=0;
      binning_phiMeson_tmp[iy]=0;
      binning_phiMeson_sq[iy]=0;

      binning_eta[iy]=0;
      binning_eta_tmp[iy]=0;
      binning_eta_sq[iy]=0;

      binning_omega[iy]=0;
      binning_omega_tmp[iy]=0;
      binning_omega_sq[iy]=0;

      binning_rho[iy]=0;
      binning_rho_tmp[iy]=0;
      binning_rho_sq[iy]=0;

      binning_p[iy]=0;
      binning_p_tmp[iy]=0;
      binning_p_sq[iy]=0;

      binning_pbar[iy]=0;
      binning_pbar_tmp[iy]=0;
      binning_pbar_sq[iy]=0;

      binning_piplus[iy]=0;
      binning_piplus_tmp[iy]=0;
      binning_piplus_sq[iy]=0;

      binning_Kplus[iy]=0;
      binning_Kplus_tmp[iy]=0;
      binning_Kplus_sq[iy]=0;

      binning_piminus[iy]=0;
      binning_piminus_tmp[iy]=0;
      binning_piminus_sq[iy]=0;

      binning_Kminus[iy]=0;
      binning_Kminus_tmp[iy]=0;
      binning_Kminus_sq[iy]=0;

      binning_hm[iy]=0;
      binning_hm_tmp[iy]=0;
      binning_hm_sq[iy]=0;

      binning_hp[iy]=0;
      binning_hp_tmp[iy]=0;
      binning_hp_sq[iy]=0;

      binning_K0[iy]=0;
      binning_K0_tmp[iy]=0;
      binning_K0_sq[iy]=0;

      binning_K0S[iy]=0;
      binning_K0S_tmp[iy]=0;
      binning_K0S_sq[iy]=0;

      binning_K0L[iy]=0;
      binning_K0L_tmp[iy]=0;
      binning_K0L_sq[iy]=0;

      binning_q[iy]=0;
      binning_g[iy]=0;
      binning_gamma[iy]=0;
      binning_gamma_sq[iy]=0;
      binning_gamma_tmp[iy]=0;
    }

  if(martini.returnTrackHistory())
    {
      for(int iy=0; iy<static_cast<int>(20/dtfm); iy++)
	{
	  for (int ind=0; ind<20; ind++)
	    binning_history[ind][iy]=0;
	}
    }
 
  for(int iy=0; iy<Lbins; iy++)
    {
      binning_distance_q[iy]=0;
      binning_distance_g[iy]=0;
      binning_distance_q_tmp[iy]=0;
      binning_distance_g_tmp[iy]=0;
      binning_distance_q_sq[iy]=0;
      binning_distance_g_sq[iy]=0;
      binning_distance_q_forward[iy]=0;
      binning_distance_g_forward[iy]=0;
      binning_distance_q_tmp_forward[iy]=0;
      binning_distance_g_tmp_forward[iy]=0;
      binning_distance_q_sq_forward[iy]=0;
      binning_distance_g_sq_forward[iy]=0;
    }

  for(int iy=0; iy<xbins; iy++)
    for(int ix=0; ix<xbins; ix++)
      {
	binning_x[iy][ix]=0.;
	binning_mom[iy][ix]=0.;
      }

  for(int ix=0; ix<corrPhiBins; ix++)
    {
      binning_corr_phi[ix]=0;
      binning_corr_phi_sq[ix]=0;
      binning_corr_phi_tmp[ix]=0;
      binning_corr_phi2[ix]=0;
      binning_corr_phi2_sq[ix]=0;
      binning_corr_phi2_tmp[ix]=0;
      binning_corr_phi3[ix]=0;
      binning_corr_phi3_sq[ix]=0;
      binning_corr_phi3_tmp[ix]=0;
      binning_corr_phi4[ix]=0;
      binning_corr_phi4_sq[ix]=0;
      binning_corr_phi4_tmp[ix]=0;
    }

  for(int ix=0; ix<ztBins; ix++)
    {
      binning_corr_zT[ix]=0;
      binning_corr_zT_sq[ix]=0;
      binning_corr_zT_tmp[ix]=0;
      binning_corr_zT2[ix]=0;
      binning_corr_zT2_sq[ix]=0;
      binning_corr_zT2_tmp[ix]=0;
      binning_corr_zT3[ix]=0;
      binning_corr_zT3_sq[ix]=0;
      binning_corr_zT3_tmp[ix]=0;
    }

  for(int ix=0; ix<pthBins; ix++)
    {
      binning_corr_pTh[ix]=0;
      binning_corr_pTh_sq[ix]=0;
      binning_corr_pTh_tmp[ix]=0;
      binning_corr_pTh2[ix]=0;
      binning_corr_pTh2_sq[ix]=0;
      binning_corr_pTh2_tmp[ix]=0;
      binning_corr_pTh3[ix]=0;
      binning_corr_pTh3_sq[ix]=0;
      binning_corr_pTh3_tmp[ix]=0;
    }

  for(int ix=0; ix<phiBins; ix++)
    {
      for(int iy=0; iy<bins; iy++)
	{
	  binning_phi_sq[ix][iy]=0.;
	  binning_phi_pt[ix][iy]=0.;
	  binning_phi_tmp[ix][iy]=0.;
	  binning_phi_forward_sq[ix][iy]=0.;
	  binning_phi_forward_pt[ix][iy]=0.;
	  binning_phi_forward_tmp[ix][iy]=0.;
	  binning_phi_sq_q[ix][iy]=0.;
	  binning_phi_pt_q[ix][iy]=0.;
	  binning_phi_tmp_q[ix][iy]=0.;
	  binning_phi_sq_g[ix][iy]=0.;
	  binning_phi_pt_g[ix][iy]=0.;
	  binning_phi_tmp_g[ix][iy]=0.;
	}
      for(int iy=0; iy<Lbins; iy++)
	{
	  binning_phi_distance_sq_q[ix][iy]=0.;
	  binning_phi_distance_q[ix][iy]=0.;
	  binning_phi_distance_tmp_q[ix][iy]=0.;
	  binning_phi_distance_sq_g[ix][iy]=0.;
	  binning_phi_distance_g[ix][iy]=0.;
	  binning_phi_distance_tmp_g[ix][iy]=0.;
	  binning_phi_distance_sq_q_forward[ix][iy]=0.;
	  binning_phi_distance_q_forward[ix][iy]=0.;
	  binning_phi_distance_tmp_q_forward[ix][iy]=0.;
	  binning_phi_distance_sq_g_forward[ix][iy]=0.;
	  binning_phi_distance_g_forward[ix][iy]=0.;
	  binning_phi_distance_tmp_g_forward[ix][iy]=0.;
	}
      totalNumber_phi_q[ix] = 0;
      totalNumber_phi_g[ix] = 0;
      totalNumber_phi_q_forward[ix] = 0;
      totalNumber_phi_g_forward[ix] = 0;
    }
  
  totalNumber_q = 0;
  totalNumber_g = 0;
  totalNumber_q_forward = 0;
  totalNumber_g_forward = 0;
  

  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  
  ofstream foutt("test.dat",ios::app); 

  partonsInQhat = 0;

  for (int j=0; j<runs; j++)      // loop over all events
    {
      plistInitial->clear();
      martini.setQtTot(0.);
      martini.setNCol(0);
      pt2=0.;
      Erun=0.;
      qhat = 0.;
      dEdx = 0.;
      ptsq = 0.;
      dE = 0.;
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
	  rn=random->genrand64_real1();
	  //	  cout << rn << endl;
	  if (rn>gluonFraction) 
	    Sid=1;
	  else
	    Sid=21;
	  jp1.id(Sid);
	  plist[0]->push_back( jp1 );
	  //plist[0]->push_back( jp2 );
	}
      
      if (martini.returnEvolution() == 1)         // evolve in medium if settings want that
	{
	  for(int i=0; i<mt; i++) // loop over all time steps 
	    {
	      counter = martini.evolve(plist, counter, i);
	      counter+=1;
	    }
	}
      if(martini.returnTrackHistory() && martini.returnFixedEnergy() == 1)
	{
	  //  cout << "size of array=" << martini.Tt->size() << endl;
	  for (int ni=0; ni<martini.Tt->size(); ni++)
	    {
	      //    cout << martini.Tt->at(ni) << " " << martini.Tx->at(ni) << " " << martini.Ty->at(ni) << " " << martini.Tz->at(ni) 
	      //     << " " << martini.TdEdt->at(ni) << " " << martini.Tdpxdt->at(ni) << " " << martini.Tdpydt->at(ni) 
	      //     << " " << martini.Tdpzdt->at(ni) << endl;
	      post = floor(martini.Tt->at(ni)/dtfm);
	      if(post>0 && post<static_cast<int>(20./dtfm))
		{
		  binning_history[13][post]+=martini.Tt->at(ni);
		  binning_history[0][post]+=martini.Tx->at(ni);
		  binning_history[1][post]+=martini.Ty->at(ni);
		  binning_history[14][post]+=martini.Tx->at(ni)*martini.Tx->at(ni);
		  binning_history[15][post]+=martini.Ty->at(ni)*martini.Ty->at(ni);
		  binning_history[2][post]+=martini.Tz->at(ni);
		  binning_history[3][post]+=martini.TdEdt->at(ni);
		  binning_history[4][post]+=martini.Tdpxdt->at(ni);
		  binning_history[5][post]+=martini.Tdpydt->at(ni);
		  binning_history[6][post]+=martini.Tdpzdt->at(ni);
		  binning_history[7][post]+=1;// counts the entries
		  binning_history[8][post]+=martini.TE->at(ni);
		  binning_history[9][post]+=martini.Tpx->at(ni);
		  binning_history[10][post]+=martini.Tpy->at(ni);
		  binning_history[11][post]+=martini.TQGPfrac->at(ni);
		}
	    }
	  martini.Tt->clear();
	  martini.Tx->clear();
	  martini.Ty->clear();
	  martini.Tz->clear();
	  martini.TE->clear();
	  martini.TQGPfrac->clear();
	  martini.Tpx->clear();
	  martini.Tpy->clear();
	  martini.TdEdt->clear();
	  martini.Tdpxdt->clear();
	  martini.Tdpydt->clear();
	  martini.Tdpzdt->clear();
	}
      
      // bin partons:
      p_imax = plist[0]->size();
      
      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	{
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
	  
	  elasticCollisions = plist[0]->at(p_i).elasticCollisions();
	  posElasticCollisions = floor(elasticCollisions*(ecbins/ecscale));

	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
	  
  	  
	  if ( pt>4. && pt<6.66 && abs(y)<=0.35 ) // bin initial positions of all partons that are within the intersting y range 
	                              //and above a certain p_t //// abs(y)<=ymax && 
	    {
	      if( (abs(id) > 0 && abs(id) < 4) || id == 21 )
		{
		  // if(posxi>0 && posxi<bins && posyi>0 && posyi<bins && abs(y)<=0.35 && abs(xini)>0.00000001) 
		  //  {
		  //    binning_x[posxi][posyi]+=1./static_cast<double>(p_imax);
		  //    binning_mom[posxi][posyi]+=plist[0]->at(p_i).pini().pAbs(); // initial energy of parton
		      //cout << "adding one at " << posxi << ", " << posyi << endl;
		  //  }
		  //if(posr>0 && posr<bins) binning_r[posr]+=1;
		  //if(posr>0 && posr<bins) binning_Eini[posr]+=plist[0]->at(p_i).pini().pAbs();
		  if(postheta>0 && postheta<bins) binning_theta[postheta]+=1;
		  if(posElasticCollisions>=0 && posElasticCollisions<ecbins) binning_ec[posElasticCollisions]+=1;
		  totalElasticCollisions += elasticCollisions;
		  partonsInElasticCollisions ++;
		}
	    }
	}
  
      p_imax = plist[0]->size();
      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	{
	  plistInitial->push_back(plist[0]->at(p_i));
	  id = plist[0]->at(p_i).id();
	  p = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
	  projector[0] = plist[0]->at(p_i).pini().px();
	  projector[1] = plist[0]->at(p_i).pini().py();
	  projector[2] = plist[0]->at(p_i).pini().pz();
	 
	  pini = sqrt(projector[0]*projector[0]+projector[1]*projector[1]+projector[2]*projector[2]);
	  ptini = sqrt(projector[0]*projector[0]+projector[1]*projector[1]);
	 	  
	  projector[0]/=pini;
	  projector[1]/=pini;
	  projector[2]/=pini;
	  
	  //cout << "projector=" << projector[0] << " " <<  projector[1] << " " <<  projector[2] << endl; 

	  // momentum component parallel to initial momentum:
	  pLongToInitial = projector[0]*plist[0]->at(p_i).p().px()
	    +projector[1]*plist[0]->at(p_i).p().py()
	    +projector[2]*plist[0]->at(p_i).p().pz();

	  //cout << "pl=" << pl << endl;

	  plong[0] = pLongToInitial*projector[0];
	  plong[1] = pLongToInitial*projector[1];
	  plong[2] = pLongToInitial*projector[2];
	  
	  ptrans[0] =  plist[0]->at(p_i).p().px()-plong[0];
	  ptrans[1] =  plist[0]->at(p_i).p().py()-plong[1];
	  ptrans[2] =  plist[0]->at(p_i).p().pz()-plong[2];

 	  pTransToInitial = sqrt(ptrans[0]*ptrans[0]+ptrans[1]*ptrans[1]+ptrans[2]*ptrans[2]); // transverse to initial momentum

	  //cout << "p=" << plist[0]->at(p_i).p().px() << " " <<  plist[0]->at(p_i).p().py() << " " << plist[0]->at(p_i).p().pz() << endl;
	  //cout << "p_long=" << plong[0] << " " << plong[1] << " " << plong[2] << endl;
	  //cout << "p_trans=" << ptrans[0] << " " << ptrans[1] << " " << ptrans[2] << endl;

	  pl = plist[0]->at(p_i).p().pz();
	  pt = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)); // transverse to beam axis


	  Erun += sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
	  pt2 += (pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
	  posy = floor(pt*(bins/scale));
	  length = plist[0]->at(p_i).tFinal()-plist[0]->at(p_i).tini();
	  //cout << "tfinal=" <<  plist[0]->at(p_i).tFinal() << ", tini=" <<  plist[0]->at(p_i).tini() << endl;
	 
	  En = sqrt(pt*pt+pl*pl);
	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
	  phi = asin(sqrt(pow(plist[0]->at(p_i).p().py(),2.))/pt);   // azimuthal angle
	  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
	  
	  xi=plist[0]->at(p_i).xini();
	  yi=plist[0]->at(p_i).yini();
	  zi=plist[0]->at(p_i).zini();

	  xf=plist[0]->at(p_i).x();
	  yf=plist[0]->at(p_i).y();
	  zf=plist[0]->at(p_i).z();


	  L = sqrt((xf-xi)*(xf-xi)+(yf-yi)*(yf-yi)+(zf-zi)*(zf-zi));
	  
	  posL = floor(L*(Lbins/Lscale));                            // bin the distance traveled
	  //	  if (pt>3.) cout << "L=" << L << ", posL=" << posL << ", posL/(Lbins/Lscale)=" << posL/(Lbins/Lscale) << endl;

// 	  if (L<0.01) 
// 	    {
// 	      cout << "for L=" << L << " we have" << endl;
// 	      cout << "x_f=" << plist[0]->at(p_i).x() << ", x_i=" <<  plist[0]->at(p_i).xini() << endl;
// 	      cout << "y_f=" << plist[0]->at(p_i).y() << ", y_i=" <<  plist[0]->at(p_i).yini() << endl;
// 	      cout << "z_f=" << plist[0]->at(p_i).z() << ", z_i=" <<  plist[0]->at(p_i).zini() << endl;
// 	    }
	  //cout << "p=" << p << endl;
	  //countSplits+=plist[0]->at(p_i).splits();
	  //avSplits+=plist[0]->at(p_i).splits()/plist[0]->size();
	  //if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[0]->at(p_i).splits();
	  // if (p>500.) countHardGluons++;
	  
	  if (martini.returnFixedEnergy()>2) // don't use a rapidity cut for the fixed energy calculation
	    {
	      if( id > 0 && id < 4 ) if(posy>=0 && posy<bins) binning_q[posy]+=1;
	      if( id == 21 ) if(posy>=0 && posy<bins) binning_g[posy]+=1;
	    }
	  else if ( abs(y)>1.2 && abs(y)<2.2 )
	    {
	      if( id > 0 && id < 4 )
		{
		  if(posL>=0 && posL<Lbins && pt>5.75)
		    {
		      totalNumber_q_forward+=1;
		      binning_distance_q_forward[posL]+=1;
		      if(posphi>=0 && posphi<phiBins)
			{
			  totalNumber_phi_q_forward[posphi]+=1;
			  binning_phi_distance_q_forward[posphi][posL] += 1;
			}
		    }
		}
	      else if( id == 21 )
		{
		  if(posL>=0 && posL<Lbins && pt>5.75)
		    {
		      totalNumber_g_forward+=1;
		      binning_distance_g_forward[posL]+=1;
		      if(posphi>=0 && posphi<phiBins)
			{
			  totalNumber_phi_g_forward[posphi]+=1;
			  binning_phi_distance_g_forward[posphi][posL] += 1;
			}
		    }
		}
	    }
	  else if ( abs(y)<=ymaxPartons )
	    {
	      if( id > 0 && id < 4 )
		{
		  if(posL>=0 && posL<Lbins && pt>5.75)
		    {
		      totalNumber_q+=1;
		      binning_distance_q[posL]+=1;
		      if(posphi>=0 && posphi<phiBins)
			{
			  totalNumber_phi_q[posphi]+=1;
			  binning_phi_distance_q[posphi][posL] += 1;
			}
		      if (length>0.)//&& ptini>5.
			{
			  partonsInQhat++;
			  qhat += pTransToInitial*pTransToInitial/length;
			  ptsq  += pTransToInitial*pTransToInitial;
			  dEdx += (pini-En)/length;
			  dE   += (pini-En);
			  //  cout << "qhat=" << qhat << ", pTransToInitial=" << pTransToInitial 
			  //  << ", distance=" << length << endl;
			  //cout << "dE=" << dE << ", dEdx=" << dEdx 
			  //  << ", distance=" << length << endl;
			}
		    }
		  if(posy>=0 && posy<bins) 
		    {
		      binning_q[posy]+=1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_pt_q[posphi][posy] += 1;
		    }
		}
	      else if( id == 21 )
		{
		  if(posL>=0 && posL<Lbins && pt>5.75)
		    {
		      totalNumber_g+=1;
		      binning_distance_g[posL]+=1;
		      if(posphi>=0 && posphi<phiBins)
			{
			  totalNumber_phi_g[posphi]+=1;
			  binning_phi_distance_g[posphi][posL] += 1;
			}
		      if (length>0.)//&& ptini>5.
			{
			  partonsInQhat++;
			  qhat += pTransToInitial*pTransToInitial/length;
			  ptsq  += pTransToInitial*pTransToInitial;
			  dEdx += (pini-En)/length;
			  dE   += (pini-En);
			  //cout << "qhat=" << qhat << ", pTransToInitial=" << pTransToInitial 
			  //   << ", distance=" << length << endl;
			}
		    }
		  if(posy>=0 && posy<bins) 
		    {
		      binning_g[posy]+=1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_pt_g[posphi][posy] += 1;
		    }
		}
	    }
	} // end loop over all partons
      
      qhatTotal+=qhat;
      dEdxTotal+=dEdx;
      ptsqTotal+=ptsq;
      dETotal+=dE;
      
      for(int ix=0; ix<phiBins; ix++)
	{
	  for(int iy=0; iy<bins; iy++)// phi-pt bins
	    {
	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy]-binning_phi_tmp_q[ix][iy];
	      binning_phi_sq_q[ix][iy]+=binning_phi_tmp_q[ix][iy]*binning_phi_tmp_q[ix][iy];
	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy];
	      
	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy]-binning_phi_tmp_g[ix][iy];
	      binning_phi_sq_g[ix][iy]+=binning_phi_tmp_g[ix][iy]*binning_phi_tmp_g[ix][iy];
	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy];
	    }
	  for(int iy=0; iy<Lbins; iy++)//phi-L bins
	    {
	      binning_phi_distance_tmp_q[ix][iy]=binning_phi_distance_q[ix][iy]-binning_phi_distance_tmp_q[ix][iy];
	      binning_phi_distance_sq_q[ix][iy]+=binning_phi_distance_tmp_q[ix][iy]*binning_phi_distance_tmp_q[ix][iy];
	      binning_phi_distance_tmp_q[ix][iy]=binning_phi_distance_q[ix][iy];
	   
	      binning_phi_distance_tmp_g[ix][iy]=binning_phi_distance_g[ix][iy]-binning_phi_distance_tmp_g[ix][iy];
	      binning_phi_distance_sq_g[ix][iy]+=binning_phi_distance_tmp_g[ix][iy]*binning_phi_distance_tmp_g[ix][iy];
	      binning_phi_distance_tmp_g[ix][iy]=binning_phi_distance_g[ix][iy];

	      binning_phi_distance_tmp_q_forward[ix][iy]=binning_phi_distance_q_forward[ix][iy]-binning_phi_distance_tmp_q_forward[ix][iy];
	      binning_phi_distance_sq_q_forward[ix][iy]+=binning_phi_distance_tmp_q_forward[ix][iy]*binning_phi_distance_tmp_q_forward[ix][iy];
	      binning_phi_distance_tmp_q_forward[ix][iy]=binning_phi_distance_q_forward[ix][iy];
	   
	      binning_phi_distance_tmp_g_forward[ix][iy]=binning_phi_distance_g_forward[ix][iy]-binning_phi_distance_tmp_g_forward[ix][iy];
	      binning_phi_distance_sq_g_forward[ix][iy]+=binning_phi_distance_tmp_g_forward[ix][iy]*binning_phi_distance_tmp_g_forward[ix][iy];
	      binning_phi_distance_tmp_g_forward[ix][iy]=binning_phi_distance_g_forward[ix][iy];
	    }
	}

      for(int iy=0; iy<Lbins; iy++)//L-bins
	{
	  binning_distance_q_tmp[iy]=binning_distance_q[iy]-binning_distance_q_tmp[iy];
       	  binning_distance_q_sq[iy]+=binning_distance_q_tmp[iy]*binning_distance_q_tmp[iy];
	  binning_distance_q_tmp[iy]=binning_distance_q[iy];

	  binning_distance_g_tmp[iy]=binning_distance_g[iy]-binning_distance_g_tmp[iy];
       	  binning_distance_g_sq[iy]+=binning_distance_g_tmp[iy]*binning_distance_g_tmp[iy];
	  binning_distance_g_tmp[iy]=binning_distance_g[iy];

	  binning_distance_q_tmp_forward[iy]=binning_distance_q_forward[iy]-binning_distance_q_tmp_forward[iy];
       	  binning_distance_q_sq_forward[iy]+=binning_distance_q_tmp_forward[iy]*binning_distance_q_tmp_forward[iy];
	  binning_distance_q_tmp_forward[iy]=binning_distance_q_forward[iy];

	  binning_distance_g_tmp_forward[iy]=binning_distance_g_forward[iy]-binning_distance_g_tmp_forward[iy];
       	  binning_distance_g_sq_forward[iy]+=binning_distance_g_tmp_forward[iy]*binning_distance_g_tmp_forward[iy];
	  binning_distance_g_tmp_forward[iy]=binning_distance_g_forward[iy];
	}

     
      pt2r+=pt2/p_imax;
      Etot+=Erun/p_imax;
      NColTot+=static_cast<double>(martini.returnNCol())/static_cast<double>(p_imax);
      mfpTot+=martini.returnMfp();
      qtTotAll+=martini.returnQtTot()/p_imax;
      int numEventsOld=0;
      // fragmentation
      if (martini.returnFragmentationSwitch() == 1)
	{
	  if ( martini.returnFullEvent() == 1 ) 
	    {  
	      fullEvent.clear();
	      martini.pythia.event.clear();
	      totalNNs = martini.returnTotalNNs();
	      for (int i=0; i<totalNNs; i++)
		{
		  numEvents+=martini.fragmentation( plist, i );
		  if (numEvents==numEventsOld) continue;
		  numEventsOld=numEvents;
		  p_imax = martini.pythia.event.size();
		  //cout << "p_imax=" << p_imax<< endl;
		  for ( int p_i=0; p_i<p_imax; p_i++) 
		    {
		      if (martini.pythia.event[p_i].isFinal() && martini.pythia.event[p_i].pT()>0.5 ) //only keep particles with pT>0.5GeV
			{
			  fullEvent.append(martini.pythia.event[p_i]);
			}
		    }
		}
	    }
	  else // if not a full event
	    {
	      numEvents+=martini.fragmentation( plist );
	      fullEvent = martini.pythia.event;
	    } 
	  
	  int gotAssoc;
	  double pl;
	  double En;
	  double y, y2; // rapidity
	  double eta, eta2; //pseudo-rapidity
	  double phi, phi2, deltaPhi; // angle with respect to the reaction plane
	  double pTgamma, pThadron;
	  p_imax = fullEvent.size();
	  //cout << "fullEvent.size=" << fullEvent.size() << endl;
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all particles
	    {
	      //cout << fullEvent[p_i].id() << " " 
	      //   << fullEvent[p_i].p() << " ";
	      
	      if (fullEvent[p_i].isFinal()) totalSum++;
	      id = fullEvent[p_i].id();
	      
	      // compute deltaPhi correlation between charged hadrons
	      if ( fullEvent[p_i].charge() != 0 && fullEvent[p_i].isHadron() && fullEvent[p_i].isFinal() )
		{
		  // if the hadron is a trigger with momentum between ptmin and ptmax (e.g. 4 and 6 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > 4. && fullEvent[p_i].pT() < 6. )
		    {
		      //cout << "trigger p_T=" << fullEvent[p_i].pT() << endl;
		      // trigger is not included in associated because of the p_T restriction pT_as<pT_tr...
		      //gotAssoc = 0;
		      pl = fullEvent[p_i].pz();                               // p_long
		      p = fullEvent[p_i].pT();       
		      eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		      if ( abs(eta)<=ymaxCorr )
			{
			  Ntrigger++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > 2. && fullEvent[p_j].pT() < fullEvent[p_i].pT() 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() )
			      {
			 	
				pl = fullEvent[p_j].pz();                               // p_long
				p = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<=ymaxCorr) //pseudorapidity |eta| < ymaxCorr
				  {
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    posphi = floor((PI+deltaPhi)*(corrPhiBins/(2.*PI)));               // bin number in phi
				    
				    //cout << "assoc. p_T=" << fullEvent[p_j].pT() << endl;
				    
				    //cout << "posphi=" << posphi << endl;
				    
				    //cout << "y_t=" << y << ", eta_t=" << eta << endl;
				    if(posphi>=0 && posphi<=corrPhiBins)
				      {
					binning_corr_phi[posphi] += 1;
					// cout << "phi=" << phi << ", phi2=" << phi2 << ", deltaPhi=" << deltaPhi << endl;
				      }
				  }
			      }		   
			}
		    }
		}
	    

                // find the distribution of initial partons (position and p_T for a given final hadron)
	        double pl,pt,theta, En, eta, xini, yini;
		id = fullEvent[p_i].id();
		pt = fullEvent[p_i].pT();
		pl = fullEvent[p_i].pz(); //p_long
		En = sqrt(pt*pt+pl*pl);
		eta = 0.5*log((sqrt(pl*pl+pt*pt)+pl)/(sqrt(pl*pl+pt*pt)-pl)); //pseudo-rapidity
		int number;
		int posE;
		number=0;
		if ( pt>3. && pt<4. && abs(eta)<=0.35 && (abs(id)==211 || abs(id)==2212 || abs(id)==321) ) 
		  // bin initial positions of all partons that are within the intersting y range 
		  //and above a certain p_t //// abs(y)<=ymax && 
		  {
		       for (int ci=1; ci<plistInitial->size(); ci++)
		      {
			//			cout << " x_ini(" << ci << ")=" << plistInitial->at(ci).xini() << endl;
			if(plistInitial->at(ci).pini().pT()>plistInitial->at(number).pini().pT() 
			   && (abs(plistInitial->at(ci).id())<4||abs(plistInitial->at(ci).id())==21))
			  { 
			    number = ci;
			    //cout << " pt_ini=" << plistInitial->at(ci).pini().pT() << endl;
			  }
		      }
		    r = sqrt(pow(plistInitial->at(number).xini(),2.)+pow(plistInitial->at(number).yini(),2.)); 
		    xini = plistInitial->at(number).xini();
		    yini = plistInitial->at(number).yini();
		    
		    double alpha, phix; 
		    alpha=atan(-fullEvent[p_i].py()/fullEvent[p_i].px());
		    //now rotate the initial position vector so that the hadron points in the -x direction.
		    phix=atan(yini/xini);
		    double xiniNew = cos(alpha)*xini-sin(alpha)*yini;
		    double yiniNew = sin(alpha)*xini+cos(alpha)*yini;
		    double newPx =  cos(alpha)*(fullEvent[p_i].px()+xini)-sin(alpha)*(fullEvent[p_i].py()+yini);
		    double newPy =  sin(alpha)*(fullEvent[p_i].px()+xini)+cos(alpha)*(fullEvent[p_i].py()+yini);
		    newPx-=xiniNew;
		    newPy-=yiniNew;
		    //cout << "xini=" << xini << ", yini=" << yini << endl;
		    if(newPx>0.)
		      {
			newPx*=-1;
			newPy*=-1;
			xiniNew*=-1.;
			yiniNew*=-1.;
		      }
		    cout << "xiniNew=" << xiniNew << ", yiniNew=" << yiniNew << endl;
		    cout << "new px=" << newPx << ", new py=" << newPy << endl;
		    //martini.pythia.event.list();
		    cout << plistInitial->at(number).id() << " pT_ini=" << plistInitial->at(number).pini().pT() 
			 << " sqrt(px^2_ini+p_y^2_ini)=" << sqrt(pow(plistInitial->at(number).pini().px(),2.)+pow(plistInitial->at(number).pini().py(),2.))
			 << " p_T_final=" << pt << " E_final=" << En << " r=" << r << endl;
		    //    sleep(1);
		    xiniNew-=martini.returnHydroTau0();// shift by distance traveled before the evolution started
		    posxi = floor((xiniNew+10)*(xbins/xscale));
		    posyi = floor((yiniNew+10)*(xbins/xscale));
		    posr = floor(r*(bins/rscale));
		    posE = floor(plistInitial->at(number).pini().pT()*(bins/rscale));
		    if(abs(plistInitial->at(number).id())<4 && abs(plistInitial->at(number).id())>0) initQuarks+=1;
		    if(abs(plistInitial->at(number).id()) == 21) initGluons+=1;
		    
		    trig=0;
		    trigTrig=0;

		    for (int i = 1; i< plistInitial->size(); i++)
 		      {
			if (plistInitial->at(i).pini().pz()<0.2 && plistInitial->at(i).pini().pT()>1.5)
			  //compute projection of parton pT on initial hadron direction.
			  {
			    px=(plistInitial->at(i).pini().px()*fullEvent[p_i].px()
				+plistInitial->at(i).pini().py()*fullEvent[p_i].py())/fullEvent[p_i].pT();
			    pxSq=px*px;
			    if(plistInitial->at(i).pini().px()/fullEvent[p_i].px()<0.)
			      {
				pxSqSum += pxSq;
				pySq += pow(plistInitial->at(i).pini().pT(),2.)-pxSq;
				entriesForPy++;
				pxTot = sqrt(pxSqSum/entriesForPy); 
				pyTot = sqrt(pySq/entriesForPy);
				cout << "sqrt(<px^2>)=" << pxTot << endl;
				cout << "sqrt(<py^2>)=" << pyTot << endl;
				if(trig==0)
				  {
				    ctrig++;
				    trig=1;
				  }
				pxAll += pxSq;
				pyAll += pow(plistInitial->at(i).pini().pT(),2.)-pxSq;
				pxAllTot = sqrt(pxAll/ctrig);
				pyAllTot = sqrt(pyAll/ctrig);
				cout << "<Sum px^2>=" << pxAllTot << endl;
				cout << "<Sum py^2>=" << pyAllTot << endl;
				cout <<"ctrig=" << ctrig << endl; 
				cout <<"entriesForPy=" << entriesForPy << endl; 
			      }
			    else
			      {
				pxSqSumTrig += pxSq;
				pySqTrig += pow(plistInitial->at(i).pini().pT(),2.)-pxSq;
				entriesForPyTrig++;
				pxTotTrig = sqrt(pxSqSumTrig/entriesForPyTrig); 
				pyTotTrig = sqrt(pySqTrig/entriesForPyTrig);
				if(trigTrig==0)
				  {
				    ctrigTrig++;
				    trigTrig=1;
				  }
				pxAllTrig += pxSq;
				pyAllTrig += pow(plistInitial->at(i).pini().pT(),2.)-pxSq;
				pxAllTotTrig = sqrt(pxAllTrig/ctrigTrig);
				pyAllTotTrig = sqrt(pyAllTrig/ctrigTrig);
				cout << "sqrt(<px^2>)_Trigger=" << pxTotTrig << endl;
				cout << "sqrt(<py^2>)_Trigger=" << pyTotTrig << endl;
				cout << "<Sum px^2>_Trigger=" << pxAllTotTrig << endl;
				cout << "<Sum py^2>_Trigger=" << pyAllTotTrig << endl;
				cout <<"ctrigTrig=" << ctrigTrig << endl; 
				cout <<"entriesForPyTrig=" << entriesForPyTrig << endl; 
			      }
			  }
		      }
		  		    
		    if(posxi>0 && posxi<xbins && posyi>0 && posyi<xbins) 
		      {
			binning_x[posxi][posyi]+=1;
			binning_mom[posxi][posyi]+=plistInitial->at(number).pini().pT();
			//binning_mom[posxi][posyi]+=plist[0]->at(p_i).pini().pAbs(); // initial energy of parton
			//cout << "adding one at " << posxi << ", " << posyi << endl;
		      }
		    if(posr>0 && posr<bins) binning_r[posr]+=1;
		    if(posr>0 && posr<bins) binning_Eini[posr]+=plistInitial->at(number).pini().pT();
		    if(posE>0 && posE<bins) binning_Eini1[posE]+=1;
		    if(martini.returnTrackHistory() && abs(alpha)<0.1 && fullEvent[p_i].px()<0.)
		      {
			cout << " ********************* adding history" << endl;
			//  cout << "size of array=" << martini.Tt->size() << endl;
			for (int ni=0; ni<martini.Tt->size(); ni++)
			  {
			    //    cout << martini.Tt->at(ni) << " " << martini.Tx->at(ni) << " " << martini.Ty->at(ni) << " " << martini.Tz->at(ni) 
			    //     << " " << martini.TdEdt->at(ni) << " " << martini.Tdpxdt->at(ni) << " " << martini.Tdpydt->at(ni) 
			    //     << " " << martini.Tdpzdt->at(ni) << endl;
			    post = floor(martini.Tt->at(ni)/dtfm+0.000001);
			    if(post>0 && post<static_cast<int>(20./dtfm) && martini.Tpx->at(ni)>0.)
			      {
			// 	cout << "t=" << martini.Tt->at(ni) << ", x=" << martini.Tx->at(ni) 
// 				     << ", y=" << martini.Ty->at(ni)
// 				     << ", z=" << martini.Tz->at(ni) << endl;
				binning_history[13][post]+=martini.Tt->at(ni);
				binning_history[0][post]+=martini.Tx->at(ni);
				binning_history[1][post]+=martini.Ty->at(ni);
				binning_history[2][post]+=martini.Tz->at(ni);
				binning_history[14][post]+=martini.Tx->at(ni)*martini.Tx->at(ni);
				binning_history[15][post]+=martini.Ty->at(ni)*martini.Ty->at(ni);
				binning_history[3][post]+=martini.TdEdt->at(ni);
				binning_history[4][post]+=martini.Tdpxdt->at(ni);
				binning_history[5][post]+=martini.Tdpydt->at(ni);
				binning_history[6][post]+=martini.Tdpzdt->at(ni);
				binning_history[7][post]+=1;// counts the entries
				binning_history[8][post]+=martini.TE->at(ni);
				binning_history[9][post]+=martini.Tpx->at(ni);
				binning_history[10][post]+=martini.Tpy->at(ni);
				binning_history[11][post]+=martini.TQGPfrac->at(ni);
			      }
			  }
			martini.Tt->clear();
			martini.Tx->clear();
			martini.Ty->clear();
			martini.Tz->clear();
			martini.TE->clear();
			martini.TQGPfrac->clear();
			martini.Tpx->clear();
			martini.Tpy->clear();
			martini.TdEdt->clear();
			martini.Tdpxdt->clear();
			martini.Tdpydt->clear();
			martini.Tdpzdt->clear();
		      }
		  }
	      
       	      // compute deltaPhi correlation between charged hadrons 2 (STAR 2)
	      if ( fullEvent[p_i].charge() != 0 && fullEvent[p_i].isHadron() && fullEvent[p_i].isFinal() )
		{
		  // if the hadron is a trigger with momentum between ptmin and ptmax (e.g. 4 and 6 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > 4. && fullEvent[p_i].pT() < 6. )
		    {
		      //cout << "trigger p_T=" << fullEvent[p_i].pT() << endl;
		      // trigger is not included in associated because of the p_T restriction pT_as<pT_tr...
		      //gotAssoc = 0;
		      pl = fullEvent[p_i].pz();                               // p_long
		      p = fullEvent[p_i].pT();       
		      eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		      if ( abs(eta)<=ymaxCorr2 )
			{
			  Ntrigger2++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > 2.5 && fullEvent[p_j].pT() < 4. 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() )
			      {
			 	
				pl = fullEvent[p_j].pz();                               // p_long
				p = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<=ymaxCorr2) //pseudorapidity |eta| < ymaxCorr
				  {
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    posphi = floor((PI+deltaPhi)*(corrPhiBins/(2.*PI)));               // bin number in phi
				    
				    //cout << "assoc. p_T=" << fullEvent[p_j].pT() << endl;
				    
				    //cout << "posphi=" << posphi << endl;
				    
				    //cout << "y_t=" << y << ", eta_t=" << eta << endl;
				    if(posphi>=0 && posphi<=corrPhiBins)
				      {
					binning_corr_phi2[posphi] += 1;
					// cout << "phi=" << phi << ", phi2=" << phi2 << ", deltaPhi=" << deltaPhi << endl;
				      }
				  }
			      }		   
			}
		    }
		}

	      // compute deltaPhi correlation between charged hadrons 3 (PHENIX 1)
	      if ( fullEvent[p_i].charge() != 0 && fullEvent[p_i].isHadron() && fullEvent[p_i].isFinal() )
		{
		  // if the hadron is a trigger with momentum between ptmin and ptmax (e.g. 4 and 6 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > 3. && fullEvent[p_i].pT() < 4. )
		    {
		      //cout << "trigger p_T=" << fullEvent[p_i].pT() << endl;
		      // trigger is not included in associated because of the p_T restriction pT_as<pT_tr...
		      //gotAssoc = 0;
		      pl = fullEvent[p_i].pz();                               // p_long
		      p = fullEvent[p_i].pT();       
		      eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		      if ( abs(eta)<ymaxCorr3 )
			{
			  Ntrigger3++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > 2. && fullEvent[p_j].pT() < 3. 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() )
			      {
			 	
				pl = fullEvent[p_j].pz();                               // p_long
				p = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<ymaxCorr3) //pseudorapidity |eta| < ymaxCorr
				  {
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    posphi = floor((PI+deltaPhi)*(corrPhiBins/(2.*PI)));               // bin number in phi
				    
				    //cout << "assoc. p_T=" << fullEvent[p_j].pT() << endl;
				    
				    //cout << "posphi=" << posphi << endl;
				    
				    //cout << "y_t=" << y << ", eta_t=" << eta << endl;
				    if(posphi>=0 && posphi<=corrPhiBins)
				      {
					binning_corr_phi3[posphi] += 1;
					// cout << "phi=" << phi << ", phi2=" << phi2 << ", deltaPhi=" << deltaPhi << endl;
				      }
				  }
			      }		   
			}
		    }
		}

	      // compute deltaPhi correlation between charged hadrons 4 (PHENIX 2)
	      if ( fullEvent[p_i].charge() != 0 && fullEvent[p_i].isHadron() && fullEvent[p_i].isFinal() )
		{
		  // if the hadron is a trigger with momentum between ptmin and ptmax (e.g. 4 and 6 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > 5. && fullEvent[p_i].pT() < 10. )
		    {
		      //cout << "trigger p_T=" << fullEvent[p_i].pT() << endl;
		      // trigger is not included in associated because of the p_T restriction pT_as<pT_tr...
		      //gotAssoc = 0;
		      pl = fullEvent[p_i].pz();                               // p_long
		      p = fullEvent[p_i].pT();       
		      eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		      if ( abs(eta)<=ymaxCorr3 )
			{
			  Ntrigger4++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > 3. && fullEvent[p_j].pT() < 5. 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() )
			      {
			 	
				pl = fullEvent[p_j].pz();                               // p_long
				p = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<=ymaxCorr3) //pseudorapidity |eta| < ymaxCorr
				  {
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    posphi = floor((PI+deltaPhi)*(corrPhiBins/(2.*PI)));               // bin number in phi
				    
				    //cout << "assoc. p_T=" << fullEvent[p_j].pT() << endl;
				    
				    //cout << "posphi=" << posphi << endl;
				    
				    //cout << "y_t=" << y << ", eta_t=" << eta << endl;
				    if(posphi>=0 && posphi<=corrPhiBins)
				      {
					binning_corr_phi4[posphi] += 1;
					// cout << "phi=" << phi << ", phi2=" << phi2 << ", deltaPhi=" << deltaPhi << endl;
				      }
				  }
			      }		   
			}
		    }
		}

	      // gamma-hadron correlations:
	      double zT;
	      double gammapTmin = 5.;
	      double gammapTmax = 7.;
	      double hadronpTmin = 1.;
	      double hadronpTmax = 7.;
	      if ( fullEvent[p_i].id() == 22 && fullEvent[p_i].isFinal() && ( fullEvent[p_i].status()<90 || fullEvent[p_i].status()>99 ))
		// if we have a direct photon (no decay photons (as the status 90-99 would indicate)
		{
		  // if the photon is a trigger with momentum between ptmin and ptmax (e.g. 8 and 16 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > gammapTmin && fullEvent[p_i].pT() < gammapTmax )
		    {
		      cout << "photon trigger p_T=" << fullEvent[p_i].pT() << endl;
		      pl = fullEvent[p_i].pz();                               // p_long
		      pTgamma = fullEvent[p_i].pT();       
		      y = 0.5*log((En+pl)/(En-pl));   
		      eta = 0.5*log((sqrt(pl*pl+pTgamma*pTgamma)+pl)/(sqrt(pl*pl+pTgamma*pTgamma)-pl)); //pseudo-rapidity
		      
		      if(abs(eta)<ymaxGammaCorr) //rapidity |eta| < ymaxGammaCorr
			{
			  NPhotontrigger++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > hadronpTmin && fullEvent[p_j].pT() < hadronpTmax 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() ) // associated charged hadrons
			      {
			   	
				pl = fullEvent[p_j].pz();                               // p_long
				pThadron = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+pThadron*pThadron)+pl)/(sqrt(pl*pl+pThadron*pThadron)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<ymaxGammaCorr) //rapidity |eta| < ymaxGammaCorr
				  {
				    
				    cout << "assoc.hadron p_T=" << pThadron << endl;
				    
				    zT=pThadron/pTgamma; // momentum fraction	
				    
				    posZT = floor((zT)*(ztBins/(ztscale)));
				    pospTh = floor((pThadron)*(pthBins/(pthscale)));
				     
				    cout << "zT=" << zT << endl;
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    if(posZT>=0 && posZT<=ztBins && abs(deltaPhi-PI)<PI/5.) 
				      {
					binning_corr_zT[posZT] += 1;
				      }
				
				    if(pospTh>=0 && pospTh<=pthBins && abs(deltaPhi-PI)<PI/5.) 
				      {
					binning_corr_pTh[pospTh] += 1;
				      }
				  }
			      }		   
			}
		      if(abs(eta)<ymaxGammaCorr3) //rapidity |eta| < ymaxGammaCorr
			{
			  NPhotontrigger3++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > hadronpTmin && fullEvent[p_j].pT() < hadronpTmax 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() ) // associated charged hadrons
			      {
			   	
				pl = fullEvent[p_j].pz();                               // p_long
				pThadron = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+pThadron*pThadron)+pl)/(sqrt(pl*pl+pThadron*pThadron)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<ymaxGammaCorr3) //rapidity |eta| < ymaxGammaCorr
				  {
				    
				    cout << "assoc.hadron p_T=" << pThadron << endl;
				    
				    zT=pThadron/pTgamma; // momentum fraction	
				    
				    posZT = floor((zT)*(ztBins/(ztscale)));
				    pospTh = floor((pThadron)*(pthBins/(pthscale)));
				     
				    cout << "zT=" << zT << endl;
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    if(posZT>=0 && posZT<=ztBins && abs(deltaPhi-PI)<PI/5.) 
				      {
					binning_corr_zT3[posZT] += 1;
				      }
				
				    if(pospTh>=0 && pospTh<=pthBins && abs(deltaPhi-PI)<PI/5.) 
				      {
					binning_corr_pTh3[pospTh] += 1;
				      }
				  }
			      }		   
			}
		    }
		}
	    

	      // gamma-hadron correlations 2 (another set of cuts):
	      double gammapTmin2 = 8.;
	      double gammapTmax2 = 16.;
	      double hadronpTmin2 = 0.;
	      double hadronpTmax2 = 16.;
	      double ymaxGammaCorr2 = 1.;
	      if ( fullEvent[p_i].id() == 22 && fullEvent[p_i].isFinal() && ( fullEvent[p_i].status()<90 || fullEvent[p_i].status()>99 ))
		// if we have a direct photon (no decay photons (as the status 90-99 would indicate)
		{
		  // if the photon is a trigger with momentum between ptmin and ptmax (e.g. 8 and 16 GeV)
		  //cout << "p_T=" << fullEvent[p_i].pT() << endl;
		  if ( fullEvent[p_i].pT() > gammapTmin2 && fullEvent[p_i].pT() < gammapTmax2 )
		    {
		      cout << "photon trigger p_T=" << fullEvent[p_i].pT() << endl;
		      pl = fullEvent[p_i].pz();                               // p_long
		      pTgamma = fullEvent[p_i].pT();       
		      y = 0.5*log((En+pl)/(En-pl));   
		      eta = 0.5*log((sqrt(pl*pl+pTgamma*pTgamma)+pl)/(sqrt(pl*pl+pTgamma*pTgamma)-pl)); //pseudo-rapidity
		      
		      if(abs(eta)<ymaxGammaCorr2) //rapidity |eta| < ymaxGammaCorr
			{
			  NPhotontrigger2++;
			  for ( int p_j=0; p_j<p_imax; p_j++) // loop over all partons again and find correlations
			    
			    if ( fullEvent[p_j].pT() > hadronpTmin2 && fullEvent[p_j].pT() < hadronpTmax2 
				 && fullEvent[p_j].isHadron() && fullEvent[p_j].charge() != 0 &&
				 fullEvent[p_j].isFinal() ) // associated charged hadrons
			      {
			   	
				pl = fullEvent[p_j].pz();                               // p_long
				pThadron = fullEvent[p_j].pT();       
				eta2 = 0.5*log((sqrt(pl*pl+pThadron*pThadron)+pl)/(sqrt(pl*pl+pThadron*pThadron)-pl)); //pseudo-rapidity
				
				if(abs(eta2)<ymaxGammaCorr2) //rapidity |eta| < ymaxGammaCorr
				  {
				    
				    cout << "assoc.hadron p_T=" << pThadron << endl;
				    
				    zT=pThadron/pTgamma; // momentum fraction	
				    
				    posZT = floor((zT)*(ztBins/(ztscale)));
				    pospTh = floor((pThadron)*(pthBins/(pthscale)));
				     
				    cout << "zT=" << zT << endl;
				    
				    // azimuthal angle of trigger
				    phi = fullEvent[p_i].phi();                     // azimuthal angle
				    // azimuthal angle of associated
				    phi2 = fullEvent[p_j].phi();                    // azimuthal angle
				    
				    deltaPhi = phi2-phi;
				    
				    deltaPhi = phi2-phi;
				    if (deltaPhi>PI)
				      deltaPhi = deltaPhi-2.*PI;
				    if (deltaPhi<-PI)
				      deltaPhi = deltaPhi+2.*PI;
				    
				    if(posZT>=0 && posZT<=ztBins && abs(deltaPhi-PI)<=0.63) 
				      {
					binning_corr_zT2[posZT] += 1;
				      }
				
				    if(pospTh>=0 && pospTh<=pthBins && abs(deltaPhi-PI)<=0.63) 
				      {
					binning_corr_pTh2[pospTh] += 1;
				      }
				    
				  }
			      }		   
			}
		    }
		}

	      // count pi_0s (111) or pi+ (211)
	      if ( id == 111 ) // pions (pi0)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  phi = asin(sqrt(pow(fullEvent[p_i].py(),2.))/p);        // azimuthal angle
		  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
		  //cout << "p_T=" << p << ", px=" << pythia.event[p_i].px() << ", py=" << pythia.event[p_i].py() 
		  //     << ", p_T=" << sqrt(pow(pythia.event[p_i].px(),2.)+pow(pythia.event[p_i].py(),2.)) 
		  //     << ", phi=" << phi << ", phi_deg=" << phi/PI*180. 
		  //     << ", posphi=" << posphi << endl;
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(eta)<=ymaxCorr3)
		    {
		      binning[posy] += 1;
		      pisum += 1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_pt[posphi][posy] += 1;
		    }
		  if(posy>=0 && posy<bins && abs(eta)>1.2 && abs(eta)<2.2)
		    {
		      binning_forward[posy] += 1;
		      if(posphi>=0 && posphi<phiBins)
			binning_phi_forward_pt[posphi][posy] += 1;
		    }
		}
	      if ( id == 2212 ) // p (p+)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_p[posy] += 1;
		      pSum += 1;
		    }
		}
	      if ( id == -2212 ) // pbar (p-)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_pbar[posy] += 1;
		      pbarSum += 1;
		    }
		}
	       if ( id == 211 ) // pions (pi+)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
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
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_piminus[posy] += 1;
		      piMinusSum += 1;
		    }
		}
	      if ( id == 113 ) // rho0
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_rho[posy] += 1;
		      rhoSum += 1;
		    }
		}
	     if ( id == 223 ) // omega (782)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_omega[posy] += 1;
		      omegaSum += 1;
		    }
		}
	     if ( id == 221 ) // eta
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_eta[posy] += 1;
		      etaSum += 1;
		    }
		}
	     if ( id == 333 ) // phi (1020)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_phiMeson[posy] += 1;
		      phiSum += 1;
		    }
		}
	     if ( id == 11 ) // e-
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_e[posy] += 1;
		      eSum += 1;
		    }
		}
	     if ( id == -11 ) // e+
	       {		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_eplus[posy] += 1;
		      eplusSum += 1;
		    }
		}
	     if ( id == 13 ) // mu-
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_mu[posy] += 1;
		      muSum += 1;
		    }
		}
	     if ( id == -13 ) // mu+
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_muplus[posy] += 1;
		      muplusSum += 1;
		    }
		}
	     if ( id == 443 ) // J/Psi
	       {		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_JPsi[posy] += 1;
		      JPsiSum += 1;
		    }
		}
	      if ( id == 321 ) // Kaons (K+)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
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
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_Kminus[posy] += 1;
		      KMinusSum += 1;
		    }
		}
	      if ( id == 311 ) // Kaons (K0)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_K0[posy] += 1;
		      K0Sum += 1;
		    }
		}
	      if ( id == 310 ) // Kaons (K0S)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_K0S[posy] += 1;
		      K0SSum += 1;
		    }
		}
	      if ( id == 130 ) // Kaons (K0L)
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_K0L[posy] += 1;
		      K0LSum += 1;
		    }
		}
	      if ( id == 22 && ( fullEvent[p_i].status()<90 || fullEvent[p_i].status()>99 ) ) // photons
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_gamma[posy] += 1;
		    }
		}
	      if ( id < 0 && fullEvent[p_i].isHadron()&& fullEvent[p_i].isFinal()
		   && fullEvent[p_i].charge()<0 ) // hadrons-
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_hm[posy] += 1;
		      hmsum += 1;
		    }
		}
	      if ( id > 0 && fullEvent[p_i].isHadron() && fullEvent[p_i].isFinal()
		   && fullEvent[p_i].charge()>0 ) // hadrons+
		{		
		  p = fullEvent[p_i].pT();                                // p_trans
		  posy = floor(p*(bins/scale));                              // bin number
		  pl = fullEvent[p_i].pz();                               // p_long
		  En = fullEvent[p_i].e();                                // energy
		  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
		  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  if(posy>=0 && posy<bins && abs(y)<=ymax)
		    {
		      binning_hp[posy] += 1;
		      hpsum += 1;
		    }
		}
	    }
	}

      for(int iy=0; iy<bins; iy++)
	{
	  binning_tmp[iy]=binning[iy]-binning_tmp[iy];
       	  binning_sq[iy]+=binning_tmp[iy]*binning_tmp[iy];
	  binning_tmp[iy]=binning[iy];

	  binning_forward_tmp[iy]=binning_forward[iy]-binning_forward_tmp[iy];
       	  binning_forward_sq[iy]+=binning_forward_tmp[iy]*binning_forward_tmp[iy];
	  binning_forward_tmp[iy]=binning_forward[iy];

	  binning_mu_tmp[iy]=binning_mu[iy]-binning_mu_tmp[iy];
       	  binning_mu_sq[iy]+=binning_mu_tmp[iy]*binning_mu_tmp[iy];
	  binning_mu_tmp[iy]=binning_mu[iy];

	  binning_muplus_tmp[iy]=binning_mu[iy]-binning_mu_tmp[iy];
       	  binning_muplus_sq[iy]+=binning_mu_tmp[iy]*binning_mu_tmp[iy];
	  binning_muplus_tmp[iy]=binning_mu[iy];

	  binning_JPsi_tmp[iy]=binning_JPsi[iy]-binning_JPsi_tmp[iy];
       	  binning_JPsi_sq[iy]+=binning_JPsi_tmp[iy]*binning_JPsi_tmp[iy];
	  binning_JPsi_tmp[iy]=binning_JPsi[iy];

	  binning_e_tmp[iy]=binning_e[iy]-binning_e_tmp[iy];
       	  binning_e_sq[iy]+=binning_e_tmp[iy]*binning_e_tmp[iy];
	  binning_e_tmp[iy]=binning_e[iy];

	  binning_eplus_tmp[iy]=binning_eplus[iy]-binning_eplus_tmp[iy];
       	  binning_eplus_sq[iy]+=binning_eplus_tmp[iy]*binning_eplus_tmp[iy];
	  binning_eplus_tmp[iy]=binning_eplus[iy];

	  binning_phiMeson_tmp[iy]=binning_phiMeson[iy]-binning_phiMeson_tmp[iy];
       	  binning_phiMeson_sq[iy]+=binning_phiMeson_tmp[iy]*binning_phiMeson_tmp[iy];
	  binning_phiMeson_tmp[iy]=binning_phiMeson[iy];

	  binning_eta_tmp[iy]=binning_eta[iy]-binning_eta_tmp[iy];
       	  binning_eta_sq[iy]+=binning_eta_tmp[iy]*binning_eta_tmp[iy];
	  binning_eta_tmp[iy]=binning_eta[iy];

	  binning_omega_tmp[iy]=binning_omega[iy]-binning_omega_tmp[iy];
       	  binning_omega_sq[iy]+=binning_omega_tmp[iy]*binning_omega_tmp[iy];
	  binning_omega_tmp[iy]=binning_omega[iy];

	  binning_rho_tmp[iy]=binning_rho[iy]-binning_rho_tmp[iy];
       	  binning_rho_sq[iy]+=binning_rho_tmp[iy]*binning_rho_tmp[iy];
	  binning_rho_tmp[iy]=binning_rho[iy];

	  binning_p_tmp[iy]=binning_p[iy]-binning_p_tmp[iy];
       	  binning_p_sq[iy]+=binning_p_tmp[iy]*binning_p_tmp[iy];
	  binning_p_tmp[iy]=binning_p[iy];

	  binning_pbar_tmp[iy]=binning_pbar[iy]-binning_pbar_tmp[iy];
       	  binning_pbar_sq[iy]+=binning_pbar_tmp[iy]*binning_pbar_tmp[iy];
	  binning_pbar_tmp[iy]=binning_pbar[iy];
	 
	  binning_piplus_tmp[iy]=binning_piplus[iy]-binning_piplus_tmp[iy];
       	  binning_piplus_sq[iy]+=binning_piplus_tmp[iy]*binning_piplus_tmp[iy];
	  binning_piplus_tmp[iy]=binning_piplus[iy];
	  
	  binning_piminus_tmp[iy]=binning_piminus[iy]-binning_piminus_tmp[iy];
       	  binning_piminus_sq[iy]+=binning_piminus_tmp[iy]*binning_piminus_tmp[iy];
	  binning_piminus_tmp[iy]=binning_piminus[iy];
	  
	  binning_Kminus_tmp[iy]=binning_Kminus[iy]-binning_Kminus_tmp[iy];
       	  binning_Kminus_sq[iy]+=binning_Kminus_tmp[iy]*binning_Kminus_tmp[iy];
	  binning_Kminus_tmp[iy]=binning_Kminus[iy];
	  
	  binning_Kplus_tmp[iy]=binning_Kplus[iy]-binning_Kplus_tmp[iy];
       	  binning_Kplus_sq[iy]+=binning_Kplus_tmp[iy]*binning_Kplus_tmp[iy];
	  binning_Kplus_tmp[iy]=binning_Kplus[iy];
	  
	  binning_K0_tmp[iy]=binning_K0[iy]-binning_K0_tmp[iy];
       	  binning_K0_sq[iy]+=binning_K0_tmp[iy]*binning_K0_tmp[iy];
	  binning_K0_tmp[iy]=binning_K0[iy];
	  
	  binning_K0S_tmp[iy]=binning_K0S[iy]-binning_K0S_tmp[iy];
       	  binning_K0S_sq[iy]+=binning_K0S_tmp[iy]*binning_K0S_tmp[iy];
	  binning_K0S_tmp[iy]=binning_K0S[iy];
	  
	  binning_K0L_tmp[iy]=binning_K0L[iy]-binning_K0L_tmp[iy];
       	  binning_K0L_sq[iy]+=binning_K0L_tmp[iy]*binning_K0L_tmp[iy];
	  binning_K0L_tmp[iy]=binning_K0L[iy];
	  
	  binning_hm_tmp[iy]=binning_hm[iy]-binning_hm_tmp[iy];
       	  binning_hm_sq[iy]+=binning_hm_tmp[iy]*binning_hm_tmp[iy];
	  binning_hm_tmp[iy]=binning_hm[iy];
	  
	  binning_hp_tmp[iy]=binning_hp[iy]-binning_hp_tmp[iy];
       	  binning_hp_sq[iy]+=binning_hp_tmp[iy]*binning_hp_tmp[iy];
	  binning_hp_tmp[iy]=binning_hp[iy];
	  
	  binning_gamma_tmp[iy]=binning_gamma[iy]-binning_gamma_tmp[iy];
       	  binning_gamma_sq[iy]+=binning_gamma_tmp[iy]*binning_gamma_tmp[iy];
	  binning_gamma_tmp[iy]=binning_gamma[iy];
	  for(int ix=0; ix<phiBins; ix++)
	    {
	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy]-binning_phi_tmp[ix][iy];
	      binning_phi_sq[ix][iy]+=binning_phi_tmp[ix][iy]*binning_phi_tmp[ix][iy];
	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy];

	      binning_phi_forward_tmp[ix][iy]=binning_phi_forward_pt[ix][iy]-binning_phi_forward_tmp[ix][iy];
	      binning_phi_forward_sq[ix][iy]+=binning_phi_forward_tmp[ix][iy]*binning_phi_forward_tmp[ix][iy];
	      binning_phi_forward_tmp[ix][iy]=binning_phi_forward_pt[ix][iy];
	    }
	}

      if(j%1000==0)
	{
	  cout << "#" << j  << endl; // write every 1000th time step
	  foutt << "#" << j << endl;
	}
    } // end loop over all events
  foutt.close();

  for(int iy=0; iy<corrPhiBins; iy++)
    {
      binning_corr_phi_tmp[iy]=binning_corr_phi[iy]-binning_corr_phi_tmp[iy];
      binning_corr_phi_sq[iy]+=binning_corr_phi_tmp[iy]*binning_corr_phi_tmp[iy];
      binning_corr_phi_tmp[iy]=binning_corr_phi[iy];
      binning_corr_phi2_tmp[iy]=binning_corr_phi2[iy]-binning_corr_phi2_tmp[iy];
      binning_corr_phi2_sq[iy]+=binning_corr_phi2_tmp[iy]*binning_corr_phi2_tmp[iy];
      binning_corr_phi2_tmp[iy]=binning_corr_phi2[iy];
      binning_corr_phi3_tmp[iy]=binning_corr_phi3[iy]-binning_corr_phi3_tmp[iy];
      binning_corr_phi3_sq[iy]+=binning_corr_phi3_tmp[iy]*binning_corr_phi3_tmp[iy];
      binning_corr_phi3_tmp[iy]=binning_corr_phi3[iy];
      binning_corr_phi4_tmp[iy]=binning_corr_phi4[iy]-binning_corr_phi4_tmp[iy];
      binning_corr_phi4_sq[iy]+=binning_corr_phi4_tmp[iy]*binning_corr_phi4_tmp[iy];
      binning_corr_phi4_tmp[iy]=binning_corr_phi4[iy];
    }

  for(int iy=0; iy<ztBins; iy++)
    {
      binning_corr_zT_tmp[iy]=binning_corr_zT[iy]-binning_corr_zT_tmp[iy];
      binning_corr_zT_sq[iy]+=binning_corr_zT_tmp[iy]*binning_corr_zT_tmp[iy];
      binning_corr_zT_tmp[iy]=binning_corr_zT[iy];
      binning_corr_zT2_tmp[iy]=binning_corr_zT2[iy]-binning_corr_zT2_tmp[iy];
      binning_corr_zT2_sq[iy]+=binning_corr_zT2_tmp[iy]*binning_corr_zT2_tmp[iy];
      binning_corr_zT2_tmp[iy]=binning_corr_zT2[iy];
      binning_corr_zT3_tmp[iy]=binning_corr_zT3[iy]-binning_corr_zT3_tmp[iy];
      binning_corr_zT3_sq[iy]+=binning_corr_zT3_tmp[iy]*binning_corr_zT3_tmp[iy];
      binning_corr_zT3_tmp[iy]=binning_corr_zT3[iy];
    }

  for(int iy=0; iy<pthBins; iy++)
    {
      binning_corr_pTh_tmp[iy]=binning_corr_pTh[iy]-binning_corr_pTh_tmp[iy];
      binning_corr_pTh_sq[iy]+=binning_corr_pTh_tmp[iy]*binning_corr_pTh_tmp[iy];
      binning_corr_pTh_tmp[iy]=binning_corr_pTh[iy];
      binning_corr_pTh2_tmp[iy]=binning_corr_pTh2[iy]-binning_corr_pTh2_tmp[iy];
      binning_corr_pTh2_sq[iy]+=binning_corr_pTh2_tmp[iy]*binning_corr_pTh2_tmp[iy];
      binning_corr_pTh2_tmp[iy]=binning_corr_pTh2[iy];
      binning_corr_pTh3_tmp[iy]=binning_corr_pTh3[iy]-binning_corr_pTh3_tmp[iy];
      binning_corr_pTh3_sq[iy]+=binning_corr_pTh3_tmp[iy]*binning_corr_pTh3_tmp[iy];
      binning_corr_pTh3_tmp[iy]=binning_corr_pTh3[iy];
    }
      
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

  fstream foutec("elasticCollisions.dat",ios::out); 
  foutec.precision(12);  
  foutec << "> Number_of_events = " << runs << endl; 
  foutec << "> quarks+gluons " << endl;
  for(int iy=0; iy<ecbins; iy++)
    {
      foutec << (iy)/(ecbins/ecscale)/2.+(iy+1)/(ecbins/ecscale)/2. << " " << ecbins*static_cast<double>(binning_ec[iy])/ecscale << endl; //devide by runs later
    }
  foutec.close();

  fstream foutx("xy.dat",ios::out); 
  foutx.precision(12);  
  foutx << "> Number_of_events = " << runs << endl; 
  foutx << "> y,x " << endl;
  for(int iy=0; iy<xbins; iy++)
    for(int ix=0; ix<xbins; ix++)
      {
	foutx << (iy)/(xbins/xscale) - 10 << " " << (ix)/(xbins/xscale) -10 << " " << xbins*xbins*binning_x[iy][ix]/xscale/xscale
	      << " " << xbins*xbins*static_cast<double>(binning_mom[iy][ix])/xscale/xscale
	      << endl; //devide by runs later
	if (ix == xbins-1) foutx << endl;
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

  fstream foutEini1("Eini1.dat",ios::out); 
  foutEini1.precision(12);  
  foutEini1 << "> Number_of_events = " << runs << endl; 
  foutEini1 << "> E_ini " << endl;
  for(int iy=0; iy<bins; iy++)
    {
      foutEini1 << (iy)/(bins/rscale) << " " << bins*static_cast<double>(binning_Eini1[iy])/rscale  << endl; //divide by runs later
    }
  foutEini1.close();

  fstream foutEini("EiniOfR.dat",ios::out); 
  foutEini.precision(12);  
  foutEini << "> Number_of_events = " << runs << endl; 
  foutEini << "> r, E_ini " << endl;
  for(int iy=0; iy<bins; iy++)
    {
      foutEini << (iy)/(bins/rscale) << " " << binning_Eini[iy]/static_cast<double>(binning_r[iy])  << endl; //divide by runs later
    }
  foutEini.close();

  fstream foutqg("initQandG.dat",ios::out); 
  foutqg.precision(12);  
  foutqg << "> Number_of_events = " << runs << endl; 
  foutqg << "> quarks, Gluons " << endl;
      foutqg << initQuarks << " " << initGluons << endl;
  foutqg.close();

  fstream foutL("distances.dat",ios::out); 
  foutL.precision(12);  
  foutL << "> Number_of_events = " << runs << endl; 
  foutL << "> distance traveled: quarks, gluons " << endl;
  for(int iy=0; iy<Lbins; iy++)
    {
      foutL << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
	   << Lbins*static_cast<double>(binning_distance_q[iy])/Lscale/(2.*ymax)/totalNumber_q << " " 
	    << Lbins*Lbins*static_cast<double>(binning_distance_q_sq[iy])/Lscale/(2.*ymax)/Lscale/(2.*ymax)/totalNumber_q/totalNumber_q
	    << " " << Lbins*static_cast<double>(binning_distance_g[iy])/Lscale/(2.*ymax)/totalNumber_g << " " 
	    << Lbins*Lbins*static_cast<double>(binning_distance_g_sq[iy])/Lscale/(2.*ymax)/Lscale/(2.*ymax)/totalNumber_g/totalNumber_g
	    << endl;
    }
  foutL.close();

  fstream foutLF("distancesForward.dat",ios::out); 
  foutLF.precision(12);  
  foutLF << "> Number_of_events = " << runs << endl; 
  foutLF << "> distance traveled: quarks, gluons, 1.2<|eta|<2.2 " << endl;
  for(int iy=0; iy<Lbins; iy++)
    {
      foutLF << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
	   << Lbins*static_cast<double>(binning_distance_q_forward[iy])/Lscale/(2.*ymax)/totalNumber_q_forward << " " 
	    << Lbins*Lbins*static_cast<double>(binning_distance_q_sq_forward[iy])/Lscale/(2.*ymax)/Lscale/(2.*ymax)/totalNumber_q_forward/totalNumber_q_forward
	    << " " << Lbins*static_cast<double>(binning_distance_g_forward[iy])/Lscale/(2.*ymax)/totalNumber_g_forward << " " 
	    << Lbins*Lbins*static_cast<double>(binning_distance_g_sq_forward[iy])/Lscale/(2.*ymax)/Lscale/(2.*ymax)/totalNumber_g_forward/totalNumber_g_forward
	    << endl;
    }
  foutLF.close();

  if(martini.returnTrackHistory())
    {
      fstream foutH("averageHistory.dat",ios::out); 
      foutH.precision(12);  
      foutH << "> Number_of_events = " << runs << endl; 
      foutH << "> t x y z dE/dt dpx/dt dpy/dt dpz/dt E px py QGPfrac counts x-error y-error" << endl;
      for(int iy=0; iy<static_cast<int>(20./dtfm); iy++)
	{
	  if(binning_history[7][iy]==0)
	    binning_history[7][iy]=1;
	  
	  foutH << ((iy)*(dtfm)+(iy+1)*dtfm)/2. << " " 
		<< binning_history[0][iy]/binning_history[7][iy] << " " 
		<< binning_history[1][iy]/binning_history[7][iy] 
		<< " " << binning_history[2][iy]/binning_history[7][iy] 
		<< " " << -binning_history[3][iy]/binning_history[7][iy] 
		<< " " << -binning_history[4][iy]/binning_history[7][iy] 
		<< " " << -binning_history[5][iy]/binning_history[7][iy] 
		<< " " << -binning_history[6][iy]/binning_history[7][iy]
		<< " " << binning_history[8][iy]/binning_history[7][iy]	
		<< " " << binning_history[9][iy]/binning_history[7][iy]	
		<< " " << binning_history[10][iy]/binning_history[7][iy]	
		<< " " << binning_history[11][iy]/binning_history[7][iy]	
		<< " " << binning_history[7][iy] 
		<< " " 
		<< (binning_history[14][iy]/binning_history[7][iy]-pow(binning_history[0][iy]/binning_history[7][iy],2.))/sqrt(binning_history[7][iy]) 
		<< " " 
		<< (binning_history[15][iy]/binning_history[7][iy]-pow(binning_history[1][iy]/binning_history[7][iy],2.))/sqrt(binning_history[7][iy]) 
		<< endl;
	}
      foutH.close();
    }
  
// 		  binning_history[13][post]+=martini.Tt->at(ni);
// 		  binning_history[0][post]+=martini.Tx->at(ni);
// 		  binning_history[1][post]+=martini.Ty->at(ni);
// 		  binning_history[2][post]+=martini.Tz->at(ni);
// 		  binning_history[3][post]+=martini.TdEdt->at(ni);
// 		  binning_history[4][post]+=martini.Tdpxdt->at(ni);
// 		  binning_history[5][post]+=martini.Tdpydt->at(ni);
// 		  binning_history[6][post]+=martini.Tdpzdt->at(ni);
// 		  binning_history[7][post]+=1;// counts the entries
// 		  binning_history[8][post]+=martini.TE->at(ni);
// 		  binning_history[9][post]+=martini.Tpx->at(ni);
// 		  binning_history[10][post]+=martini.Tpy->at(ni);
//		  binning_history[11][post]+=martini.TQGPfrac->at(ni);
  
  fstream foutphiLq("distancesQuarkPhi.dat",ios::out); 
  foutphiLq << "> quarks. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
  foutphiLq.precision(12);  
  foutphiLq << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<Lbins; iy++)
    {
      foutphiLq << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[0][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[0] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[0][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[0]/totalNumber_phi_q[0]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[1][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[1] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[1][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[1]/totalNumber_phi_q[1]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[2][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[2] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[2][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[2]/totalNumber_phi_q[2]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[3][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[3] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[3][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[3]/totalNumber_phi_q[3]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[4][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[4] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[4][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[4]/totalNumber_phi_q[4]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q[5][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q[5] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q[5][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q[5]/totalNumber_phi_q[5]
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphiLq.close();

  fstream foutphiLqF("distancesQuarkPhiForward.dat",ios::out); 
  foutphiLqF << "> quarks. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90, 1.2<|eta|<2.2" << endl;
  foutphiLqF.precision(12);  
  foutphiLqF << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<Lbins; iy++)
    {
      foutphiLqF << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[0][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[0] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[0][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[0]/totalNumber_phi_q_forward[0]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[1][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[1] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[1][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[1]/totalNumber_phi_q_forward[1]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[2][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[2] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[2][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[2]/totalNumber_phi_q_forward[2]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[3][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[3] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[3][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[3]/totalNumber_phi_q_forward[3]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[4][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[4] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[4][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[4]/totalNumber_phi_q_forward[4]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_q_forward[5][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_q_forward[5] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_q_forward[5][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_q_forward[5]/totalNumber_phi_q_forward[5]
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphiLqF.close();


  fstream foutphiLg("distancesGluonPhi.dat",ios::out); 
  foutphiLg << "> gluons. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
  foutphiLg.precision(12);  
  foutphiLg << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<Lbins; iy++)
    {
      foutphiLg << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[0][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[0] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[0][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[0]/totalNumber_phi_g[0]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[1][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[1] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[1][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[1]/totalNumber_phi_g[1]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[2][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[2] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[2][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[2]/totalNumber_phi_g[2]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[3][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[3] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[3][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[3]/totalNumber_phi_g[3]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[4][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[4] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[4][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[4]/totalNumber_phi_g[4]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g[5][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g[5] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g[5][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g[5]/totalNumber_phi_g[5]
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphiLg.close();

  fstream foutphiLgF("distancesGluonPhiForward.dat",ios::out); 
  foutphiLgF << "> gluons. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90, 1.2<|eta|<2.2" << endl;
  foutphiLgF.precision(12);  
  foutphiLgF << "> Number_of_events = " << numEvents << endl; 
  for(int iy=0; iy<Lbins; iy++)
    {
      foutphiLgF << ((iy)/(Lbins/Lscale)+(iy+1)/(Lbins/Lscale))/2. << " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[0][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[0] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[0][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[0]/totalNumber_phi_g_forward[0]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[1][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[1] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[1][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[1]/totalNumber_phi_g_forward[1]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[2][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[2] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[2][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[2]/totalNumber_phi_g_forward[2]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[3][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[3] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[3][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[3]/totalNumber_phi_g_forward[3]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[4][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[4] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[4][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[4]/totalNumber_phi_g_forward[4]<< " " 
		<< Lbins*static_cast<double>(binning_phi_distance_g_forward[5][iy])/Lscale/(2.*ymaxPartons)/totalNumber_phi_g_forward[5] << " " 
		<< Lbins*Lbins*static_cast<double>(binning_phi_distance_sq_g_forward[5][iy])/Lscale/(2.*ymaxPartons)/Lscale/(2.*ymaxPartons)
	/totalNumber_phi_g_forward[5]/totalNumber_phi_g_forward[5]
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphiLgF.close();


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
		<< bins*static_cast<double>(binning_phi_pt_q[0][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[0][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[1][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[1][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[2][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[2][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[3][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[3][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[4][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[4][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_q[5][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_q[5][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)
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
		<< bins*static_cast<double>(binning_phi_pt_g[0][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[0][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[1][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[1][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[2][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[2][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[3][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[3][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[4][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[4][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)<< " " 
		<< bins*static_cast<double>(binning_phi_pt_g[5][iy])/scale/(2.*ymaxPartons) << " " 
		<< bins*bins*static_cast<double>(binning_phi_sq_g[5][iy])/scale/(2.*ymaxPartons)/scale/(2.*ymaxPartons)
		<< endl; //devide by numEvents later 
      // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
      // second output is the sum of squares to calculate the error in the end. 
    }
  foutphig.close();

  cout << endl;
  if (martini.returnFragmentationSwitch() == 1)
    {
      fstream foutcp("corrPhi.dat",ios::out); 
      foutcp << "> deltaPhi:" << endl;
      foutcp.precision(12);  
      foutcp << "> Number_of_events = " << numEvents << endl; 
      foutcp << "> Number_of_triggers = " << Ntrigger << endl; 
      for(int iy=0; iy<corrPhiBins; iy++)
	{
	  foutcp << ((iy)/(corrPhiBins/(2.*PI))+(iy+1)/(corrPhiBins/(2.*PI)))/2.-PI << " " 
		 << corrPhiBins*static_cast<double>(binning_corr_phi[iy])/(2.*PI) << " " 
		 << corrPhiBins*corrPhiBins*static_cast<double>(binning_corr_phi_sq[iy])/(2.*PI)/(2.*PI) << " "
		 << static_cast<double>(binning_corr_phi[iy]) << " " 
		 << endl; //devide by numEvents later 
	  // (scale/bins)=\Delta phi = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutcp.close();

      fstream foutcp2("corrPhi2.dat",ios::out); 
      foutcp2 << "> deltaPhi:" << endl;
      foutcp2.precision(12);  
      foutcp2 << "> Number_of_events = " << numEvents << endl; 
      foutcp2 << "> Number_of_triggers = " << Ntrigger2 << endl; 
      for(int iy=0; iy<corrPhiBins; iy++)
	{
	  foutcp2 << ((iy)/(corrPhiBins/(2.*PI))+(iy+1)/(corrPhiBins/(2.*PI)))/2.-PI << " " 
		 << corrPhiBins*static_cast<double>(binning_corr_phi2[iy])/(2.*PI) << " " 
		 << corrPhiBins*corrPhiBins*static_cast<double>(binning_corr_phi2_sq[iy])/(2.*PI)/(2.*PI) << " "
		 << static_cast<double>(binning_corr_phi2[iy]) << " " 
		 << endl; //devide by numEvents later 
	  // (scale/bins)=\Delta phi = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutcp2.close();

      fstream foutcp3("corrPhi3.dat",ios::out); 
      foutcp3 << "> deltaPhi:" << endl;
      foutcp3.precision(12);  
      foutcp3 << "> Number_of_events = " << numEvents << endl; 
      foutcp3 << "> Number_of_triggers = " << Ntrigger3 << endl; 
      for(int iy=0; iy<corrPhiBins; iy++)
	{
	  foutcp3 << ((iy)/(corrPhiBins/(2.*PI))+(iy+1)/(corrPhiBins/(2.*PI)))/2.-PI << " " 
		 << corrPhiBins*static_cast<double>(binning_corr_phi3[iy])/(2.*PI) << " " 
		 << corrPhiBins*corrPhiBins*static_cast<double>(binning_corr_phi3_sq[iy])/(2.*PI)/(2.*PI) << " "
		 << static_cast<double>(binning_corr_phi3[iy]) << " " 
		 << endl; //devide by numEvents later 
	  // (scale/bins)=\Delta phi = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutcp3.close();

      fstream foutcp4("corrPhi4.dat",ios::out); 
      foutcp4 << "> deltaPhi:" << endl;
      foutcp4.precision(12);  
      foutcp4 << "> Number_of_events = " << numEvents << endl; 
      foutcp4 << "> Number_of_triggers = " << Ntrigger4 << endl; 
      for(int iy=0; iy<corrPhiBins; iy++)
	{
	  foutcp4 << ((iy)/(corrPhiBins/(2.*PI))+(iy+1)/(corrPhiBins/(2.*PI)))/2.-PI << " " 
		 << corrPhiBins*static_cast<double>(binning_corr_phi4[iy])/(2.*PI) << " " 
		 << corrPhiBins*corrPhiBins*static_cast<double>(binning_corr_phi4_sq[iy])/(2.*PI)/(2.*PI) << " "
		 << static_cast<double>(binning_corr_phi4[iy]) << " " 
		 << endl; //devide by numEvents later 
	  // (scale/bins)=\Delta phi = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutcp4.close();

      fstream foutzt("zT.dat",ios::out); 
      foutzt << "> zT:" << endl;
      foutzt.precision(12);  
      foutzt << "> Number_of_events = " << numEvents << endl; 
      foutzt << "> Number_of_triggers = " << NPhotontrigger << endl; 
      for(int iy=0; iy<ztBins; iy++)
	{
	  foutzt << ((iy)/(ztBins/(ztscale))+(iy+1)/(ztBins/(ztscale)))/2. << " " 
		 << ztBins*static_cast<double>(binning_corr_zT[iy])/(ztscale) << " " 
		 << ztBins*ztBins*static_cast<double>(binning_corr_zT_sq[iy])/(ztscale)/(ztscale) << " "
		 << static_cast<double>(binning_corr_zT[iy]) << " " 
		 << endl; //devide by numEvents later 
	}
      foutzt.close();

      fstream foutpth("pTh.dat",ios::out); 
      foutpth << "> zT:" << endl;
      foutpth.precision(12);  
      foutpth << "> Number_of_events = " << numEvents << endl; 
      foutpth << "> Number_of_triggers = " << NPhotontrigger << endl; 
      for(int iy=0; iy<pthBins; iy++)
	{
	  foutpth << ((iy)/(pthBins/(pthscale))+(iy+1)/(pthBins/(pthscale)))/2. << " " 
		  << pthBins*static_cast<double>(binning_corr_pTh[iy])/(pthscale) << " " 
		  << pthBins*pthBins*static_cast<double>(binning_corr_pTh_sq[iy])/(pthscale)/(pthscale) << " "
		  << static_cast<double>(binning_corr_pTh[iy]) << " " 
		  << endl; //devide by numEvents later 
	}
      foutpth.close();

      fstream foutzt2("zT2.dat",ios::out); 
      foutzt2 << "> zT:" << endl;
      foutzt2.precision(12);  
      foutzt2 << "> Number_of_events = " << numEvents << endl; 
      foutzt2 << "> Number_of_triggers = " << NPhotontrigger2 << endl; 
      for(int iy=0; iy<ztBins; iy++)
	{
	  foutzt2 << ((iy)/(ztBins/(ztscale))+(iy+1)/(ztBins/(ztscale)))/2. << " " 
		 << ztBins*static_cast<double>(binning_corr_zT2[iy])/(ztscale) << " " 
		 << ztBins*ztBins*static_cast<double>(binning_corr_zT2_sq[iy])/(ztscale)/(ztscale) << " "
		 << static_cast<double>(binning_corr_zT2[iy]) << " " 
		 << endl; //devide by numEvents later 
	}
      foutzt2.close();

      fstream foutpth2("pTh2.dat",ios::out); 
      foutpth2 << "> zT:" << endl;
      foutpth2.precision(12);  
      foutpth2 << "> Number_of_events = " << numEvents << endl; 
      foutpth2 << "> Number_of_triggers = " << NPhotontrigger2 << endl; 
      for(int iy=0; iy<pthBins; iy++)
	{
	  foutpth2 << ((iy)/(pthBins/(pthscale))+(iy+1)/(pthBins/(pthscale)))/2. << " " 
		  << pthBins*static_cast<double>(binning_corr_pTh2[iy])/(pthscale) << " " 
		  << pthBins*pthBins*static_cast<double>(binning_corr_pTh2_sq[iy])/(pthscale)/(pthscale) << " "
		  << static_cast<double>(binning_corr_pTh2[iy]) << " " 
		  << endl; //devide by numEvents later 
	}
      foutpth2.close();

      fstream foutzt3("zT3.dat",ios::out); 
      foutzt3 << "> zT:" << endl;
      foutzt3.precision(12);  
      foutzt3 << "> Number_of_events = " << numEvents << endl; 
      foutzt3 << "> Number_of_triggers = " << NPhotontrigger3 << endl; 
      for(int iy=0; iy<ztBins; iy++)
	{
	  foutzt3 << ((iy)/(ztBins/(ztscale))+(iy+1)/(ztBins/(ztscale)))/2. << " " 
		 << ztBins*static_cast<double>(binning_corr_zT3[iy])/(ztscale) << " " 
		 << ztBins*ztBins*static_cast<double>(binning_corr_zT3_sq[iy])/(ztscale)/(ztscale) << " "
		 << static_cast<double>(binning_corr_zT3[iy]) << " " 
		 << endl; //devide by numEvents later 
	}
      foutzt3.close();

      fstream foutpth3("pTh3.dat",ios::out); 
      foutpth3 << "> zT:" << endl;
      foutpth3.precision(12);  
      foutpth3 << "> Number_of_events = " << numEvents << endl; 
      foutpth3 << "> Number_of_triggers = " << NPhotontrigger3 << endl; 
      for(int iy=0; iy<pthBins; iy++)
	{
	  foutpth3 << ((iy)/(pthBins/(pthscale))+(iy+1)/(pthBins/(pthscale)))/2. << " " 
		  << pthBins*static_cast<double>(binning_corr_pTh3[iy])/(pthscale) << " " 
		  << pthBins*pthBins*static_cast<double>(binning_corr_pTh3_sq[iy])/(pthscale)/(pthscale) << " "
		  << static_cast<double>(binning_corr_pTh3[iy]) << " " 
		  << endl; //devide by numEvents later 
	}
      foutpth3.close();

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

      fstream foutF("pi0-forward.dat",ios::out); 
      foutF << "> pi0:" << endl;
      foutF.precision(12);  
      foutF << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutF << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_forward[iy])/scale/(2.*ymax) << " " 
	       << bins*bins*static_cast<double>(binning_forward_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
	       << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutF.close();

      fstream foutphi("pi0phi.dat",ios::out); 
      foutphi << "> pi0. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
      foutphi.precision(12);  
      foutphi << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutphi << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_phi_pt[0][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[0][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[1][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[1][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[2][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[2][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[3][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[3][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[4][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[4][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_pt[5][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_sq[5][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)
	       << endl; //devide by numEvents later 
	  // divide by 2*ymaxCorr3 to get 1 unit of rapidity (\Delta y = 2*ymaxCorr3). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutphi.close();
      
      fstream foutphiF("pi0phi-forward.dat",ios::out); 
      foutphiF << "> pi0. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
      foutphiF.precision(12);  
      foutphiF << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutphiF << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[0][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[0][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[1][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[1][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[2][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[2][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[3][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[3][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[4][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[4][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)<< " " 
	       << bins*static_cast<double>(binning_phi_forward_pt[5][iy])/scale/(2.*ymaxCorr3) << " " 
	       << bins*bins*static_cast<double>(binning_phi_forward_sq[5][iy])/scale/(2.*ymaxCorr3)/scale/(2.*ymaxCorr3)
	       << endl; //devide by numEvents later 
	  // divide by 2*ymaxCorr3 to get 1 unit of rapidity (\Delta y = 2*ymaxCorr3). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutphiF.close();
      
      fstream fout2("hminus.dat",ios::out); 
      fout2 << endl << "> h-:" << endl;
      fout2.precision(12);  
      fout2 << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fout2 << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		<< bins*static_cast<double>(binning_hm[iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_hm_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		<< endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	}
      fout2.close();

      fstream fouthp("hplus.dat",ios::out); 
      fouthp << endl << "> h+:" << endl;
      fouthp.precision(12);  
      fouthp << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fouthp << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		 << bins*static_cast<double>(binning_hp[iy])/scale/(2.*ymax) << " " 
		 << bins*bins*static_cast<double>(binning_hp_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		 << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	}
      fouthp.close();

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
		     << bins*static_cast<double>(binning_piplus[iy])/scale/(2.*ymax) << " " 
		     << bins*bins*static_cast<double>(binning_piplus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		     << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutpiplus.close();

      fstream foutrho("rho0.dat",ios::out); 
      foutrho << "> rho0:" << endl;
      foutrho.precision(12);  
      foutrho << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutrho << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_rho[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_rho_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutrho.close();

      fstream fouteta("eta.dat",ios::out); 
      fouteta << "> eta:" << endl;
      fouteta.precision(12);  
      fouteta << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fouteta << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_eta[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_eta_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      fouteta.close();

      fstream foutmu("mu.dat",ios::out); 
      foutmu << "> mu-:" << endl;
      foutmu.precision(12);  
      foutmu << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutmu << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		 << bins*static_cast<double>(binning_mu[iy])/scale/(2.*ymax) << " " 
		 << bins*bins*static_cast<double>(binning_mu_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		 << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutmu.close();

      fstream foutmuplus("muplus.dat",ios::out); 
      foutmuplus << "> mu+:" << endl;
      foutmuplus.precision(12);  
      foutmuplus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutmuplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		     << bins*static_cast<double>(binning_muplus[iy])/scale/(2.*ymax) << " " 
		     << bins*bins*static_cast<double>(binning_muplus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		     << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutmuplus.close();

      fstream foutphim("phi.dat",ios::out); 
      foutphim << "> phi(1020):" << endl;
      foutphim.precision(12);  
      foutphim << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutphim << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_phiMeson[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_phiMeson_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutphim.close();

      fstream foutJPsi("JPsi.dat",ios::out); 
      foutJPsi << "> JPsi(1S):" << endl;
      foutJPsi.precision(12);  
      foutJPsi << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutJPsi << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		   << bins*static_cast<double>(binning_JPsi[iy])/scale/(2.*ymax) << " " 
		   << bins*bins*static_cast<double>(binning_JPsi_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		   << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutJPsi.close();

      fstream foute("electrons.dat",ios::out); 
      foute << "> e-:" << endl;
      foute.precision(12);  
      foute << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foute << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_e[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_e_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foute.close();

      fstream fouteplus("positrons.dat",ios::out); 
      fouteplus << "> e+:" << endl;
      fouteplus.precision(12);  
      fouteplus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  fouteplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		    << bins*static_cast<double>(binning_eplus[iy])/scale/(2.*ymax) << " " 
		    << bins*bins*static_cast<double>(binning_eplus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		    << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      fouteplus.close();

      fstream foutomega("omega.dat",ios::out); 
      foutomega << "> omega:" << endl;
      foutomega.precision(12);  
      foutomega << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutomega << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		    << bins*static_cast<double>(binning_omega[iy])/scale/(2.*ymax) << " " 
		    << bins*bins*static_cast<double>(binning_omega_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		    << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutomega.close();

      fstream foutp("protons.dat",ios::out); 
      foutp << "> p:" << endl;
      foutp.precision(12);  
      foutp << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutp << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		<< bins*static_cast<double>(binning_p[iy])/scale/(2.*ymax) << " " 
		<< bins*bins*static_cast<double>(binning_p_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		<< endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutp.close();

      fstream foutpbar("anti-protons.dat",ios::out); 
      foutpbar << "> pbar:" << endl;
      foutpbar.precision(12);  
      foutpbar << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutpbar << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		   << bins*static_cast<double>(binning_pbar[iy])/scale/(2.*ymax) << " " 
		   << bins*bins*static_cast<double>(binning_pbar_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		   << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutpbar.close();

      fstream foutpiminus("pi-.dat",ios::out); 
      foutpiminus << "> pi-:" << endl;
      foutpiminus.precision(12);  
      foutpiminus << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutpiminus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		      << bins*static_cast<double>(binning_piminus[iy])/scale/(2.*ymax) << " " 
		      << bins*bins*static_cast<double>(binning_piminus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		      << endl; //devide by numEvents later 
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
		    << bins*static_cast<double>(binning_Kplus[iy])/scale/(2.*ymax) << " " 
		    << bins*bins*static_cast<double>(binning_Kplus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		    << endl; //devide by numEvents later 
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
		     << bins*static_cast<double>(binning_Kminus[iy])/scale/(2.*ymax) << " " 
		     << bins*bins*static_cast<double>(binning_Kminus_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		     << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutKminus.close();

      fstream foutK0("K0.dat",ios::out); 
      foutK0 << "> K0:" << endl;
      foutK0.precision(12);  
      foutK0 << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutK0 << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		     << bins*static_cast<double>(binning_K0[iy])/scale/(2.*ymax) << " " 
		     << bins*bins*static_cast<double>(binning_K0_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		     << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutK0.close();

      fstream foutK0S("K0S.dat",ios::out); 
      foutK0S << "> K0S:" << endl;
      foutK0S.precision(12);  
      foutK0S << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutK0S << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_K0S[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_K0S_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutK0S.close();

      fstream foutK0L("K0L.dat",ios::out); 
      foutK0L << "> K0L:" << endl;
      foutK0L.precision(12);  
      foutK0L << "> Number_of_events = " << numEvents << endl; 
      for(int iy=0; iy<bins; iy++)
	{
	  foutK0L << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
		  << bins*static_cast<double>(binning_K0L[iy])/scale/(2.*ymax) << " " 
		  << bins*bins*static_cast<double>(binning_K0L_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
		  << endl; //devide by numEvents later 
	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
	  // second output is the sum of squares to calculate the error in the end. 
	}
      foutK0L.close();
    }  

  cout.precision(6);

  qhatTotal/=static_cast<double>(partonsInQhat);
  dEdxTotal/=static_cast<double>(partonsInQhat);
  ptsqTotal/=static_cast<double>(partonsInQhat);
  dETotal/=static_cast<double>(partonsInQhat);
  
  if (martini.returnFragmentationSwitch() == 1)
    {
      cout << "pi+/pi-=" << static_cast<double>(piPlusSum)/static_cast<double>(piMinusSum) <<endl;
      cout << "K+/pi+=" << static_cast<double>(KPlusSum)/static_cast<double>(piPlusSum) <<endl;
      cout << "K-/pi-=" << static_cast<double>(KMinusSum)/static_cast<double>(piMinusSum) <<endl;
      cout << "K+/K-=" << static_cast<double>(KPlusSum)/static_cast<double>(KMinusSum) <<endl;
      cout << "K+=" << static_cast<double>(KPlusSum) << endl;
      cout << "K-=" << static_cast<double>(KMinusSum) << endl;
      cout << "pi+=" << static_cast<double>(piPlusSum) << endl;
      cout << "pi-=" << static_cast<double>(piMinusSum) << endl;
      cout << "pisum=" << pisum << endl;
      cout << "sum of negatively charged hadrons=" << hmsum << endl;
    }

  cout << "sigmaGen=" << martini.pythia.info.sigmaGen() << endl;
  cout << "used events=" << numEvents << endl;
  cout << "time steps=" << mt << endl;
  cout << "total number of particles=" << totalSum << endl; 
  cout << "total number of partons != u,d,s,g (mainly diquarks with p_t < 1 GeV):" << otherPartons << endl; 
  cout << "<q_T^2>/mfp =" << qtTotAll/maxTime << " GeV^2/fm" << endl;
  cout << "number of collisions =" << NColTot << endl;
  cout << "partons: " << p_imax << endl;
  cout << "mfp=" << mfpTot/runs  << " fm" << endl;
  cout << "average distance travelled=" << mfpTot/runs*NColTot  << " fm" << endl;
  cout << "average elastic collisions per parton=" 
       << static_cast<double>(totalElasticCollisions)/static_cast<double>(partonsInElasticCollisions) << endl;
  cout << "pt^2 =" << ptsqTotal << " GeV^2" << endl;
  cout << "delta E =" << dETotal << " GeV" << endl;
  cout << "qhat =" << qhatTotal << " GeV^2/fm" << endl;
  cout << "dE/dx =" << dEdxTotal << " GeV/fm" << endl;

  cout << "sqrt(<px^2>)=" << pxTot << endl;
  cout << "sqrt(<py^2>)=" << pyTot << endl;
  cout << "sqrt(<px^2>)_Trigger=" << pxTotTrig << endl;
  cout << "sqrt(<py^2>)_Trigger=" << pyTotTrig << endl;
  cout << "sqrt(Sum<px^2>)=" << pxAllTot << endl;
  cout << "sqrt(Sum<py^2>)=" << pyAllTot << endl;
  cout << "sqrt(Sum<px^2>)_Trigger=" << pxAllTotTrig << endl;
  cout << "sqrt(Sum<py^2>)_Trigger=" << pyAllTotTrig << endl;
  ofstream foutqhat("./output/qhatvstimeb7.5x.dat",ios::app); 
  foutqhat.precision(4);  
  foutqhat << maxTime-0.6 << " " << qhatTotal << " " << dEdxTotal << " " << ptsqTotal << " " << dETotal << endl; 
  foutqhat.close();

  delete plist[0];
  delete plist[1];
  delete plist;
  delete plistInitial;
  delete random;
}

