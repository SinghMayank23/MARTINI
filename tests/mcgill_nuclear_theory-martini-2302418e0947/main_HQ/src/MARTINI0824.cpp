// MARTINI.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class
// MARTINI: provide the main user interface to everything else.
// it also contains the main evolution routine

#define AMYpCut 4.01
#include "LorentzBooster.cxx"
#include "MARTINI.h"

// constructor
MARTINI::MARTINI()
{
  import  = new Import();
  importLRates  = new ImportLRates();
  random  = new Random();
  testing = new Testing();
  hydroSetup = new HydroSetup();
  glauber = new Glauber();
  settings.init();
  gsl_rand = gsl_rng_alloc(gsl_rng_taus);
}

// destructor
MARTINI::~MARTINI()
{
  delete import;
  delete importLRates;
  delete random;
  delete rates;
  delete testing;
  delete hydroSetup;
  delete glauber;
  delete gsl_rand;
}

// evolve every parton by one time step. This is the core of MARTINI.
int MARTINI::evolve(vector<Parton>  **plist, int counter, int it)
{
  //ofstream foutt("test.dat",ios::app); 
  ofstream fouttracks, fouttracksg, fouttracks2, fouttracksg2, foutpovray, foutpovrayg;
  if ( trackPartons )
    {
      fouttracks.open("./output/qtracks.dat",ios::app); 
      fouttracksg.open("./output/gtracks.dat",ios::app); 
      fouttracks2.open("./output/qtracks2.dat",ios::app); 
      fouttracksg2.open("./output/gtracks2.dat",ios::app); 
      foutpovray.open("./output/povraytracks.dat",ios::app); 
      foutpovrayg.open("./output/povraytracksg.dat",ios::app); 
    }
  cout.precision(6);

  HydroInfo hydroInfo;
  HydroInfo hydroInfoEarlier;
  ReturnValue f;                      // will contain \Delta p, the change in momentum due to a process
  Norms n;                            // holds the integrals \int dGamma_i dp, currently three for all radiative process
  double qqRate, qgRate;              // hold the integrals \int dGamma_i domega for the elastic processes
  double gqRate, ggRate;
  double conversionqg;                // hold the total conversion rate q->g
  double conversionqgamma;            // hold the total conversion rate q->gamma
  double conversiongq;                // hold the total conversion rate g->q

  Vec4 vecp, newvecp, vecet;          // momentum four-vectors
  Vec4 vecqt;                         // transverse momentum four-vector
  Vec4 vecpRest;
  Parton newOne;                      // parton object for the additionally created parton

  int imax = plist[0]->size();        // number of partons in the list
  int id, newID;                      // original and new parton's ID
  int col, acol, newCol, newAcol;     // number of the color string attached to the old, and new parton
  int ix, iy, iz, itau;               // parton's position in cell coordinates
  int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction
  int radiate;                        // switches to hold what process will happen
  int radiatePhoton;
  int radiate_ggg;
  int radiate_gqq;
  int elastic_qq;
  int elastic_qg;
  int elastic_gq;
  int elastic_gg;
  int convert_quark;
  int convert_quark_to_gamma;
  int convert_gluon;
  int GluonToQuark;                   // holds the decision whether a gluon converts to a quark (=1) or anti-quark (=0)

  double dt = dtfm/hbarc;             // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
                                      // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
  double T;                           // temperature
  double QGPfrac;                     // fraction of QGP in the mixed phase
  double p, q, qt;                    // parton momentum in [GeV]
  double x, y, z;                     // parton's position in [fm]
  double px, py, pz;                  // parton's momentum
  double t, tau;                      // lab time and lab tau
  double vx, vy, vz, veta;            // flow velocity of the cell
  double vetaZ, vetaZTau;             // flow velocity of the cell
  double beta;                        // absolute value of the flow velocity
  double gamma;                       // gamma factor 
  double cosPhi, cosPhiRest;          // angle between pvec and v, and angle between pvec_restframe and v
  double cosPhiRestEl;                // angle between pvec_restframe and v in case of elastic collisions with transv. mom. transfer
  double boostBack;                   // boost factor to boost back from rest- to lab-frame
  double pRest;                       // rest frame abs. value of momentum
  double pxRest, pyRest, pzRest;      // rest frame momentum components
  double newpx, newpy, newpz;         // temporary momentum components
  double omega;                       // energy transfered in elastic collision
  double tempAbs;                     // temporary storage for an absolute value of a momentum vector
  double totalQuarkProb;              // total probability that the quark undergoes some interaction
  double totalGluonProb;              // total probability that the gluon undergoes some interaction
  double totalDeltaE;                 // total energy loss in one time step 
  double totalDeltapx;                // total momentum px change in one time step  
  double totalDeltapy;       
  double totalDeltapz;
  const double rd=12.3;               // ratio of degrees of freedom g_QGP/g_H (for Kolb hydro only)
  
  ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);

  totalDeltaE = 0.;
  totalDeltapx = 0.;
  totalDeltapy = 0.;
  totalDeltapz = 0.;

  if ( hydroWhichHydro == 3 || hydroWhichHydro == 4 ) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);

  for ( int i=0; i<imax; i++)                                   // loop over all partons
    { 
      elastic_qq = 0;                                           // in the beginning assume nothing will happen
      elastic_qg = 0;                                           // then decide according to rates if processes occur 
      elastic_gq = 0;                                           // then decide according to rates if processes occur 
      elastic_gg = 0;
      radiate = 0;
      radiatePhoton = 0;
      radiate_ggg = 0;
      radiate_gqq = 0;
      convert_quark_to_gamma = 0;
      convert_quark = 0;
      convert_gluon = 0;


      counter++;                                                // counter that characterizes every step uniquely
                                                                // it is used to give new partons a unique color index
     
      if (plist[0]->at(i).frozen()==1 && !trackHistory) continue; 
                                                                // if a particle is frozen out, don't evolve it any further

      id = plist[0]->at(i).id();                                // id of parton i 
      vecp = plist[0]->at(i).p();                               // four-vector p of parton i
      p = vecp.pAbs();                                          // |p| of parton i 

      x = plist[0]->at(i).x();                                  // x value of position in [fm]
      y = plist[0]->at(i).y();                                  // y value of position in [fm]
      z = plist[0]->at(i).z();                                  // z value of position in [fm]

      // boost - using original position...
      ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
      iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
                                                                // note that x and y run from -hydroXmax to +hydroXmax
	                                                        // and ix and iy from 0 to 2*hydroXmax
	                                                        // hence the (hydroXmax+x or y) for both
      if (hydroWhichHydro == 3 || hydroWhichHydro == 4 ) 
	{	
	  iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
	}


      if (fixedTemperature == 0)
	{
	  t = it*dtfm;        //had a +tau0 here. now let partons travel doing their vacuum shower stuff first...
	  tau = 0.;           //just give tau some value - will change below

	  plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
	  // +++ if you want to move the partons before tau0 do the position update here.
	  if ( moveBeforeTau0 == 1 )
	    {
	      plist[0]->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
	      plist[0]->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
	      plist[0]->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
	    }

	  if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
	  if (tau<hydroTau0||tau>hydroTauMax-1) continue;       // do not evolve if tau<tau0 or tau>tauMax

	  // +++ if you DO NOT want to move the partons before tau0 do the position update here.
	  if ( moveBeforeTau0 == 0 )
	    {
	      plist[0]->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
	      plist[0]->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
	      plist[0]->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
	    }

	  // get the temperature and flow velocity at the current position and time:
	  hydroInfo = hydroSetup->getHydroValues(x, y, z, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
						 hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
	
	  T = hydroInfo.T;

	  if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 1 || hydroWhichHydro == 5) 
	    {
	      QGPfrac = hydroInfo.QGPfrac;
	    }
	  else QGPfrac = 1.;


	  if(pCutPropToT)
	    {
	      pCut = 4.*T; //standard: 4*T
	    }

	  //testing:
	  //cout << "pCut=" << pCut << endl;
	  //QGPfrac=1.;

	  vx = hydroInfo.vx;
	  vy = hydroInfo.vy;
	  vz = hydroInfo.vz;

	  if (!trackHistory)
	    {
	      // warn if out of range and move on to next parton (should not happen!)
	      if ( ix < 0 || ix >= ixmax ) 
		{
		  cout << "WARNING - x out of range" << endl;
		  plist[0]->at(i).frozen(1);
		  continue;
		}
	      if ( iy < 0 || iy >= ixmax )
		{
		  cout << "WARNING - y out of range" << endl;
		  plist[0]->at(i).frozen(1);
		  continue;
		}
	      if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4) && ( iz < 0 || iz >= izmax ) ) 
		{
		  cout << "WARNING - z out of range" << endl;
		  plist[0]->at(i).frozen(1);
		  continue;
		}
	    }

	  // ************************************************************************************************* output T
	  //if ( i == 0 ) cout << t << " " << tau << " " << x << " " << y << " " << z << " " 
	  //	     << T << " " << QGPfrac << " " << vecp.e() << endl;
	  // ************************************************************************************************* output T
	  
	  if (T<hydroTfinal && hydroWhichHydro == 2)                 // for 2D hydro evolutions stop at T_final (can change that)
	    {
	      plist[0]->at(i).frozen(1);
	      plist[0]->at(i).tFinal(t);
	      lengthTraveled = sqrt((x-plist[0]->at(i).xini())*(x-plist[0]->at(i).xini())
				     +(y-plist[0]->at(i).yini())*(y-plist[0]->at(i).yini())
				     +(z-plist[0]->at(i).zini())*(z-plist[0]->at(i).zini()));
	      mfp = lengthTraveled/NCol;
	      continue;                                             // do not evolve if T<T_c
	    }

	  // check the fraction of QGP at the current position (and time)
	  if (QGPfrac == 0. && (hydroWhichHydro == 3 || hydroWhichHydro == 1 || hydroWhichHydro == 4 || hydroWhichHydro == 5) )
	    {
	      plist[0]->at(i).frozen(1);
	      plist[0]->at(i).tFinal(t);
	      lengthTraveled = sqrt((x-plist[0]->at(i).xini())*(x-plist[0]->at(i).xini())
				     +(y-plist[0]->at(i).yini())*(y-plist[0]->at(i).yini())
				     +(z-plist[0]->at(i).zini())*(z-plist[0]->at(i).zini()));
	      mfp = lengthTraveled/NCol;
	      if(trackHistory && abs(z)<0.5)
		{
		  Tt->push_back(t);
		  Tx->push_back(x);
		  Ty->push_back(y);
		  Tz->push_back(z);
		  TE->push_back(p);
		  TQGPfrac->push_back(QGPfrac);
		  Tpx->push_back(vecp.px());
		  Tpy->push_back(vecp.py());
		  TdEdt->push_back(totalDeltaE/dtfm);
		  Tdpxdt->push_back(totalDeltapx/dtfm);
		  Tdpydt->push_back(totalDeltapy/dtfm);
		  Tdpzdt->push_back(totalDeltapz/dtfm);
		}
   	      continue;                                             // do not evolve if there is no QGP here
	    }
	  else if (QGPfrac < 1. && (hydroWhichHydro == 3 || hydroWhichHydro == 1 || hydroWhichHydro == 4 || hydroWhichHydro == 5))
	      // if QGP fraction is < 1, randomly decide to go on 
	    {
	      double rn=random->genrand64_real1();
	      //cout << "QGPfrac=" << QGPfrac << " random#=" << rn << endl;
	      if ( rn > QGPfrac )            // skip evolution with probability (1-QGPfrac)
		{
		  if(trackHistory  && abs(z)<0.5)
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
   		  continue;
		}
	    }

	  px = vecp.px();                                           // the initial momentum vector
	  py = vecp.py();
	  pz = vecp.pz();

	  beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
	  cosPhi = (px*vx + py*vy + pz*vz)/(vecp.pAbs()*beta);      // angle between p and flow vel.
	  gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
	  
	  //cout << "gamma=" << gamma << endl;
	  pRest = p * gamma * (1.-beta*cosPhi);                     // boost p to fluid cell's rest frame
	  
	  if(pRest/T<=AMYpCut)
	    {
	      if ( trackPartons )
		{
		  if (id==21) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg2 << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }		  
		  if (abs(id)<4) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracks2 << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		}
	      if(trackHistory && abs(z)<0.5  )
		{
		  Tt->push_back(t);
		  Tx->push_back(x);
		  Ty->push_back(y);
		  Tz->push_back(z);
		  TE->push_back(p);
		  TQGPfrac->push_back(QGPfrac);
		  Tpx->push_back(vecp.px());
		  Tpy->push_back(vecp.py());
		  TdEdt->push_back(totalDeltaE/dtfm);
		  Tdpxdt->push_back(totalDeltapx/dtfm);
		  Tdpydt->push_back(totalDeltapy/dtfm);
		  Tdpzdt->push_back(totalDeltapz/dtfm);
		}
	      continue;
	    }
	  
	  pxRest = -vx*gamma*p 
	    + (1.+(gamma-1.)*vx*vx/(beta*beta))*px 
	    + (gamma-1.)*vx*vy/(beta*beta)*py
	    + (gamma-1.)*vx*vz/(beta*beta)*pz;
	  pyRest = -vy*gamma*p 
	    + (1.+(gamma-1.)*vy*vy/(beta*beta))*py 
	    + (gamma-1.)*vx*vy/(beta*beta)*px
	    + (gamma-1.)*vy*vz/(beta*beta)*pz;
	  pzRest = -vz*gamma*p 
	    + (1.+(gamma-1.)*vz*vz/(beta*beta))*pz 
	    + (gamma-1.)*vx*vz/(beta*beta)*px
	    + (gamma-1.)*vy*vz/(beta*beta)*py;
	  
	  
	  vecpRest.px(pxRest);
	  vecpRest.py(pyRest);
	  vecpRest.pz(pzRest);
	  
	  cosPhiRest = (pxRest*vx + pyRest*vy + pzRest*vz)/(pRest*beta); // angle between pRest and flow vel.
	}
      else // the fixed temperature case:
	{
	  t = it*dtfm;        //had a +tau0 here. now let partons travel doing their vacuum shower stuff first...
	  plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
	  beta = 0.;
	  gamma = 1.; 
	  cosPhiRest=1.;
	  pRest = p;
	  T = fixedT;
	  vecpRest=vecp;
	  plist[0]->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
	  plist[0]->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
	  plist[0]->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
	  x = plist[0]->at(i).x();                                  // x value of position in [fm]
	  y = plist[0]->at(i).y();                                  // y value of position in [fm]
	  z = plist[0]->at(i).z();                                  // z value of position in [fm]
	  if(pRest/T<=AMYpCut) continue; // do not evolve partons with momenta below this scale 
	}	  
      
      boostBack = gamma*(1.+beta*cosPhiRest);                   // boost factor back from rest- to lab-frame 

      lengthTraveled = sqrt((x-plist[0]->at(i).xini())*(x-plist[0]->at(i).xini())
			    +(y-plist[0]->at(i).yini())*(y-plist[0]->at(i).yini())
			    +(z-plist[0]->at(i).zini())*(z-plist[0]->at(i).zini()));
      
      // get total probabilities for radiative processes
      if (doRadiative == 1) 
	{
	  n = rates->integrate(pRest/T,T,alpha_s, Nf, Ldependence, lengthTraveled);          // n stores the areas under the prob. distr.
      	  // Warn if probability for any process is larger than 1:
	  if (dt*n.Gamma/gamma>1. && pRest/T>AMYpCut) 
	    cout << "WARNING: probability to emit during one time step " << dt*n.Gamma/gamma 
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*n.Gamma_ggg/gamma >1. && pRest/T>AMYpCut) 
	    cout << "WARNING: probability ggg during one time step " << dt*n.Gamma_ggg/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << ", T=" << T << endl;
	  if (dt*n.Gamma_gqq/gamma>1. && pRest/T>AMYpCut) 
	    cout << "WARNING: probability gqq during one time step " << dt*n.Gamma_gqq/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*n.Gamma_em/gamma>1. && pRest/T>AMYpCut) 
	    cout << "WARNING: probability qqgamma during one time step " << dt*n.Gamma_em/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	}
      
      // get total probabilities for elastic processes
      if (doElastic == 1)
	{
	  qqRate = elastic->totalRate(pRest, T, alpha_s, Nf, 1, multElasticRate); 
	  gqRate = qqRate*9./4.;
	  qgRate = elastic->totalRate(pRest, T, alpha_s, Nf, 3, multElasticRate); 
	  ggRate = qgRate*9./4.; 
	  conversionqg = rates->Gammaqg( pRest, T, alpha_s );
	  conversiongq = rates->Gammagq( pRest, T, alpha_s, Nf );
	  if ( photonSwitch == 1 ) 
	    {
	      conversionqgamma = rates->Gammaqgamma( pRest, T, alpha_s );
	      if ( abs(id) == 1 ) conversionqgamma*=(4./9.); // multiplying by (ef/e)^2
	      else if ( abs(id) == 2 || abs(id) == 3 ) conversionqgamma*=(1./9.); // multiplying by (ef/e)^2
	    }
	  else 
	    conversionqgamma = 0.;
	  // Warn if probability for any process is larger than 1:
	  if (dt*qqRate/gamma>1.) 
	    cout << "WARNING: probability qq during one time step " << dt*qqRate/gamma 
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*gqRate/gamma >1.) 
	    cout << "WARNING: probability gq during one time step " << dt*gqRate/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*qgRate/gamma>1.) 
	    cout << "WARNING: probability qg during one time step " << dt*qgRate/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*ggRate/gamma>1.) 
	    cout << "WARNING: probability gg during one time step " << dt*ggRate/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*conversionqg/gamma>1.) 
	    cout << "WARNING: probability conversion q->g during one time step " << dt*conversionqg/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*conversionqgamma/gamma>1.) 
	    cout << "WARNING: probability conversion q->gamma during one time step " << dt*conversionqgamma/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	  if (dt*conversiongq/gamma>1.) 
	    cout << "WARNING: probability conversion g->q during one time step " << dt*conversiongq/gamma
		 << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
	}

      // for testing what happens if quark rates are similar to gluon rates
      if (quarksEqualGluons==2)
	{
	  n.Gamma_ggg*=4./9.;
	  n.Gamma_gqq*=4./9.;
	  gqRate = gqRate*4./9.;
	  ggRate = ggRate*4./9.;
	}
      else if (quarksEqualGluons==1)
	{
	  n.Gamma*=9./4.;
	  qqRate = qqRate*9./4.;
	  qgRate = qgRate*9./4.;
	}
      // this is to set the quark rates to be similar to gluon rates.
  
      if ( abs(id) > 0 && abs(id) < 4 ) // if parton is a quark (u, d, or s), let it evolve like a quark
	{
	  // in the following decide what to do in this time step
	  if (doRadiative == 1 && doElastic == 1 && pRest/T>AMYpCut ) //radiative + elastic 
	    {
	      if ( abs(id) == 1 ) // multiply by (e_f/e)^2
		n.Gamma_em*=4./9.;
	      else
		n.Gamma_em*=1./9.;

	      totalQuarkProb = dt/gamma*(n.Gamma)
		+dt/gamma*(qqRate)+dt/gamma*(qgRate)
		+dt/gamma*conversionqg+dt/gamma*conversionqgamma; // dt/gamma is dt_rest-frame

	      if ( photonSwitch == 1 ) 
		{
		  totalQuarkProb += dt/gamma*(n.Gamma_em);
		}

	      if (totalQuarkProb>1) cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
					 << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
	      
	      if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
		{
		  NCol+=1;
		  if ( random->genrand64_real1() < dt/gamma*(n.Gamma)/totalQuarkProb ) 
		    {
		      radiate = 1;
		    }
		  else if ( photonSwitch == 1 && random->genrand64_real1() < dt/gamma*n.Gamma_em/(totalQuarkProb-dt/gamma*(n.Gamma)))
		    {
		      radiatePhoton = 1;
		    }
		  else 
		    {
		      if ( random->genrand64_real1() < conversionqg/(qqRate+qgRate+conversionqg+conversionqgamma) )
			{
			  convert_quark = 1;
			}
		      else if ( random->genrand64_real1() < conversionqgamma/(qqRate+qgRate+conversionqgamma) )
			{
			  convert_quark_to_gamma = 1;
			}
		      else
			{
			  if ( random->genrand64_real1() < qqRate/(qqRate+qgRate) )
			    {
			      elastic_qq = 1;
			      elastic_qg = 0;
			    }
			  else
			    {
			      elastic_qq = 0;
			      elastic_qg = 1;
			    }
			}
		    }
		}
	      else 
		{
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			} 
			  foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
   		  continue;
		}
	    }
	  else if (doRadiative == 1 && doElastic == 0 && pRest/T>AMYpCut) //radiative
	    {
	      totalQuarkProb = dt/gamma*(n.Gamma);
	      
	      if ( photonSwitch == 1 ) 
		{
		  totalQuarkProb += dt/gamma*(n.Gamma_em);
		}
	      
	      if (totalQuarkProb>1) cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
					 << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
	      
	      if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
		{
		  NCol+=1;
		  if ( random->genrand64_real1() < dt/gamma*(n.Gamma)/totalQuarkProb ) 
		    {
		      radiate = 1;
		    }
		  else if ( photonSwitch == 1 )
		    {
		      radiatePhoton = 1;
		    }
		}
	      else 
		{
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
		  continue;
		}
	    }
	  else if (doRadiative == 0 && doElastic == 1 && pRest/T>AMYpCut) //elastic
	    {
	      totalQuarkProb = dt/gamma*(qqRate)+dt/gamma*(qgRate)+dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
	      // dt/gamma is dt_rest-frame
	      
	      if (totalQuarkProb>1) cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
					 << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
	      if ( random->genrand64_real1() < totalQuarkProb )  // check if something happens
		{
		  NCol+=1;
		  if ( random->genrand64_real1() < conversionqg/(qqRate+qgRate+conversionqg+conversionqgamma) )
		    {
		      convert_quark = 1;
		    }
		  else if ( random->genrand64_real1() < conversionqgamma/(qqRate+qgRate+conversionqgamma) )
		    {
		      convert_quark_to_gamma = 1;
		    }
		  else 
		    {
		      if ( random->genrand64_real1() < qqRate/(qqRate+qgRate) )
			{
			  elastic_qq = 1;
			  elastic_qg = 0;
			}
		      else
			{
			  elastic_qq = 0;
			  elastic_qg = 1;
			}
		    }
		}
	      else 
		{
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
		      foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
		  continue;
		}
	    }
	  else continue;

	  // now do what has been decided before
	  if( radiate == 1 ) // see if emission happens
	    {
	      if(pRest/T>AMYpCut) // do not evolve partons with momenta below this scale 
		{
		  // do process 1, q->qg
		  f = rates->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 1); 
		  p = (pRest - f.y*T)*boostBack;                // quark's new momentum in the lab frame
		  newvecp=p/vecp.pAbs()*vecp;                   // new quark momentum
		  if(trackHistory && f.y*T<pCut) 
		    {
		      totalDeltaE+=newvecp.pAbs()-vecp.pAbs();  // add to total energy change in this time step 
		      totalDeltapx+=newvecp.px()-vecp.px();  // add to total energy change in this time step 
		      totalDeltapy+=newvecp.py()-vecp.py();  // add to total energy change in this time step 
		      totalDeltapz+=newvecp.pz()-vecp.pz();  // add to total energy change in this time step 
		    }
		  
		  if (p>=0) plist[0]->at(i).p(newvecp);         // change quark's momentum to the new one
		  //plist[0]->at(i).splits(plist[0]->at(i).splits()+1);
		  if (f.y*T>pCut)                      // if new parton is kept (f.y*T > threshold [in GeV])
		    {
		      newOne.p(f.y*T*vecp.px()/vecp.pAbs()*boostBack,  // emitted gluon's momentum
			       f.y*T*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
			       f.y*T*vecp.pz()/vecp.pAbs()*boostBack);
		      newOne.pini(newOne.p());
		      col = plist[0]->at(i).col();              // color of original quark  
		      acol = plist[0]->at(i).acol();            // anti-color
		      if (col!=0)                               // if we had a quark
			{
			  plist[0]->at(i).col(10000000+counter);   // set color to new color
			  //plist[0]->at(i).acol(0);              // set anti-color to zero
			  newOne.col(col);                      // set new gluon's color to quark's original color
			  newOne.acol(10000000+counter);           // set new gluon's anti-color to quark's new color
			}
		      else if (acol!=0)                         // if we had an anti-quark
			{
			  plist[0]->at(i).col(0);               // set color to zero
			  //plist[0]->at(i).acol(10000000+counter);  // set anti-color to new color
			  newOne.col(10000000+counter);            // set new particle's color
			  newOne.acol(acol);                    // set new particle's anti-color   
			}
		      newOne.id(21);                            // emitted parton is a gluon
		      //newOne.splits(0);
		      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
		      newOne.mass(0.);
		      newOne.frozen(0);
		      newOne.x(x);                              // set the new parton's initial position
		      newOne.y(y);
		      newOne.z(z);
		      newOne.tini(t);
		      newOne.xini(x);                              // set the new parton's initial position
		      newOne.yini(y);
		      newOne.zini(z);
		      newOne.eventNumber(plist[0]->at(i).eventNumber());
		      newOne.source(2);                         // AMY parton
		      plist[0]->push_back(newOne);              // add the gluon to the list of partons if (f.y>AMYpCut)
		      if ( trackPartons) 
			{
			  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
			  if(abs(rapidity)<=1000.)
			    {
			      fouttracks << t << " " << x << " " << y << " " << z << " " << f.y*T << endl;
			    }
			  foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << f.y*T << ", " << endl;
			}
		    }
		}
	    }
	  else if( radiatePhoton == 1 ) // see if photon emission happens
	    {
	      if(pRest/T>0.01) // do not evolve partons with momenta below this scale 
		{
		  cout << "The loop on line 752 of MARTINI.cpp has been entered..." << endl;
		  cout << pRest/T << endl;
		  // do process 4, q->qgamma
		  f = rates->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 4); 
		  newOne.p(f.y*T*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
			   f.y*T*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
			   f.y*T*vecp.pz()/vecp.pAbs()*boostBack);
		  cout << "f.y = " << f.y << endl;
		  newOne.pini(newOne.p());
		  newOne.id(22);                            // emitted parton is a photon
		  //newOne.splits(0);
		  newOne.elasticCollisions(0);              // set initial no. of el colls to zero
		  newOne.mass(0.);
		  newOne.frozen(1);                         // photons do not interact anymore
		  newOne.x(x);                              // set the new parton's initial position
		  newOne.y(y);
		  newOne.z(z);
		  newOne.tini(t);
		  newOne.xini(x);                              // set the new parton's initial position
		  newOne.yini(y);
		  newOne.zini(z);
		  newOne.eventNumber(plist[0]->at(i).eventNumber());
		  newOne.source(2);                         // AMY photon
		  plist[0]->push_back(newOne);              // add the photon to the list of partons
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		}
	    }
	  else if ( convert_quark == 1 && pRest/T>AMYpCut ) // do q->g
	    {
	      col = plist[0]->at(i).col();              // color of original quark  
	      acol = plist[0]->at(i).acol();            // anti-color
	      if (col!=0)                               // if we had a quark
		{
		  plist[0]->at(i).acol(30000000+counter);  // add a new anti-color for the gluon
		  newOne.col(30000000+counter);            // create new thermal gluon
		  newOne.acol(40000000+counter);                        
		}
	      else if (acol!=0)
		{
		  plist[0]->at(i).col(30000000+counter);   // add a new color for the gluon
		  newOne.acol(30000000+counter);           // create new thermal gluon
		  newOne.col(40000000+counter);                        
		}

	      plist[0]->at(i).id(21);                   // convert to gluon
	      plist[0]->at(i).source(1);                // Conversion gluon
              newOne.id(21);                            // add a thermal gluon
 	      newOne.mass(0.); 
 	      newOne.frozen(0);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
	      newOne.tini(t);
 	      newOne.xini(x);                              // set the new parton's initial position
 	      newOne.yini(y);
 	      newOne.zini(z);
 	      newOne.p(random->thermal(T, -1));         // sample momentum from Bose distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
	      plist[0]->push_back(newOne); 

	      if (col!=0)                               // if we had a quark
		{
		  newOne.id(id);                        // add a thermal quark of same flavor
		  newOne.col(40000000+counter);            // attach the color string to the thermal quark
		  newOne.acol(0);                        
		}
	      else if (acol!=0)                         // if we had an anti-quark
		{
		  newOne.id(id);                        // add a thermal quark of same flavor
		  newOne.acol(40000000+counter);           // attach the color string to the thermal anti-quark
		  newOne.col(0);                        
		}
	      if ( abs(id) < 3) 
		newOne.mass(0.33);

	      else 
		newOne.mass(0.5);
	      newOne.frozen(0);
	      newOne.x(x);                              // set the new parton's initial position
	      newOne.y(y);
	      newOne.z(z);
	      newOne.tini(t);
	      newOne.xini(x);                              // set the new parton's initial position
	      newOne.yini(y);
	      newOne.zini(z);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
	      plist[0]->push_back(newOne); 
	      if ( trackPartons) 
		{
		  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		  if(abs(rapidity)<=1000.)
		    {
		      fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
		    }
		      foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		}
	      if(trackHistory  && abs(z)<0.5  )
		{
		  Tt->push_back(t);
		  Tx->push_back(x);
		  Ty->push_back(y);
		  Tz->push_back(z);
		  TE->push_back(p);
		  TQGPfrac->push_back(QGPfrac);
		  Tpx->push_back(vecp.px());
		  Tpy->push_back(vecp.py());
		  TdEdt->push_back(totalDeltaE/dtfm);
		  Tdpxdt->push_back(totalDeltapx/dtfm);
		  Tdpydt->push_back(totalDeltapy/dtfm);
		  Tdpzdt->push_back(totalDeltapz/dtfm);
		}
	      continue;                                 //prevent the new gluon from interacting again in this time step 
	    }
	  else if ( convert_quark_to_gamma == 1 && pRest/T>0.01 ) // do q->gamma
	    {
	      cout << "The if on line 879 of MARTINI.cpp has been entered..." << endl;
	      cout << pRest/T << endl;
	      col = plist[0]->at(i).col();              // color of original quark  
	      acol = plist[0]->at(i).acol();            // anti-color
	      if (col!=0)                               // if we had a quark
		{
		  plist[0]->at(i).acol(0);              // make color neutral photon
		  plist[0]->at(i).col(0);               // make color neutral photon
		  newOne.col(col);                      // attach the color string to the thermal gluon
		  newOne.acol(60000000+counter);                        
		}
	      else if (acol!=0)
		{
		  plist[0]->at(i).acol(0);              // make color neutral photon
		  plist[0]->at(i).col(0);               // make color neutral photon
		  newOne.acol(acol);                    // attach the color string to the thermal gluon
		  newOne.col(60000000+counter);                        
		}

	      plist[0]->at(i).id(22);                   // convert to photon
	      plist[0]->at(i).source(1);                // Conversion photon
              newOne.id(21);                            // add a thermal gluon
 	      newOne.mass(0.); 
 	      newOne.frozen(0);
 	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
	      newOne.tini(t);
 	      newOne.xini(x);                              // set the new parton's initial position
 	      newOne.yini(y);
 	      newOne.zini(z);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.p(random->thermal(T, -1));         // sample momentum from Bose distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
 	      plist[0]->push_back(newOne); 

	      if (col!=0)                               // if we had a quark
		{
		  newOne.id(id);                        // add a thermal quark of same flavor as initial quark
		  newOne.col(60000000+counter);            // attach the color string to the thermal quark
		  newOne.acol(0);                        
		}
	      else if (acol!=0)                         // if we had an anti-quark
		{
		  newOne.id(id);                        // add a thermal anti-quark of same flavor as initial anti-quark
		  newOne.acol(60000000+counter);           // attach the color string to the thermal anti-quark
		  newOne.col(0);                        
		}
	      if ( abs(id) < 3) 
		newOne.mass(0.33);
	      else 
		newOne.mass(0.5);
	      newOne.frozen(0);
	      newOne.x(x);                              // set the new parton's initial position
	      newOne.y(y);
	      newOne.z(z);
	      newOne.tini(t);
	      newOne.xini(x);                           // set the new parton's initial position
	      newOne.yini(y);
	      newOne.zini(z);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
	      plist[0]->push_back(newOne); 
	      if ( trackPartons) 
		{
		  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		  if(abs(rapidity)<=1000.)
		    {
		      fouttracks << t << " " << x << " " << y << " " << z << " " << 0. << endl;
		    }
		      foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		}
	      if(trackHistory  && abs(z)<0.5  )
		{
		  Tt->push_back(t);
		  Tx->push_back(x);
		  Ty->push_back(y);
		  Tz->push_back(z);
		  TE->push_back(p);
		  TQGPfrac->push_back(QGPfrac);
		  Tpx->push_back(vecp.px());
		  Tpy->push_back(vecp.py());
		  TdEdt->push_back(totalDeltaE/dtfm);
		  Tdpxdt->push_back(totalDeltapx/dtfm);
		  Tdpydt->push_back(totalDeltapy/dtfm);
		  Tdpzdt->push_back(totalDeltapz/dtfm);
		}
	      continue;                                 //prevent the new gluon from interacting again in this time step 
	    }
	  else if ( elastic_qq == 1 && pRest/T>AMYpCut ) // do qq->qq
	    {
	      omega = elastic->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 1); // this is omega/T
	      if( transferTransverseMomentum == 1 && omega < 800.)
		q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s, Nf, random, import, 1); // this is q/T
	      else 
		q=omega;
	    }
	  else if ( elastic_qg == 1 && pRest/T>AMYpCut ) // do qg->qg
	    {
	      omega = elastic->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 3); // this is omega/T 
	      if( transferTransverseMomentum == 1 && omega < 800.)
		{ 
		  q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s, Nf, random, import, 3); // this is q/T
		} 
	      else 
		q = omega;
	    }
	  if ( elastic_qq == 1 || elastic_qg == 1 ) // in both cases ( qq->qq and qg->qg ) change the momentum
	    {
	      double cosTheta_pq = (-omega*omega*T*T+2.*pRest*omega*T+q*q*T*T)/(2.*pRest*q*T);
	      //cout << "cosTheta_pq=" << cosTheta_pq << endl;
	      double qt = q*T*sqrt(1.-cosTheta_pq*cosTheta_pq);                // transverse momentum transfer
	      qtTot+=qt*qt;
	      //cout << "p=" << pRest << endl;
	      //cout << "omega=" << omega << endl;
	      //cout << "q=" << q << endl;
	      //cout << "qt^2=" << qtTot << endl;
	      
	      newvecp=elastic->getNewMomentum(vecpRest, omega*T, q*T, random); // takes everything in GeV

	      //NCol+=1;
	      //cout << "MARTINI got the newvecp=" << newvecp;
	      
	      // now I need to boost the full vector back to the lab frame
	      // angle between new pRest and flow vel.:
	      if (beta!=0) cosPhiRestEl = (newvecp.px()*(vx) + newvecp.py()*(vy) + newvecp.pz()*(vz))/(newvecp.pAbs()*beta); 
	      else cosPhiRestEl=1.;
	   
	      //cout << "cosPhiRestEl=" << cosPhiRestEl << endl;
	      
	      p = newvecp.pAbs() * gamma * (1.+beta*cosPhiRestEl);    // boost p back to lab frame
	      //cout << "p=" << p  << endl;
	      
	      if (beta!=0)
		{
		  newpx = vx*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecp.px() 
		    + (gamma-1.)*vx*vy/(beta*beta)*newvecp.py()
		    + (gamma-1.)*vx*vz/(beta*beta)*newvecp.pz();
		  newpy = vy*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecp.py()
		    + (gamma-1.)*vx*vy/(beta*beta)*newvecp.px()
		    + (gamma-1.)*vy*vz/(beta*beta)*newvecp.pz();
		  newpz = vz*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecp.pz()
		    + (gamma-1.)*vx*vz/(beta*beta)*newvecp.px()
		    + (gamma-1.)*vy*vz/(beta*beta)*newvecp.py();
		  
		//   double newpx2, newpy2, newpz2;
		  
// 		  newpx2 = vx*gamma*vecpRest.pAbs() 
// 		    + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecpRest.px() 
// 		    + (gamma-1.)*vx*vy/(beta*beta)*vecpRest.py()
// 		    + (gamma-1.)*vx*vz/(beta*beta)*vecpRest.pz();
// 		  newpy2 = vy*gamma*vecpRest.pAbs() 
// 		    + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecpRest.py()
// 		    + (gamma-1.)*vx*vy/(beta*beta)*vecpRest.px()
// 		    + (gamma-1.)*vy*vz/(beta*beta)*vecpRest.pz();
// 		  newpz2 = vz*gamma*vecpRest.pAbs() 
// 		    + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecpRest.pz()
// 		    + (gamma-1.)*vx*vz/(beta*beta)*vecpRest.px()
// 		    + (gamma-1.)*vy*vz/(beta*beta)*vecpRest.py();

		  
		//   if ( newpy>0.01 || newpy<-0.01  )
// 		    {
// 		      cout << "p: " << px << " " << py << " " << pz << endl;
// 		      cout << "v: " << vx << " " << vy << " " << vz << endl;
// 		      cout << "r: " << pxRest << " " << pyRest << " " << pzRest << endl;
// 		      cout << "ar: " << newvecp.px() << " " << newvecp.py() << " " << newvecp.pz() << endl;
// 		      cout << "ap: " << newpx << " " << newpy << " " << newpz << endl;
// 		      cout << "vy=" << vy << ", gamma=" << gamma << ", newvecp.pAbs()=" << newvecp.pAbs() << endl;
// 		      cout << " vy*gamma*newvecp.pAbs() = " <<  vy*gamma*newvecp.pAbs() << endl;
// 		      cout << "(1.+(gamma-1.)*vy*vy/(beta*beta))*newvecp.py() = " << (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecp.py() << endl;
// 		    }

		  newvecp.px(newpx);
		  newvecp.py(newpy);
		  newvecp.pz(newpz);
		  newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
			  // cout << "p=" << p << " e=" << newvecp.e() << endl; 

		}
	      if (p>=0) plist[0]->at(i).p(newvecp);           // change quark's momentum to the new one
	      plist[0]->at(i).elasticCollisions(plist[0]->at(i).elasticCollisions()+1); // increase number of elastic processes by one 
		
	      if(trackHistory && abs(newvecp.pAbs()-vecp.pAbs())<pCut) 
		{
		  totalDeltaE+=newvecp.pAbs()-vecp.pAbs();          // add to total energy change in this time step 
		  totalDeltapx+=newvecp.px()-vecp.px();     
		  totalDeltapy+=newvecp.py()-vecp.py();    
		  totalDeltapz+=newvecp.pz()-vecp.pz();     
		}
	
	      if ( trackPartons) 
		{
		  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		  if(abs(rapidity)<=1000.)
		    {
		      fouttracks << t << " " << x << " " << y << " " << z << " " << omega*T << endl;
		    }      
		      foutpovray << t << ", " << x << ", " << y << ", " << z << ", " << omega*T << ", " << endl;
		}
	    }
	}

      id = plist[0]->at(i).id();

      if ( id == 21 ) // if parton is a gluon, let it evolve like a gluon
	{
	  // in the following decide what to do in this time step
	  if (doRadiative == 1 && doElastic == 1 && pRest/T>AMYpCut )  //radiative + elastic
	    {
	      totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq)
		+dt/gamma*(gqRate)+dt/gamma*(ggRate)+dt/gamma*conversiongq; 
	
	      if (totalGluonProb>1) cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
					 << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
	      if ( random->genrand64_real1() < totalGluonProb ) // check if something happens with the quark
		{
		  if ( random->genrand64_real1() < (dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq))/totalGluonProb ) 
		    {
		      if ( random->genrand64_real1() < n.Gamma_ggg/(n.Gamma_ggg+n.Gamma_gqq) ) 
			radiate_ggg=1;
		      else 
			radiate_gqq=1;
		    }
		  else 
		    {
		      if ( random->genrand64_real1() < conversiongq/(conversiongq+gqRate+ggRate) )
			{
			  convert_gluon = 1;
			}
		      else
			{
			  if ( random->genrand64_real1() < gqRate/(gqRate+ggRate) )
			    {
			      elastic_gq = 1;
			    }
			  else
			    {
			      elastic_gg = 1;
			    }
			}
		    }
		}
	      else 
		{
		  if ( trackPartons ) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
		  continue;
		}
	    }
	  else if (doRadiative == 1 && doElastic == 0 && pRest/T>AMYpCut) //radiative
	    {
	      totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq); 
	
	      if (totalGluonProb>1) cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
					 << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
	      if ( random->genrand64_real1() < totalGluonProb ) // check if something happens with the quark
		{
		  if ( random->genrand64_real1() < n.Gamma_ggg/(n.Gamma_ggg+n.Gamma_gqq) ) 
		    radiate_ggg=1;
		  else 
		    radiate_gqq=1;
		}
	      else 
		{
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}		
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }  
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
		  continue;
		}
	    }
	  else if (doRadiative == 0 && doElastic == 1 && pRest/T>AMYpCut) //elastic
	    {
	      totalGluonProb = dt/gamma*(gqRate)+dt/gamma*(ggRate)+dt/gamma*conversiongq; // dt/gamma is dt_rest-frame
	      
	      if (totalGluonProb>1) cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
					 << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
	      if ( random->genrand64_real1() < totalGluonProb )  // check if something happens
		{
		  if ( random->genrand64_real1() < conversiongq/(conversiongq+gqRate+ggRate) )
		    {
		      convert_gluon = 1;
		    }
		  else
		    {
		      if ( random->genrand64_real1() < gqRate/(gqRate+ggRate) )
			{
			  elastic_gq = 1;
			}
		      else
			{
			  elastic_gg = 1;
			}
		    }
		}
	      else 
		{
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg << t << " " << x << " " << y << " " << z << " " << 0. << endl;
			}
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		    }
		  if(trackHistory  && abs(z)<0.5  )
		    {
		      Tt->push_back(t);
		      Tx->push_back(x);
		      Ty->push_back(y);
		      Tz->push_back(z);
		      TE->push_back(p);
      		      TQGPfrac->push_back(QGPfrac);
		      Tpx->push_back(vecp.px());
		      Tpy->push_back(vecp.py());
		      TdEdt->push_back(totalDeltaE/dtfm);
		      Tdpxdt->push_back(totalDeltapx/dtfm);
		      Tdpydt->push_back(totalDeltapy/dtfm);
		      Tdpzdt->push_back(totalDeltapz/dtfm);
		    }
		  continue;
		}
	    }
	  else continue;
	  
	  // now do what has been decided before
	  if( radiate_ggg == 1) // g -> gg
	    {
	      if(pRest/T>AMYpCut) 
		{
		  f = rates->findValuesRejection
		    (pRest/T, T, alpha_s, Nf, random, import, 3);   // do process 3
		  
		  cout.precision(6);
		  //cout << "------------------------------" << endl;
		  //double oldp = p;
		  //cout << "p_before=" << p << endl;
		  p = (pRest - f.y*T)*boostBack;
				
		  /*
		    cout << "p_after/p=" << p/oldp << endl;
		    cout << "pRest_after/pRest=" << (pRest-f.y*T)/pRest << endl;
		    
		    cout << "gamma=" << gamma << endl;
		    cout << "beta=" << beta << endl;
		    cout << "cosPhi=" << cosPhi << endl;
		    cout << "ratio cor/fake=" << p/((pRest - f.y*T)/(gamma*(1.-beta*cosPhi))) << endl;
		    cout << "pt^2=" << pow2(p*vecp.px()/vecp.pAbs())+pow2(p*vecp.py()/vecp.pAbs()) << endl;
		    cout << "pt_alt^2=" << pow2(vecp.px()/vecp.pAbs()*((pRest - f.y*T)/(gamma*(1.-beta*cosPhi))))+
		    pow2((vecp.py()/vecp.pAbs()*((pRest - f.y*T)/(gamma*(1.-beta*cosPhi))))) << endl;
		  */
		  
		  //cout << "ratio cor/fake=" << (gamma*(1.+beta*cosPhi))/( 1./(gamma*(1.-beta*cosPhi))) << endl;
		  newvecp=p/vecp.pAbs()*vecp;                       // new gluon momentum
		  if (p>=0) plist[0]->at(i).p(newvecp);                       // change gluon's momentum to the new one
		  //plist[0]->at(i).splits(plist[0]->at(i).splits()+1); // count the splits (one more)
		  //cout << "f.y=" << f.y << endl;
		  //cout << "f.y*gamma*(1+beta*cosPhi)=" << f.y*gamma*(1.+beta*cosPhi) << endl;
		  //cout << "f.y/(gamma*(1-beta*cosPhi))=" << f.y/(gamma*(1.-beta*cosPhi)) << endl;
		  if(trackHistory && f.y*T<pCut) 
		    {
		      totalDeltaE+=newvecp.pAbs()-vecp.pAbs();  // add to total energy change in this time step 
		      totalDeltapx+=newvecp.px()-vecp.px();  // add to total energy change in this time step 
		      totalDeltapy+=newvecp.py()-vecp.py();  // add to total energy change in this time step 
		      totalDeltapz+=newvecp.pz()-vecp.pz();  // add to total energy change in this time step 
		    }
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg << t << " " << x << " " << y << " " << z << " " << f.y*T << endl;
			}
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << f.y*T << ", " << endl;
		    }
		  if(f.y*T>pCut)
		    {
		      newOne.p(f.y*T*vecp.px()/vecp.pAbs()*boostBack, // emitted gluon's momentum (coll) 
			       f.y*T*vecp.py()/vecp.pAbs()*boostBack,
			       f.y*T*vecp.pz()/vecp.pAbs()*boostBack); 
		      newOne.pini(newOne.p());
		      //col = plist[0]->at(i).col();                  // the gluon's color
		      acol = plist[0]->at(i).acol();                // the gluon's anti-color
		      plist[0]->at(i).acol(20000000+counter);          // set the first new gluon's anti-color to a new one 
		      newOne.id(21);                                // second new particle is a gluon too
		      newOne.col(20000000+counter);                    // set the second gluon's color to the new color 
		      newOne.acol(acol);                            // set second gluon's anti-color to the original one
		      //newOne.splits(0);                         // hasn't split yet => zero
		      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
		      newOne.x(x);                                  // set the new parton's initial position
		      newOne.y(y);
		      newOne.z(z);
		      newOne.tini(t);
		      newOne.xini(x);                                  // set the new parton's initial position
		      newOne.yini(y);
		      newOne.zini(z);
		      newOne.mass(0.);                              // set second gluon's mass to zero
		      newOne.frozen(0);
		      newOne.eventNumber(plist[0]->at(i).eventNumber());
		      newOne.source(1);                         // All photons produced in evolve(...) have itsSource=1
		      plist[0]->push_back(newOne);                  // add the second gluon to the list of partons
		    }
		}
	    }
	  pRest = p * gamma * (1.-beta*cosPhi); // boost p to fluid cell's rest frame
	  if( radiate_gqq == 1 ) // g -> qq
	    {
	      if(pRest/T>AMYpCut) 
		{
		  f = rates->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 2); // do process 2
		  p = (pRest - f.y*T)*boostBack;
		 		  //  if(trackHistory) totalDeltaE-=f.y*T;          // add to total energy change in this time step 
		  // choose if it is a u-ubar, d-dbar or s-sbar pair:
		  double r = random->genrand64_real1();
		  // note that when using general N_f this has to be changed
		  if (r<0.33) newID=1;
		  else if  (r<0.66) newID=2;
		  else newID=3;
		  double mass;
		  if (newID<3) mass=0.33;                           // set the quark's mass
		  else mass = 0.5;                                  // (pythia needs that for fragmentation) 
		  newvecp=p/vecp.pAbs()*vecp;                       // new quark momentum
		  plist[0]->at(i).id(newID);                        // turn gluon into quarks 
		  // ########## now anti-quark has usually smaller momentum - take care of that
		  col = plist[0]->at(i).col();                      // gluon's color
		  acol = plist[0]->at(i).acol();                    // gluon's anti-color 
		  plist[0]->at(i).col(col);                         // set quark's color
		  plist[0]->at(i).acol(0);                          // set quark's anti-color to zero
		  plist[0]->at(i).mass(mass);                       // set quark's mass
		  //plist[0]->at(i).splits(plist[0]->at(i).splits()+1); // count the splits (one more)
		  if (p>=0) plist[0]->at(i).p(newvecp);             // change quark's momentum to the new one
		  newOne.p(f.y*T*vecp.px()/vecp.pAbs()*boostBack,   // second quark's momentum (collinear)
			   f.y*T*vecp.py()/vecp.pAbs()*boostBack,
			   f.y*T*vecp.pz()/vecp.pAbs()*boostBack); 
		  newOne.pini(newOne.p());
		  newOne.id(-newID);                                // second new particle is an anti-quark (minus-sign)
		  newOne.mass(mass);                                // set anti-quark's mass  
		  newOne.col(0);                                    // set anti-quark's color to zero
		  newOne.acol(acol);                                // set anti-quark's anti-color to gluon's anti-color
		  //newOne.splits(0);                               // this is new and hasn't split yet => zero  
		  newOne.elasticCollisions(0);              // set initial no. of el colls to zero
		  newOne.x(x);                                      // set the new parton's initial position
		  newOne.y(y);
		  newOne.z(z);
		  newOne.tini(t);
		  newOne.xini(x);                                      // set the new parton's initial position
		  newOne.yini(y);
		  newOne.zini(z);
		  newOne.frozen(0);
		  newOne.eventNumber(plist[0]->at(i).eventNumber());
		  newOne.source(1);                         // All photons produced in evolve(...) have itsSource=1
		  plist[0]->push_back(newOne);                      // add the second quark to the list of partons
		  if ( trackPartons) 
		    {
		      double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		      if(abs(rapidity)<=1000.)
			{
			  fouttracksg << t << " " << x << " " << y << " " << z << " " << f.y*T << endl;
			}
			  foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << f.y*T << ", " << endl;
		    }
		}
	    }
	  if ( convert_gluon == 1 && pRest/T>AMYpCut ) // do g->q
	    {
	      if ( random->genrand64_real1() < 0.5 ) 
		GluonToQuark = 1;
	      else 
		GluonToQuark = 0;

	      col = plist[0]->at(i).col();              // color of original gluon  
	      acol = plist[0]->at(i).acol();            // anti-color

	      if (GluonToQuark == 1)
		{
		  plist[0]->at(i).acol(0);              // make the gluon a quark
		  newOne.acol(acol);                    // let the thermal gluon have the gluon's previous anti-color
		  newOne.col(50000000+counter);
		  plist[0]->at(i).id(1);                  // convert to quark
		}
	      else
		{
		  plist[0]->at(i).col(0);               // make the gluon an anti-quark
		  newOne.col(col);                      // let the thermal gluon have the gluon's previous color
		  newOne.acol(50000000+counter);
		  plist[0]->at(i).id(-1);                 // convert to anti-quark
		}
	      
	      plist[0]->at(i).mass(0.33);               
	      	
	      newOne.id(21);                            // add a thermal gluon
 	      newOne.mass(0.); 
 	      newOne.frozen(0);
 	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
	      newOne.tini(t);
 	      newOne.xini(x);                              // set the new parton's initial position
 	      newOne.yini(y);
 	      newOne.zini(z);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.p(random->thermal(T, -1));         // sample momentum from Bose distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
	      newOne.source(1);                         // All photons produced in evolve(...) have itsSource=1
	      plist[0]->push_back(newOne); 
	 
	      double r = random->genrand64_real1();
	      if (r<0.33) newID=1;
	      else if  (r<0.66) newID=2;
	      else newID=3;
	      double mass;
	      if (newID<3) mass=0.33;                           // set the quark's mass
	      else mass = 0.5;                                  // (pythia needs that for fragmentation) 

	      if (GluonToQuark == 1)
		{
		  newOne.acol(50000000+counter);           // let the thermal anti-quark have the right anti-color
		  newOne.col(0);
		  newOne.id(-newID);                        // add a thermal anti-quark
		}
	      else
		{
		  newOne.col(50000000+counter);            // let the thermal quark have the right color
		  newOne.acol(0);
		  newOne.id(newID);                         // add a thermal quark
		}

	      newOne.mass(mass); 
 	      newOne.frozen(0);
 	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
	      newOne.tini(t);
 	      newOne.xini(x);                              // set the new parton's initial position
 	      newOne.yini(y);
 	      newOne.zini(z);
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.eventNumber(plist[0]->at(i).eventNumber());
	      newOne.source(1);                         // All photons produced in evolve(...) have itsSource=1
	      plist[0]->push_back(newOne); 
	      if ( trackPartons) 
		{
		  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		  if(abs(rapidity)<=1000.)
		    {
		      fouttracksg << t << " " << x << " " << y << " " << z << " " << 0. << endl;
		    }
		      foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << 0. << ", " << endl;
		}
	      if(trackHistory  && abs(z)<0.5  )
		{
		  Tt->push_back(t);
		  Tx->push_back(x);
		  Ty->push_back(y);
		  Tz->push_back(z);
		  TE->push_back(p);
		  TQGPfrac->push_back(QGPfrac);
		  Tpx->push_back(vecp.px());
		  Tpy->push_back(vecp.py());
		  TdEdt->push_back(totalDeltaE/dtfm);
		  Tdpxdt->push_back(totalDeltapx/dtfm);
		  Tdpydt->push_back(totalDeltapy/dtfm);
		  Tdpzdt->push_back(totalDeltapz/dtfm);
		}
	      continue;                                 //prevent the new gluon from interacting again in this time step 
	    }
	  else if ( elastic_gq == 1 && pRest/T>AMYpCut ) // do gq->gq
	    {
	      omega = elastic->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 2); // this is omega/T
	      if( transferTransverseMomentum == 1 && omega < 800.)
		q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s, Nf, random, import, 2); // this is q/T
	      else
		q = omega;
	    }
	  else if ( elastic_gg == 1 && pRest/T>AMYpCut ) // do gg->gg
	    {
	      omega = elastic->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 4); 
	      if( transferTransverseMomentum == 1 && omega < 800.)
		q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s, Nf, random, import, 4); // this is q/T
	      else
		q = omega;
	    }
	  if ( elastic_gq == 1 || elastic_gg == 1 ) // in both cases ( gq->gq and gg->gg ) change the momentum 
	    {
	      newvecp=elastic->getNewMomentum(vecpRest, omega*T, q*T, random);
	      //NCol+=1;
	      
	      if(trackHistory && abs(newvecp.pAbs()-vecp.pAbs())<pCut ) 
		{
		  totalDeltaE+=newvecp.pAbs()-vecp.pAbs();          // add to total energy change in this time step 
		  totalDeltapx+=newvecp.px()-vecp.px();     
		  totalDeltapy+=newvecp.py()-vecp.py();    
		  totalDeltapz+=newvecp.pz()-vecp.pz();     
		}
	      
	      if (beta!=0) cosPhiRestEl = (newvecp.px()*vx + newvecp.py()*vy + newvecp.pz()*vz)/(newvecp.pAbs()*beta); 
	      else cosPhiRestEl=1.;
	   
	      p = newvecp.pAbs() * gamma * (1.+beta*cosPhiRestEl);    // boost p back to lab frame
	      
	      if (beta!=0)
		{
		  newpx = vx*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecp.px() 
		    + (gamma-1.)*vx*vy/(beta*beta)*newvecp.py()
		    + (gamma-1.)*vx*vz/(beta*beta)*newvecp.pz();
		  newpy = vy*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecp.py()
		    + (gamma-1.)*vx*vy/(beta*beta)*newvecp.px()
		    + (gamma-1.)*vy*vz/(beta*beta)*newvecp.pz();
		  newpz = vz*gamma*newvecp.pAbs() 
		    + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecp.pz()
		    + (gamma-1.)*vx*vz/(beta*beta)*newvecp.px()
		    + (gamma-1.)*vy*vz/(beta*beta)*newvecp.py();
		  
		  newvecp.px(newpx);
		  newvecp.py(newpy);
		  newvecp.pz(newpz);
		  newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
		}
	      if (p>=0) plist[0]->at(i).p(newvecp);           // change quark's momentum to the new one
	      plist[0]->at(i).elasticCollisions(plist[0]->at(i).elasticCollisions()+1); // increase number of elastic processes by one 
	      if ( trackPartons ) 
		{
		  double rapidity = 0.5*log((p+vecp.pz())/(p-vecp.pz()));
		  if(abs(rapidity)<=1000.)
		    {
		      fouttracksg << t << " " << x << " " << y << " " << z << " " << omega*T << endl;
		    }
		      foutpovrayg << t << ", " << x << ", " << y << ", " << z << ", " << omega*T << ", " << endl;
		}
	    }
	}
      if(trackHistory && abs(z)<0.5 )
	{
	  Tt->push_back(t);
	  Tx->push_back(x);
	  Ty->push_back(y);
	  Tz->push_back(z);
	  TE->push_back(p);
	  TQGPfrac->push_back(QGPfrac);
	  Tpx->push_back(vecp.px());
	  Tpy->push_back(vecp.py());
	  TdEdt->push_back(totalDeltaE/dtfm);
	  Tdpxdt->push_back(totalDeltapx/dtfm);
	  Tdpydt->push_back(totalDeltapy/dtfm);
	  Tdpzdt->push_back(totalDeltapz/dtfm);
	}
    } // loop over partons
  if ( trackPartons) 
    {
      fouttracks.close();
      fouttracksg.close();
      fouttracks2.close();
      fouttracksg2.close();
      foutpovray.close();
      foutpovrayg.close();
    }
  //foutt.close();
  return counter;
}

void MARTINI::plotHydroData()
{
  HydroInfo hydroInfo;
  cout.precision(6);
  //int ixmax = 80;
  //int izmax = 80;
  int ixmax = 180;
  int izmax = 180;
  int ix, iy, iz;
  double x,y;
  double z, denominator, numerator, ecc;
  double tau, t;
  double T, Tx0, Ty0, vyx0, vxy0, Vx, Vy, Veta;
  
  ofstream foutT1("./output/Txy.dat",ios::out); 
  ofstream foutT("./output/TxN.dat",ios::out); 
  ofstream foutTy("./output/TyN.dat",ios::out); 
  //foutT <<  " t  " << " x " << " y " <<  " T(x,y=0) " << " T(x=0,y) " << endl;
  
  numerator = 0.;
  denominator = 0.;
  for(int it=0; it<1; it+=1)
      {
	t=0.4+0.2*it;
	for(int ix=1; ix<ixmax; ix++)
	  for(int iy=1; iy<ixmax; iy++)
	    {
	      x=ix*0.1-9.;
	      y=iy*0.1-9.;
	     
	      //z=iz*hydroDz-hydroZmax+hydroDz;
	      //t=sqrt(tau*tau+z*z);
	      z=0.;
	      //cout << "x=" << x << endl;
	      //cout << "y=" << y << endl;
	      //cout << "z=" << z << endl;
	      //cout << "tau=" << sqrt(t*t-z*z) << endl;
	      if (t*t<z*z || sqrt(t*t-z*z)<=0.1) T=0.;
	      else
		{
		  hydroInfo = hydroSetup->getHydroValues(x, y, z, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
							 hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice);
		  T = hydroInfo.T;
		}
	      foutT1 << t << ", " << x << ", " << y << ", " << z << ", " << T << ", " 
		     << hydroInfo.vx << ", " << hydroInfo.vy << ", " <<  endl;
	      if (y>2.99 && y<3.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (y>3.99 && y<4.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (y>4.99 && y<5.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (y>5.99 && y<6.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (y>6.99 && y<7.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (y>7.99 && y<8.01)
		{
		  foutT << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		  cout << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			<< hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (x<6.71 && x>6.69)
		{
		  foutTy << t << " " << x << " " << y << " " << z << " " << T << " " << pow(T,4.) << " "  
			 << hydroInfo.vx << " " << hydroInfo.vy << " " << hydroInfo.QGPfrac << " " <<  endl;
		}
	      if (it==0)
		{
		  numerator += (pow(T,4.))*(y*y-x*x);
		  //cout << "x=" << x << ", y=" << y << ", T=" << T << ", adding " << (pow(T,4.))*(y*y-x*x) << endl;
		  denominator += (pow(T,4.))*(y*y+x*x);
		}
	    }
      }
  ecc = numerator/denominator;
  cout << "initial eccentricity=" << ecc << endl;
  cout << "num=" << numerator << endl;
  cout << "den=" << denominator << endl;
  foutT.close();
  foutTy.close();
  foutT1.close();


  /* 
  // this part is to convert the full data file into one quarter of it - due to symmetry this is all we need
  string name, name2;
  ostringstream temp;
  string istr;

  for(int it=0; it<170; it++)
      {
	temp.str("");
	temp << it;
	istr = temp.str();
	name="../hydro/nonaka/T"+istr+".dat";
	name2="../hydro/nonaka/V"+istr+".dat";
	ofstream foutT(name.c_str(),ios::out); 
	ofstream foutV(name2.c_str(),ios::out); 
	tau=0.6+0.1*it; // 0 corresponds to 0.6 fm/c
	cout << "converting time step tau=" << tau << " fm/c" << endl;
	for(int iz=0; iz<izmax; iz++)
	  {
	    for(int iy=ixmax/2; iy<ixmax; iy++)
	      {
		for(int ix=ixmax/2; ix<ixmax; ix++)
		  {
		    x=ix*hydroDx-hydroXmax;
		    y=iy*hydroDx-hydroXmax;
		    z=iz*hydroDz-hydroZmax;
		    t=sqrt(tau*tau+z*z);
		    //cout << "tau=" << tau << endl;
		    //cout << "t=" << t << endl;
		    //cout << "x=" << x << endl;
		    //cout << "y=" << y << endl;
		    //cout << "z=" << z << endl;
		    
		    hydroInfo = hydroSetup->getHydroValues2(x, y, z, t, hydroXmax, hydroZmax, hydroTau0, 
		    hydroDx, hydroDz, hydroDtau, hydroWhichHydro, lattice);
		    
		    T = hydroInfo.T;
		    Vx = hydroInfo.vx;
		    Vy = hydroInfo.vy;
		    Veta = hydroInfo.veta; // the v_eta stored in hydroInfo only for this routine

		    //foutT << x << " " << y << " " << z << " " << T << " " << hydroInfo.QGPfrac << endl;
		    foutT << T << " " << hydroInfo.QGPfrac << endl;
		    //foutV << x << " " << y << " " << z << " " << Vx << " " << Vy << " " << Veta << endl;
		    foutV << Vx << " " << Vy << " " << Veta << endl;
		    //if ( iy == ixmax-1 ) foutT << endl;
		  }
	      }
	  }
	foutT.close();
	foutV.close();
      }
  */          
  exit(1);
}

void MARTINI::plotHydroDataXYZ()
{
  HydroInfo hydroInfo;
  cout.precision(6);
  //int ixmax = 80;
  //int izmax = 80;
  int ixmax = 38;
  int izmax = 160;
  int ix, iy, iz;
  double x,y;
  double z, denominator, numerator, ecc;
  double tau, t;
  double T, Tx0, Ty0, vyx0, vxy0, Vx, Vy, Veta;
  
  ofstream foutT("./output/Txyz.dat",ios::out); 
  //foutT <<  " t  " << " x " << " y " <<  " T(x,y=0) " << " T(x=0,y) " << endl;
  
  numerator = 0.;
  denominator = 0.;
  for(int it=0; it<60; it+=1)
      {
	t=0.4+0.25*it;
	cout << "t=" << t << endl;
	for(int ix=0; ix<ixmax; ix++)
	  for(int iy=0; iy<ixmax; iy++)
	    for(int iz=0; iz<izmax; iz++)
	    {
	      x=ix*0.5-9.5;
	      y=iy*0.5-9.5;
	      z=iz*0.25-20.;
	      //t=sqrt(tau*tau+z*z);
	      //cout << "x=" << x << endl;
	      //cout << "y=" << y << endl;
	      //cout << "z=" << z << endl;
	      //cout << "tau=" << sqrt(t*t-z*z) << endl;
	      if (t*t<z*z || sqrt(t*t-z*z)<=0.4) T=0.;
	      else
		{
		  hydroInfo = hydroSetup->getHydroValues(x, y, z, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
							 hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice);
		  T = hydroInfo.T;
		}
	      foutT << t << " " << x << " " << y << " " << z << " " << T << " " << endl;
	    }
      }
  foutT.close();
  exit(1);
}



void MARTINI::sampleTA()
{
  int A = static_cast<int>(glauber->nucleusA());

  //cout << "A=" << A << endl;
  ReturnValue rv, rv2;

  for (int i = 0; i < A; i++) // get all nucleon coordinates
    {
      rv = glauber->SampleTARejection(random);
      rv2 = glauber->SampleTARejection(random);
      nucleusA.push_back(rv);
      nucleusB.push_back(rv2);
    }
}  

int MARTINI::generateEvent(vector<Parton> *plist)
{
  if (fullEvent==1)
    {
      // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
      // then lay them on top of each other and sample the number of jet events in a tube with area=inelasticXSec
      
      //cout << " number of binary collisions = " << glauber->TAB() << endl;
      
      //cout << " A=" << glauber->nucleusA() << endl;
      double A;
      double Z;
      Parton parton;                               // new Parton object
      double b = glauberImpactParam;
      int posx;
      int posy;
      int n1=0;
      int n2=0;
      double r = 1.2*pow(glauber->nucleusA(),1./3.);

      nucleusA.clear();
      nucleusB.clear();
      sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons
      
      double cellLength = sqrt(inelasticXSec*0.1); // compute cell length in fm (1 mb = 0.1 fm^2)
      double latticeLength = 4.*r;                 // spread the lattice 2r in both (+/-) directions (total length = 4r)
      int ixmax = ceil(latticeLength/cellLength);
      ixmax*=2;                                    // number of cells in the x (and y) direction 
      double xmin = -latticeLength;                // minimal x (and y) on the lattice
      
      int nucALat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus A
      int nucBLat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus B
      int collLat[ixmax][ixmax];                   // lattice with the positions of the interactions -> for information only
      
      int countCollisionsA[100];
      int countCollisionsB[100];

      for (int i=0; i<100; i++)
	{
	  countCollisionsA[i]=0;
	  countCollisionsB[i]=0;
	}

      for (int i = 0; i < ixmax; i++)              // initialize cells: zero nucleons in every cell
	for (int j = 0; j < ixmax; j++)
	  {
	    nucALat[i][j] = 0;
	    nucBLat[i][j] = 0;
	    collLat[i][j] = 0;
	  }

      
      //cout <<"size=" << nucleusA.size() << endl;
      for (int i = 0; i < nucleusA.size(); i++)    // fill nucleons into cells
	{
	  posx = floor((nucleusA.at(i).x-xmin-b/2.)/cellLength);
	  posy = floor((nucleusA.at(i).y-xmin)/cellLength);
	  nucALat[posx][posy] += 1;
	  
	  posx = floor((nucleusB.at(i).x-xmin+b/2.)/cellLength);
	  posy = floor((nucleusB.at(i).y-xmin)/cellLength);
	  nucBLat[posx][posy] += 1;
	}
      
      int done;
      int numberOfpp;
      int eventNumber=0;
      numberOfpp=0;
      //cout << "jetXSec=" << jetXSec << endl;
      //cout << "inelasticXSec=" << inelasticXSec << endl;
      A=glauber->nucleusA();
      if (A==197.) //gold - Au 
	Z=79.; 
      else if (A==63.) //copper - Cu
	Z=29.;
      else if (A==208.) //lead - Pb
	Z=82.;
      else
	{
	  Z=0;
	}
      //cout << " Z=" << Z <<endl;
      //cout << " probability to generate jet event =" << jetXSec/inelasticXSec << endl;
      for (int ix = 0; ix < ixmax; ix++) 
	for (int iy = 0; iy < ixmax; iy++)
	  {
	    for (int i = 0; i < nucALat[ix][iy]; i++)
	      for (int j = 0; j < nucBLat[ix][iy]; j++)
		{
		  //cout << "ix+b/2=" << static_cast<int>(ix+b/2.) << " # there=" << nucALat[static_cast<int>(ix+b/2.)][iy] << endl;
		  if ( random->genrand64_real1() < jetXSec/inelasticXSec )
		    {
		      countCollisionsA[i]+=1;
		      countCollisionsB[j]+=1;
		      double xPositionInCell = (random->genrand64_real1())*cellLength;
		      double yPositionInCell = (random->genrand64_real1())*cellLength;
		      // see if colliding nucleons are neutrons
		      if ( random->genrand64_real1() > Z/A ) n1=1;
		      else n1=0;
		      if ( random->genrand64_real1() > Z/A ) n2=1;
		      else n2=0;
		      //cout << "n1=" << n1 << " , n2=" << n2 << endl;
		      pythia.next(n1,n2);                                                // generate event with pythia
		      done = 0;
		      //pythia.event.list();
		      for (int ip = 0; ip < pythia.event.size(); ++ip) 
			{
			  if (pythia.event[ip].isFinal())// && (pythia.event[i].id()<4 || pythia.event[i].id()==21))
			    // if the parton is final, i.e., present after the showers, then put it in the list
			    // currently heavy quarks are taken but will not lose energy! 
			    // only g and u,d,s lose energy while the 3 quarks are assumed to be massless
			    {
			      parton.id(pythia.event[ip].id());                     // set parton id
			      parton.status(pythia.event[ip].status());             // set parton status
			      parton.mass(pythia.event[ip].m());	                // set mass
			      if (fixedTemperature==0)
				{
				  //parton.x(xmin+(ix)*cellLength+cellLength/2.);       // set position
				  //parton.y(xmin+(iy)*cellLength+cellLength/2.);
				  parton.x(xmin+(ix)*cellLength+xPositionInCell);       // set position
				  parton.y(xmin+(iy)*cellLength+yPositionInCell);
				  parton.xini(xmin+(ix)*cellLength+xPositionInCell);       // set position
				  parton.yini(xmin+(iy)*cellLength+yPositionInCell);
				  parton.tini(hydroTau0);
				  if ( done == 0 )
				    {
				      //cout << "ix=" << ix << " iy=" << iy << " nucli=" << i 
				      //     << " nuclj=" << j << " ip=" << ip << " x=" << parton.x()<< " y=" << parton.y() << endl;
				      posx = floor((parton.x()-xmin)/cellLength);
				      posy = floor((parton.y()-xmin)/cellLength);
				      collLat[posx][posy] += 1;
				      
				      //cout << numberOfpp << endl;
				      numberOfpp++;
				      done = 1;
				    }
				}
			      else
				{
				  parton.x(0.);     
				  parton.y(0.);
				  parton.xini(0.);     
				  parton.yini(0.);
				  parton.tini(0.);
				}
			      parton.z(0.);
			      parton.zini(0.);
			      parton.col(pythia.event[ip].col());         	        // set color 
			      parton.acol(pythia.event[ip].acol());                 // set anti-color
			      parton.frozen(0);                                     // parton is not frozen (will evolve)
			      parton.p(pythia.event[ip].p());                       // set momentum
			      parton.pini(pythia.event[ip].p());                    // set initial momentum
			      parton.eventNumber(eventNumber);
			      parton.elasticCollisions(0);                          // set initial no. of el colls to zero
			      parton.source(0);                                     // All initial partons have itsSource=0
			      plist->push_back(parton);                             // add the parton to the main list
			    }
			} 
		      eventNumber++;
		    }
		}
	  }

//       for (int i=0; i<100; i++)
// 	{
// 	  cout << "collisions for parton i=" << i << " =" << countCollisionsA[i] << endl;
// 	  cout << "collisions for parton j=" << i << " =" << countCollisionsB[i] << endl;
// 	}

      //cout << "eventNumber = " << eventNumber << endl;
      //cout << "total number of NN =" << numberOfpp << endl;
      //cout << "eventNumber =" << eventNumber << endl;
      totalNNs = eventNumber;
      //output for plot
      ofstream fout1("./output/density.dat",ios::out); 
      ofstream fout2("./output/densityB.dat",ios::out); 
      ofstream fout3("./output/colls.dat",ios::out); 
      ofstream fout4("./output/NN.dat",ios::app); 
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
	    if ( j == ixmax-1 ) fout1 << endl;
	  }
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
	    if ( j == ixmax-1 ) fout2 << endl;
	  }
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
	    if ( j == ixmax-1 ) fout3 << endl;
	  }
      
      fout4 << numberOfpp << endl;
      
      fout1.close();
      fout2.close();
      fout3.close();
      fout4.close();
      //cout << " done generating initial hard collisions" << endl;
      return eventNumber;
    }
  else if (fullEvent==0) // sample only one hard collision per heavy-ion event
    {
      Parton parton;                                                // new Parton object
      double A;
      double Z;
      int n1;
      int n2;
      A=glauber->nucleusA();
      if (A==197.) //gold - Au 
	Z=79.; 
      else if (A==63.) //copper - Cu
	Z=29.;
      else if (A==208.) //lead - Pb
	Z=82.;
      else
	{
	  Z=0;
	}
      if ( random->genrand64_real1() > Z/A ) n1=1;
      else n1=0;
      if ( random->genrand64_real1() > Z/A ) n2=1;
      else n2=0;
      pythia.next(n1,n2);                                                // generate event with pythia
      //pythia.event.list(); 
      ReturnValue rv;
      //if (evolution && fixedTemperature==0) 
      if (!allFromCenter)
	{
	  if(nbinFromFile){
	    //Sample a random element of the list of number of collisions. -CFY 11/2/2010
	    double randomr = random->genrand64_real2();
	    int randomint = (int)(randomr*((double)Nbin));
	    rv.x = binary[randomint][0];
	    rv.y = binary[randomint][1];
	    cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
	  }
      
	  else if (glauberEnvelope)
	    rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
	  else
	    rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
	  //  cout << rv.x << " " << rv.y << endl;
	}
      else
	{
	  rv.x=0.;
	  rv.y=0.;
	}
      if (trackHistory)
	{
	  rv.x=initialXjet;
	  rv.y=initialYjet;
	}
      for (int i = 0; i < pythia.event.size(); ++i) 
	{
	  if (pythia.event[i].isFinal())// && (pythia.event[i].id()<4 || pythia.event[i].id()==21))
	    // if the parton is final, i.e., present after the showers, then put it in the list
	    // currently heavy quarks are taken but will not lose energy! 
	    // only g and u,d,s lose energy and the 3 quarks are assumed to be massless
	    {
	      parton.id(pythia.event[i].id());                      // set parton id
	      parton.status(pythia.event[i].status());              // set parton status
	      parton.mass(pythia.event[i].m());	                // set mass
	      
	      parton.x(rv.x);                                   // set position
	      parton.y(rv.y);
	      if (fixedTemperature==0)
		parton.tini(hydroTau0);
	      else
		parton.tini(0.);
	      parton.xini(rv.x);                                // set initial position to remember
	      parton.yini(rv.y);
	      
	      
	      /*
		if (evolution && fixedTemperature==0)
		{
		//parton.x(0.);                                   // set position
		//parton.y(0.);
		//parton.xini(0.);                                // set initial position to remember
		//parton.yini(0.);
		//pythia.event[i].p(20.,10.,30.,37.42);
		parton.x(rv.x);                                   // set position
		parton.y(rv.y);
		parton.xini(rv.x);                                // set initial position to remember
		parton.yini(rv.y);
		//cout << "x=" << rv.x << ", y=" << rv.y << endl;
		//cout << "p=" << pythia.event[i].p() << endl;
		
		}
		else
		{
		parton.x(0.);     
		parton.y(0.);
		}
	      */
	      parton.z(0.);
	      parton.zini(0.);
	      parton.col(pythia.event[i].col());         	        // set color 
	      parton.acol(pythia.event[i].acol());                  // set anti-color
	      parton.frozen(0);                                     // parton is not frozen (will evolve)
	      parton.p(pythia.event[i].p());                        // set momentum
	      parton.pini(pythia.event[i].p());                     // set momentum
	      parton.elasticCollisions(0);                          // set initial no. of el colls to zero
	      //parton.splits(0);                                   // set initial no. of splittings to zero
	      parton.source(0);                                     // All intial partons have itsSource=0
	      plist->push_back(parton);                             // add the parton to the main list
	      // check jettiness of di-quarks:
	      //if (abs(parton.id())>2000) cout << "id=" << parton.id() << " px=" << parton.p().px() 
	      //   	   << " py=" << parton.p().py() << " pz=" <<  parton.p().pz() << endl;
	    }
	} 
      return 1;
    }
}

int MARTINI::fragmentation(vector<Parton> ** plist, int currentEvent )
{
  if (fragmentationMethod == 1)
    {
      if (fullEvent == 0 )
	{
	  //pythia.event.list(); 
	  pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
	  int j = 0;
	  int id, col, acol, sid, scol, sacol, status;
	  int p_imax = plist[0]->size();
	  int p_jmax = plist[1]->size();
	  int numEvents = 0;
	  //int check = 0;
	  double mass, smass;
	  Vec4 pvec, spvec; 
	  // add all the higher momentum partons:
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons in the main list (those above momentum threshold)
	    {
	      id = plist[0]->at(p_i).id();
	      status = plist[0]->at(p_i).status();
	      if (status<1 || status>100) status = 1; // if status was not set, cure this here
	      //if (id==22){
	      //cout << "status=" << status << endl;
	      //check++;
	      //}
	      pvec = plist[0]->at(p_i).p();
	      col = plist[0]->at(p_i).col();
	      acol = plist[0]->at(p_i).acol();
	      mass = plist[0]->at(p_i).mass();
	      // append( id, status, col, acol, p, mass )
	      pythia.event.append(id,status,col,acol,pvec,mass); // copy existing parton into event record
	    }
	  //if(check>0) pythia.event.list(); 
	  if(pythia.forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest (if success add 1 to total events)
	  //if(check>0) pythia.event.list(); 
	  return numEvents;
	}
      else if (fullEvent == 1)
	{
	  //pythia.event.list(); 
	  pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
	  int j = 0;
	  int id, col, acol, sid, scol, sacol, status, eventNumber;
	  int p_imax = plist[0]->size();
	  int p_jmax = plist[1]->size();
	  int numEvents = 0;
	  //int check = 0;
	  double mass, smass;
	  Vec4 pvec, spvec; 
	  // add all the higher momentum partons:
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons in the main list (those above momentum threshold)
	    {
	      id = plist[0]->at(p_i).id();
	      status = plist[0]->at(p_i).status();
	      //if (status<1 || status>100) status = 1; // if status was not set, cure this here
	      if (status<1 || status>100) status = 62; // if status was not set, cure this here
	      //if (id==22){
	      //cout << "status=" << status << endl;
	      //check++;
	      //}
	      pvec = plist[0]->at(p_i).p();
	      col = plist[0]->at(p_i).col();
	      acol = plist[0]->at(p_i).acol();
	      mass = plist[0]->at(p_i).mass();
	      eventNumber = plist[0]->at(p_i).eventNumber();
	      // append( id, status, col, acol, p, mass )
	      if (eventNumber==currentEvent){
		if(id == 22){
		  cout << status << " " << pvec.px() << " " << pvec.py() << " " << pvec.pz() << endl;
		}
		pythia.event.append(id,status,col,acol,pvec,mass); // copy existing parton into event record
	      }
	    }
	  if(pythia.forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest (if success add 1 to total events)
	  //if(check>0) pythia.event.list(); 
	  //if(check>0) pythia.event.list(); 
	  return numEvents;
	}
    } 
  else if (fragmentationMethod == 2)
    {
      pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
      int j = 0;
      int i, is, k, l, ix, iy, iz, ixmax, izmax;
      int Ncells;
      int id, col, acol, sid, scol, sacol, status;
      int p_imax = plist[0]->size();
      int p_jmax = plist[1]->size();
      int numEvents = 0;
      //int check = 0;
      HydroInfo hydroInfo;
      double mass, smass, T, x, y, z, t, V;
      Vec4 pvec, spvec; 
      double Ngluons, Nquarks;
      int Igluons, Iquarks;
      int stringCounter, rni;
      double DX, DZ;
      double rn;
      Parton newOne, temp;
      int countStrange = 0;
      int countAStrange = 0;
      int countq = 0;
      int countaq = 0;
      int countQuarkJets;
      int countAntiQuarkJets;

      //cout << "doing fragmentation, method 2" << endl;
      if (hydroWhichHydro != 3 && hydroWhichHydro != 4) 
	{	
	  hydroZmax=hydroTauMax;
	  hydroDz=hydroDx;
	}
      
      DX = 1.; //fm
      DZ = 1.; //fm

      ixmax = floor((2.*hydroXmax)/DX+0.0001);
      izmax = floor((2.*hydroZmax)/DZ+0.0001);

      Ncells = ixmax*ixmax*izmax;
      vector<Parton> ** psublist;               // pointer to array of vector<Parton> objects
      vector<double> ** cvec;
      psublist = new vector<Parton> *[Ncells];  // pointer to array of vector<Parton> objects
      
      cvec = new vector<double> *[Ncells]; 
      
      for(i=0; i<Ncells; i++)
	{
	  psublist[i] = new vector<Parton>;    // psublist[i] is the ith sublist that holds the high momentum partons that were evolved
	  cvec[i] = new vector<double>;
	}
      
      countQuarkJets = 0;
      countAntiQuarkJets = 0;

      
      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons in the main list (those above momentum threshold) 
	{                                 // to sort them into sub-lists
 	  id = plist[0]->at(p_i).id();
	  
	  if (abs(id)>22) continue; // throw out weird things

	  plist[0]->at(p_i).status(1);
	  
	  x = plist[0]->at(p_i).x();                                // x value of position in [fm]
	  y = plist[0]->at(p_i).y();                                // y value of position in [fm]
	  z = plist[0]->at(p_i).z();                                // z value of position in [fm]
	  t = plist[0]->at(p_i).tFinal();
	  
	
	  ix = floor((hydroXmax+x)/DX+0.0001);                 // x-coordinate of the cell we are in now
	  iy = floor((hydroXmax+y)/DX+0.0001);                 // y-coordinate of the cell we are in now
	                                                            // note that x and y run from -hydroXmax to +hydroXmax
	                                                            // and ix and iy from 0 to 2*hydroXmax
	                                                            // hence the (hydroXmax+x or y) for both
	  iz = floor((hydroZmax+z)/DZ+0.0001);                 // z-coordinate of the cell we are in now
  
	
	  //cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;
	  i=ix+ixmax*(iy+iz*(ixmax)); //determine sublist number by position
	  //cout << "ix=" << ix << ", iy=" << iy << ", iz=" << iz << ", i=" << i << " out of " << Ncells << endl;
	  if (plist[0]->at(p_i).p().pAbs()<2.) continue; // dont add jet under 2 GeV
	  
	  if (abs(id)<21)
	    {
	      if (plist[0]->at(p_i).id()>0) countQuarkJets++;
	      if (plist[0]->at(p_i).id()<0) countAntiQuarkJets++;
	    }

	  if (t==0) continue;
	  
	  psublist[i]->push_back(plist[0]->at(p_i)); // add the parton to the correct sublist.
	  
	  x=ix*DX-hydroXmax+DX/2.;
	  y=iy*DX-hydroXmax+DX/2.;
	  z=iz*DZ-hydroXmax+DZ/2.;

	  if (fabs(z)>t-DZ) z=fabs(t-hydroDz)*z/fabs(z);

	  cvec[i]->push_back(x); // middle of the cell
	  cvec[i]->push_back(y);
	  cvec[i]->push_back(z);
	  cvec[i]->push_back(t);
	}

      stringCounter = 1; // gives the string a unique number

      for(i=0; i<Ncells; i++) //loop over all sublists (cells)
	{
	  if (psublist[i]->size() == 0) continue;
	  
	  //	  cout << "There is/are " << psublist[i]->size() << " hard parton(s) in list " << i << endl;
	  
	  x = cvec[i]->at(0);
	  y = cvec[i]->at(1);
	  z = cvec[i]->at(2);
	  t = cvec[i]->at(3);

	  //cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;

          //get the temperature and flow velocity at the current position and time:
	  if(fixedTemperature==25) //==0
	    {
	      hydroInfo = hydroSetup->getHydroValues(x, y, z, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
						     hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
	      
	      T = hydroInfo.T;
	    }
	  else T = fixedT;
	  
	  if (T==0.) T=0.1; // give it a minimum T
	  //cout << " T in this cell =" << T << " GeV" << endl;

	  V = DX*DX*DZ;

	  Ngluons = (V*16./(PI*PI)*T*T*T*1.202056903/pow(hbarc,3.)); // 1.202056903 is Riemann zeta(3).
	  Nquarks = (V*9./(2.*PI*PI)*T*T*T*1.202056903/pow(hbarc,3.));

	  //cout << "N_gluons in this cell = " << Ngluons << endl;
	  //cout << "N_quarks and anti-quarks in this cell = " << 6*Nquarks << endl;

	  Igluons = round(Ngluons);
	  Iquarks = round(6*Nquarks); // three flavors

	  for (k=0; k<Ngluons; k++)
	    {
	      newOne.id(21);                            // add a thermal gluon
 	      newOne.mass(0.); 
 	      newOne.frozen(1);
 	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
 	      newOne.p(random->thermal(T, -1));         // sample momentum from Bose distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.status(0);
 	      psublist[i]->push_back(newOne); 
	    }

	  for (k=0; k<Nquarks; k++)
	    {
	      // here decide which kind of quark / antiquark it is
	      rn=random->genrand64_real1();	      
	      if (rn<0.3333)
		{
		  id = 1;
		  mass = 0.33;
		}
	      else if (rn<0.6666)
		{
		  id = 2;
		  mass = 0.33;
		}
	      else 
		{
		  id = 3;
		  mass = 0.5;
		}
	      rn=random->genrand64_real1();	      
	      if (rn<0.5) // anti-quark
		id*=-1; 
	      
	      newOne.id(id); 
 	      newOne.mass(mass); 
 	      newOne.frozen(1);
 	      newOne.x(x);                              // set the new parton's initial position
 	      newOne.y(y);
 	      newOne.z(z);
 	      newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
	      newOne.pini(newOne.p());
	      newOne.elasticCollisions(0);              // set initial no. of el colls to zero
	      newOne.status(0);
	      psublist[i]->push_back(newOne); 
	    }
	  
	  
	  //  shuffle order (so we don't always start with gluons but treat every parton on equal grounds)
	  for (is=0; is<psublist[i]->size(); is++) 
	    {
	      rni = static_cast<int>(random->genrand64_real1()*10000) % psublist[i]->size();  // generate a random position
	      temp = psublist[i]->at(is);                   // switch parton at position 'is' with that at 'rni' ...
	      psublist[i]->at(is) = psublist[i]->at(rni); 
	      psublist[i]->at(rni) = temp;
	    }
		 
	  
// 	  is = 0; 
// 	  while (psublist[i]->begin()+is<psublist[i]->end())  // remove heavy quarks (keep photons)
// 	    {
// 	      if (psublist[i]->at(is).id()>22) 
// 		{
// 		  psublist[i]->erase(psublist[i]->begin()+is);
// 		  is--;
// 		}
// 	      is++;
// 	    }

	  for (k=0; k<psublist[i]->size(); k++)  // loop over all partons in the sublist and reset color links
	    {
	      psublist[i]->at(k).col(0);
	      psublist[i]->at(k).acol(0);
	      //cout << "id=" << psublist[i]->at(k).id() << endl;
	    }

	  stringCounter++;
	  
	  for (int tries = 0 ; tries<4; tries++) // try to connect among thermal partons (only if there aren't enough add one later)
	    {
	      for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist
		{
		  //jets have status 1, medium partons have status 0
		  if (psublist[i]->at(k).status()==1) //  it is a jet
		    {
		      if (psublist[i]->at(k).id()==21) //   it is a gluon
			{
			  if (psublist[i]->at(k).acol()!=0 && psublist[i]->at(k).col()!=0) continue; // is already connected
			  
			  for (l=0; l<psublist[i]->size(); l++)
			    {
			      if (l==k) continue; // do not connect to itself
			      if (psublist[i]->at(l).id()!=21 && psublist[i]->at(l).id()!=22) // if quark or anti-quark
				{
				  if (psublist[i]->at(l).acol()!=0 || psublist[i]->at(l).col()!=0) continue; // is already connected
				  
				  if (psublist[i]->at(l).id()<0 && psublist[i]->at(k).col()==0) //if anti-quark
				    {
				      psublist[i]->at(l).acol(stringCounter); 
				      psublist[i]->at(k).col(stringCounter);  // connect to gluon's color
				      stringCounter++;
				      // cout << " 1 doing jg and a-q" << endl;
				    }
				  else if(psublist[i]->at(l).id()>0 && psublist[i]->at(k).acol()==0)// if quark
				    {
				      psublist[i]->at(l).col(stringCounter);
				      psublist[i]->at(k).acol(stringCounter); // connect to gluon's anti-color
				      stringCounter++;
				      //cout << " 2 doing jg and q" << endl;
				    }
				}
			      else // if gluon 
				{
				  if (psublist[i]->at(l).acol()!=0 && psublist[i]->at(l).col()!=0) continue; // is already connected
				  if(psublist[i]->at(k).status()>=0) // for jets and thermal (unnecessary test)
				    {
				      rn=random->genrand64_real1();	
				      if (rn<0.3333) // connect color end
					{
					  if (psublist[i]->at(l).acol()==0 && psublist[i]->at(k).col()==0)
					    {
					      psublist[i]->at(k).col(stringCounter);
					      psublist[i]->at(l).acol(stringCounter);
					      stringCounter++;
					      //cout << " 3 doing jg and g color end" << endl;
					    }
					}
				      else if (rn<0.6666) // connect anti-color end
					{
					  if (psublist[i]->at(l).col()==0 && psublist[i]->at(k).acol()==0)
					    {
					      psublist[i]->at(k).acol(stringCounter); 
					      psublist[i]->at(l).col(stringCounter);
					      stringCounter++;
					      //cout << " 4 doing jg and g anti-color end" << endl;
					    }
					}
				      else //connect both ends if possible
					{
					  if (psublist[i]->at(l).acol()==0 && psublist[i]->at(k).col()==0)
					    {
					      psublist[i]->at(k).col(stringCounter);
					      psublist[i]->at(l).acol(stringCounter);
					      stringCounter++;
					      //cout << " 5 doing jg and g color end when both" << endl;
					    }
					  if (psublist[i]->at(l).col()==0 && psublist[i]->at(k).acol()==0)
					    {
					      psublist[i]->at(k).acol(stringCounter); 
					      psublist[i]->at(l).col(stringCounter);
					      stringCounter++;
					      //cout << " 6 doing jg and g anti-color end when both" << endl;
					    }
					}
				    }
				}
			    }
			}
		      else // it is a jet quark or anti-quark
			{
			  if (psublist[i]->at(k).acol()!=0 || psublist[i]->at(k).col()!=0) continue; // is already connected
			  for (l=0; l<psublist[i]->size(); l++)
			    {
			      if (l==k) continue; // do not connect to itself
			      if (psublist[i]->at(l).id()!=21 && psublist[i]->at(l).id()!=22) // if quark or anti-quark
				{
				  if (psublist[i]->at(l).acol()!=0 || psublist[i]->at(l).col()!=0) continue; // is already connected
				  if (psublist[i]->at(k).id()>0 && psublist[i]->at(l).id()<0 && psublist[i]->at(k).col()==0)
				    //if quark and thermal anti-quark
				    {
				      psublist[i]->at(l).acol(stringCounter); 
				      psublist[i]->at(k).col(stringCounter);  
				      stringCounter++;
				      //cout << " 7 doing jq and a-q" << endl;
				    }
				  else if (psublist[i]->at(k).id()<0 && psublist[i]->at(l).id()>0 && psublist[i]->at(k).acol()==0) 
				    //if anti-quark and thermal quark
				    {
				      psublist[i]->at(l).col(stringCounter);
				      psublist[i]->at(k).acol(stringCounter);
				      stringCounter++;
				      //cout << " 8 doing ja-q and q" << endl;
				    }
				}
			      else // if thermal gluon
				{
				  if (psublist[i]->at(k).id()>0 && psublist[i]->at(l).acol()!=0) continue;// is already connected
				  if (psublist[i]->at(k).id()<0 && psublist[i]->at(l).col()!=0) continue; // is already connected
				  
				  if (psublist[i]->at(k).id()>0 && psublist[i]->at(k).col()==0) // if jet quark
				    {
				      psublist[i]->at(k).col(stringCounter);
				      psublist[i]->at(l).acol(stringCounter);
				      stringCounter++;
				      //cout << " 9 doing jq and g" << endl;
				    }
				  else if (psublist[i]->at(k).id()<0 && psublist[i]->at(k).acol()==0) // if jet anti-quark
				    {
				      psublist[i]->at(k).acol(stringCounter);
				      psublist[i]->at(l).col(stringCounter);
				      //cout << " 10 doing ja-q and g" << endl;
				      stringCounter++;
				    }
				}
			    }
			}
		    }
		} // end loop over all partons in the sublist
	    } // end loop over tries

	  // finally connect open ends of gluons - careful here

	  for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist
	    {
	      if (psublist[i]->at(k).id()==21 && ( (psublist[i]->at(k).col()!=0 && psublist[i]->at(k).acol()==0)
						   ||
						   (psublist[i]->at(k).acol()!=0 && psublist[i]->at(k).col()==0) )
		  )
		{
		  if (psublist[i]->at(k).status()==1)
		    //cout << " gluon jet needs extra string attached!" << endl;
		  rn=random->genrand64_real1();	     // add a thermal quark or anti-quark
		  if (rn<0.4)
		    {
		      id = 1;
		      mass = 0.33;
		    }
		  else if (rn<0.8)
		    {
		      id = 2;
		      mass = 0.33;
		    }
		  else 
		    {
		      id = 3;
		      mass = 0.5;
		    }
		  
		  newOne.mass(mass); 
		  newOne.frozen(1);
		  newOne.x(x);                              // set the new parton's initial position
		  newOne.y(y);
		  newOne.z(z);
		  newOne.p(random->thermal(T, 1)); //T         // sample momentum from Fermi distribution at temperature T
		  newOne.elasticCollisions(0);              // set initial no. of el colls to zero
		  newOne.pini(newOne.p());
		  newOne.status(0);
		  
		  if (psublist[i]->at(k).acol()==0) // add a thermal quark and connect to it
		    {
		      //cout << "adding id " << id << " parton" << endl;
		      newOne.id(id); 
		      newOne.col(stringCounter);
		      newOne.acol(0);
		      psublist[i]->at(k).acol(stringCounter);
		      stringCounter++;
		      psublist[i]->push_back(newOne); 
		    }

		  if (psublist[i]->at(k).col()==0) // add a thermal anti-quark and connect to it
		    {
		      //cout << "adding id " << -id << " parton" << endl;
		      newOne.id(-id); 
		      newOne.acol(stringCounter);
		      newOne.col(0);
		      psublist[i]->at(k).col(stringCounter);
		      stringCounter++;
		      psublist[i]->push_back(newOne); 
		    }
 		}
	    }
	  
	  is = 0;
	  while (psublist[i]->begin()+is<psublist[i]->end())  // remove unconnected thermal partons
	    {
	      if (psublist[i]->at(is).id()!=22 && psublist[i]->at(is).col()==0 && psublist[i]->at(is).acol()==0)
		{
		  //cout << "erasing position " << is << " with id=" << psublist[i]->at(is).id() << " and col/acol=" 
		  //   << psublist[i]->at(is).col() << "/" << psublist[i]->at(is).acol() << endl;
		  psublist[i]->erase(psublist[i]->begin()+is);
		  is--;
		}
	      is++;
	    }
	  

	//   for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist and print content
// 	    {
// 	      if (psublist[i]->at(k).status() == 1)
// 		cout << " cell " << i << ": jet " << k << " has ID " << psublist[i]->at(k).id() 
// 		     << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
// 	      else if (psublist[i]->at(k).status() == 0)
// 		cout << " cell " << i << ": thermal parton " << k << " has ID " << psublist[i]->at(k).id() 
// 		     << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
// 	      else
// 		cout << " cell " << i << ": weird parton " << k << " has ID " << psublist[i]->at(k).id() 
// 		     << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
// 	    }


	  for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist and append to pythia list
	    {
	      status = psublist[i]->at(k).status();
	      if (status<1 || status>100) status = 1; // if status was not set, cure this here
	      //if (id==22){
	      //cout << "status=" << status << endl;
	      //check++;
	      //}
	      id = psublist[i]->at(k).id();
	      if (abs(id)<4)
		{
		  if (id<0) countaq++;
		  if (id>0) countq++;
		  if (id==3) countStrange++;
		  if (id==-3) countAStrange++;
		}
	      pvec = psublist[i]->at(k).p();
	      col = psublist[i]->at(k).col();
	      acol = psublist[i]->at(k).acol();
	      mass = psublist[i]->at(k).mass();
	      //cout << k << " ID=" << psublist[i]->at(k).id() 
	      //   << " col=" << col << " acol=" << acol << " mass=" << mass << " p=" << pvec << endl;
	      pythia.event.append(id,status,col,acol,pvec,mass); // copy existing parton into event record
	    }
	} // end loop over all sublists


      //if(check>0) pythia.event.list(); 
      if(pythia.forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest (if success add 1 to total events)
      //if(check>0) pythia.event.list(); 
      //cout << "q-jets " << countQuarkJets << endl;
      //cout << "qbar-jets " << countAntiQuarkJets << endl;
      //cout << "s=" << countStrange << endl;
      //cout << "sbar=" << countAStrange << endl;
      //cout << "q's=" << countq << endl;
      //cout << "qbar's=" << countaq << endl;
      

      for(i=0; i<Ncells; i++)
	{
	  delete psublist[i];    // psublist[i] is the ith sublist that holds the high momentum partons that were evolved
	  delete cvec[i];
	}
      delete psublist;
      delete cvec;
     
      return numEvents;
    }
}

bool MARTINI::readFile(string fileName, bool warn) 
{
  // Open file with updates.
  const char* cstring = fileName.c_str();
  ifstream is(cstring);  
  if (!is) {
    info.errorMsg("Error in MARTINI::readFile: did not find file", fileName);
    return false;
  }
  // Read in one line at a time.
  string line;
  bool accepted = true;
  while ( getline(is, line) ) 
    {
      // Process the line.
      if ( !readString( line, false ) ) // if it isn't found do not warn but pass on to PYTHIA - PYTHIA can warn if wanted
	if (!pythia.readString( line, warn ))
	  accepted = false;
      // Reached end of input file.
    };
  return accepted;
}

bool MARTINI::readString(string line, bool warn) 
{
  // If empty line then done.
  if (line.find_first_not_of(" ") == string::npos) return true;

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isalnum(line[firstChar])) return true; 

  // Send on to Settings.
  return settings.readString(line, warn);
}

void MARTINI::initPythia()
{
  // let PYTHIA know my desired infrared cutoff on the jet momentum for the calculation of the jet cross section
  if (jetpTmin==0.) jetpTmin=0.2; // just make sure there is no divergence
  pythia.setJetpTmin(jetpTmin);
  
  // display the modified PYTHIA settings:
  pythia.settings.listChanged();
  
  // initialize PYTHIA
  pythia.init( 2212, 2212, cmEnergy ); // 2 protons, (if 14000 each 7000 GeV) /2212 is proton
  cout << endl << "PYTHIA initialized by MARTINI with sqrt(s)=" << cmEnergy << " GeV." << endl;
 
  // set the relevant cross sections that PYTHIA was so kind to compute
  elasticXSec = pythia.elasticCrossSection();
  totalXSec = pythia.totalCrossSection();
  jetXSec = pythia.jetCrossSection();
  inelasticXSec = totalXSec-elasticXSec;

  // output the cross sections
  cout << endl << "[MARTINI::initPythia]:" << endl;
  cout << "Jet cross section = " << jetXSec << " mb, with a pTmin = " << jetpTmin << " GeV." << endl;
  cout << "Total cross section = " << totalXSec << " mb." << endl;
  cout << "Total inelastic cross section = " << inelasticXSec << " mb." << endl;
  
}

bool MARTINI::init(int argc, char** argv) 
{
  /// set all parameters as they were initialized
  string sseed;
  int seed = 0;
  if (argc>1)
    {
      sseed = argv[1];
      seed = atoi(sseed.c_str()); 
    }
  // For FIC:
  file_number = seed%5 +1;
  cout << "file_number = " << file_number << endl;
  seed *= 10000;
  fstream foutt("test.dat",ios::out); 
  foutt << "testing .. " << endl;
  pCut = settings.parm("General:pCut");
  fixedT = settings.parm("General:Temperature");
  maxTime = settings.parm("General:MaxTime");
  dtfm = settings.parm("General:TimeStep");
  alpha_s = settings.parm("General:AlphaS");
  runs = settings.mode("General:Events");
  Nf = settings.mode("General:NumberOfFlavors");
  jetpTmin = settings.parm("General:JetPTMin");
  fullEvent = settings.mode("General:FullEvent");
  moveBeforeTau0 = settings.mode("General:MoveBeforeTau0");
  cmEnergy = settings.parm("General:cmEnergy");
  fullVacuumShower = settings.flag("General:FullVacuumShower");
  Ldependence = settings.mode("General:Ldependence");
  fragmentationMethod = settings.mode("General:FragmentationMethod");
  rateSelector = settings.mode("General:RadiativeRateSet");
  multElasticRate = settings.parm("General:MultiplyElasticRateBy");
  examineHQ = settings.mode("General:examineHQ");
  cout << "examineHQ = " << examineHQ << endl;
  setInitialPToZero = settings.mode("General:setInitialPToZero");
  cout << "setInitialPToZero = " << setInitialPToZero << endl;
  charmWidth = settings.parm("General:charmWidth");
  bottomWidth = settings.parm("General:bottomWidth");
  T_C_HQ = settings.parm("General:T_C_HQ");
  TwoPiTD_HQ = settings.parm("General:TwoPiTD_HQ");
  totalHQXSec = settings.parm("General:totalHQXSec");

  if (rateSelector == 2 && alpha_s!=0.3)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.3 at the moment. Exiting." << endl;
      exit(1);
    }
  if (rateSelector == 3 && alpha_s!=0.36)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.36 at the moment. Exiting." << endl;
      exit(1);
    }
  if (rateSelector == 4 && alpha_s!=0.35)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.35 at the moment. Exiting." << endl;
      exit(1);
    }
  if (rateSelector == 5 && alpha_s!=0.3)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.3 at the moment. Exiting." << endl;
      exit(1);
    }
  if (rateSelector == 6 && alpha_s!=0.4)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.4 at the moment. Exiting." << endl;
      exit(1);
    }
  if (rateSelector == 7 && alpha_s!=0.38)
    {
      cout << "The chosen radiative rates with cut are only available for alpha_s=0.38 at the moment. Exiting." << endl;
      exit(1);
    }
  
  if (settings.flag("General:RadiativeProcesses")) doRadiative=1;
  else doRadiative=0;
  if (settings.flag("General:ElasticCollisions")) doElastic=1;
  else doElastic=0;
  if (settings.flag("General:TransferTransverseMomentum")) transferTransverseMomentum=1;
  else transferTransverseMomentum=0;
  if (settings.flag("General:Evolution")) evolution=1;
  else evolution=0;
  if (settings.flag("General:FixedEnergy")) fixedEnergy=1;
  else fixedEnergy=0;
  if (settings.flag("General:FixedTemperature")) fixedTemperature=1;
  else fixedTemperature=0;
  if (settings.flag("General:Fragmentation")) fragmentationSwitch=1;
  else fragmentationSwitch=0;
  if (settings.flag("General:PhotonProduction")) photonSwitch=1;
  else photonSwitch=0;
  if (settings.flag("PDF:NuclearEffects")) nuclearEffects=1;
  else nuclearEffects=0;
  useLHAPDF = settings.flag("PDF:useLHAPDF");
  trackPartons = settings.flag("General:TrackPartons"); 
  trackHistory = settings.flag("General:TrackHistory"); 
  initialXjet = settings.parm("General:InitialXjet");
  initialYjet = settings.parm("General:InitialYjet");
  initialPXjet = settings.parm("General:InitialPXjet");
  initialPYjet = settings.parm("General:InitialPYjet");
  allFromCenter = settings.flag("General:allFromCenter");
  pCutPropToT = settings.flag("General:pCutPropToT");

  hydroTau0 = settings.parm("Hydro:tau0");
  hydroTauMax = settings.parm("Hydro:taumax");
  hydroDtau = settings.parm("Hydro:dtau");
  hydroDx = settings.parm("Hydro:dx");
  hydroDz = settings.parm("Hydro:dz");
  hydroXmax = settings.parm("Hydro:xmax");
  hydroZmax = settings.parm("Hydro:zmax");
  hydroTfinal = settings.parm("Hydro:Tfinal");
  hydroWhichHydro = settings.mode("Hydro:WhichHydro");  
  hydroSubset = settings.mode("Hydro:Subset");
  hydroViscous = settings.flag("Hydro:Viscous");  
  fixedDistribution = settings.flag("Hydro:fixedDistribution");
 
  glauberTarget = settings.word("Glauber:Target");
  glauberProjectile = settings.word("Glauber:Projectile");
  glauberImpactParam = settings.parm("Glauber:b");
  glauberIMax = settings.mode("Glauber:MaxInterpolationPoints");
  glauberEnvelope = settings.flag("Glauber:Envelope");
  
  PDFname = settings.word("PDF:LHAPDFset");
  PDFmember = settings.mode("PDF:LHAPDFmember");
  quarksEqualGluons = settings.mode("General:quarksEqualGluons");

  //Reading of these parameters added. -CFY 11/2/2010
  cout << "Is any of this being called?" << endl;
  nbinFromFile = settings.mode("General:Nbin_from_File");
  cout << "nbinFromFile = " << nbinFromFile << endl;
  //Either the name of the evolution file, or a beginning tag for the 
  //file for the case nbinFromFile == 1. -CFY
  evolution_name = settings.word("General:evolution_name");
  //The number of collisions in the file, counted up from zero
  cout << "OK up to here?" << endl;
  cout << "evolution_name = " << evolution_name << endl;
  Nbin = 0;

  stringstream s;
  s << argv[1];
  
  //If nbinFromFile = 1, read this list in.
  if(nbinFromFile == 1){
    string file_name = "x_y_binary_"+evolution_name+"_"+s.str()+".dat";
    ifstream binfile(file_name.c_str());
    ifstream justcounting(file_name.c_str());
    //const char* nbstring = nbinName.c_str();
    //ifstream binfile(nbstring);
    //Count the number of lines in the file of collision coordinates.
    //ifstream justcounting(nbstring);
    double ddummy;
    while(!justcounting.eof()){
      justcounting >> ddummy;
      justcounting >> ddummy;
      if(!justcounting.eof())Nbin++;
    }
    justcounting.close();
    
    binary = new double*[Nbin];
    double dvalue;
    for(int ib=0; ib<Nbin; ib++){
      binary[ib] = new double[2];
      binfile >> dvalue;
      binary[ib][0] = dvalue;
      binfile >> dvalue;
      binary[ib][1] = dvalue;
    }
    cout << "The number of collisions in this event = " << Nbin << endl;
    cout << "The coordinates of the last collision are " << binary[Nbin-1][0] << ", " << binary[Nbin-1][1] << endl;
  }

  // remove old NuTInST.dat and NuPInSP.dat 
  // in case things changed (only takes a second to generate the new ones
  string AinFile; // variable for input value
  string path     = "";
  string subfolder = "/main_HQ/data/";
  string file;
  int ik;
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
      path += subfolder;
    }
  else path = ".."+subfolder;
  
  file = path + "NuTInST.dat";
  cout << "reading target data from" << file.c_str() << endl;
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin)
    {
      cerr << "This file does not yet exist: " << file << endl;
    }
  else
    {
      ik = 0;
      while ( !fin.eof() && ik<2 )
	{
	  fin >> AinFile;
	  ik++;
	}
      fin.close();
      cout << "Target in file: " << AinFile << ", Target to do: " << glauberTarget << " - ";
      if(AinFile.compare(glauberTarget) != 0)
	{
	  remove( file.c_str() ); 
	  cout << "old NuTInST.dat removed" << endl;
	}
      else
	{
	  cout << "no need to generate new files." << endl;
	}
    }

  file = path + "NuPInSP.dat";
  cout << "reading projectile data from" << file.c_str() << endl;
  
  fin.open(file.c_str());
  if(!fin)
    {
      cerr << "This file does not yet exist: " << file << endl;
    }
  else
    {
      ik = 0;
      while ( !fin.eof() && ik<2 )
	{
	  fin >> AinFile;
	  ik++;
	}
      fin.close();
      
      cout << "Projectile in file: " << AinFile << ", Projectile to do: " << glauberProjectile << " - ";
      if(AinFile.compare(glauberTarget) != 0)
	{
	  remove( file.c_str() );
	  cout << "old NuPInSP.dat removed" << endl;
	}
      else
	{
	  cout << "no need to generate new files." << endl;
	}
    }

  ///read in data files with transition rates
  import->init(rateSelector);
  rates   = new Rates(rateSelector);

  //importLRates->init();
  
  //cout << "rate=" << importLRates->getRate(10., 5., 2.) << endl;

  // output of the rates into a file:
  //import->show_dGamma();

  /// initialize random number generator with current time as seed
  long long rnum;

  rnum=time(0)+seed;
  
  foutt << "random seed used=" << rnum << " made from time " << rnum -seed << " and argument " << seed << endl;
  random->init_genrand64(rnum);

  // PYTHIA settings:
  pythia.readString("Random:setSeed = on");

  // set random seed as current time
  srand(rnum);
  int now=rand()/100;
  stringstream nowstr;
  nowstr << "Random:seed = " << now;
  string timestring = nowstr.str();

  pythia.readString(timestring);    

  // turn on all hard processes
  if(examineHQ == 0){
    pythia.readString("HardQCD:all = on");
  }
  else{
    pythia.readString("HardQCD:gg2ccbar = on");
    pythia.readString("HardQCD:qqbar2ccbar = on");
    pythia.readString("HardQCD:gg2bbbar = on");
    pythia.readString("HardQCD:qqbar2bbbar = on");
  }

  // Let PYTHIA's timelike shower know about the medium (for the correct p_T/virtuality cutoff)
  if (evolution && fullVacuumShower==0)
    pythia.readString("TimeShower:doMedium = on");

  stringstream tauHydrostr;
  tauHydrostr << "TimeShower:tauHydro = " << hydroTau0;
  string tHstring = tauHydrostr.str();
  pythia.readString(tHstring);

  // finish after doing hard process
  pythia.readString("PartonLevel:all = on"); // off to only get hard process partons
  pythia.readString("HadronLevel:all = off");

  //LHAPDF::setlhaparm("EPS08");

  if (evolution)
    {
      // stop shower at scale p_T given here:
      //pythia.readString("SpaceShower:pTmin = 0.33"); // default = 0.2 - 0.33 from 1/0.6fm *0.197 GeV fm
      //pythia.readString("SpaceShower:pTminChgQ = 0.33"); // default = 0.2 - 0.33 from 1/0.6fm *0.197 GeV fm
      //pythia.readString("TimeShower:pTmin = 0.5"); // default = 0.5
    }

  if ( fullEvent == 1) // use the jetPTMin only in a full event.
    {
      stringstream ptminstr;
      ptminstr << "PhaseSpace:pTHatMin =  " << jetpTmin;
      string ptminstring = ptminstr.str();
      pythia.readString(ptminstring);
    }
  
  pythia.readString("Check:event = off");

  // set the PDF and possible nuclear effects. note: if nuclear effects are chosen, the use of LHAPDF is enforced!
  if ( nuclearEffects!=0 )
    {    
      if (!useLHAPDF) cout << "[MARTINI]:WARNING: using nuclear effects - turned on LHAPDF against initial settings." << endl;
      // init glauber to get the atomic number stored in glauber->nucleusA 
      glauber->preInit(inelasticXSec,glauberTarget,glauberProjectile,glauberImpactParam,glauberIMax);
      pythia.readString("PDF:useLHAPDF = 1");
      pythia.readString("PDF:nuclearEffects = 1");
      //pythia.readString("PDF:LHAPDFset = cteq5l.LHgrid");
      //pythia.readString("PDF:LHAPDFset = H12000ms.LHgrid");
      stringstream Astr;
      Astr << "PDF:atomicNumber = " << glauber->nucleusA();
      string Astring = Astr.str();
      pythia.readString(Astring);
      stringstream PDFnameStr;
      PDFnameStr << "PDF:LHAPDFset = " << PDFname;
      string PDFnameString = PDFnameStr.str();
      pythia.readString(PDFnameString);
      stringstream PDFmemberStr;
      PDFmemberStr << "PDF:LHAPDFmember = " << PDFmember;
      string PDFmemberString = PDFmemberStr.str();
      pythia.readString(PDFmemberString);
    }
  else
    {
      if(useLHAPDF) 
	{
	  pythia.readString("PDF:useLHAPDF = 1");
	  stringstream PDFnameStr;
	  PDFnameStr << "PDF:LHAPDFset = " << PDFname;
	  string PDFnameString = PDFnameStr.str();
	  pythia.readString(PDFnameString);
	  stringstream PDFmemberStr;
	  PDFmemberStr << "PDF:LHAPDFmember = " << PDFmember;
	  string PDFmemberString = PDFmemberStr.str();
	  pythia.readString(PDFmemberString);
	}
      pythia.readString("PDF:nuclearEffects = 0");
    }
 
  lattice = new vector<HydroCell>; // create new lattice object - list of cells with hydro information

  // read the hydro data from the data file(s)
  if (evolution && fixedTemperature==0) 
    {
      if(nbinFromFile == 1){
	hydroSetup->readHydroData(glauber->nucleusA(),hydroTau0, hydroTauMax, hydroDtau, 
				  hydroXmax, hydroZmax, hydroDx, hydroDz, 
				  hydroWhichHydro, file_number, hydroSubset, hydroViscous, glauberImpactParam, lattice, evolution_name); 
      }
      else{
	hydroSetup->readHydroData(glauber->nucleusA(),hydroTau0, hydroTauMax, hydroDtau, 
				  hydroXmax, hydroZmax, hydroDx, hydroDz, 
				  hydroWhichHydro, hydroSubset, hydroViscous, glauberImpactParam, lattice, evolution_name); 
      }
    }

  // init PYTHIA with given center of mass energy. will be changed for full AA collision.
  initPythia();

  //do it again 
  glauber->init(inelasticXSec,glauberTarget,glauberProjectile,glauberImpactParam,glauberIMax,glauberEnvelope);

  //testing->sampleTA(random, glauber);
  foutt.close();

  if (trackPartons)
    {
      ofstream fouttracks("./output/qtracks.dat",ios::out); 
      ofstream fouttracksg("./output/gtracks.dat",ios::out); 
      fouttracks <<  " t  " << " x " << " y " << " z " << " deltaE " << endl;
      fouttracks.close();
      fouttracksg <<  " t  " << " x " << " y " << " z " << " deltaE " << endl;
      fouttracksg.close();
      
      ofstream fouttracks2("./output/qtracks2.dat",ios::out); 
      ofstream fouttracksg2("./output/gtracks2.dat",ios::out); 
      fouttracks2 <<  " t  " << " x " << " y " << " z " << " deltaE " << endl;
      fouttracks2.close();
      fouttracksg2 <<  " t  " << " x " << " y " << " z " << " deltaE " << endl;
      fouttracksg2.close();

      ofstream foutpovray("./output/povraytracks.dat",ios::out);
      foutpovray <<  "  "  << endl;
      foutpovray.close();

      ofstream foutpovrayg("./output/povraytracksg.dat",ios::out);
      foutpovrayg <<  "  "  << endl;
      foutpovrayg.close();
    }

  if (trackHistory)
    {
      Tt = new vector<double>; // stores the history of the trajectory: t, x(t), y(t), z(t), dE/dt(t), dpx/dt(t), dpy/dt(t), dpz/dt(t) *for one particle*
      Tx = new vector<double>;
      Ty = new vector<double>;
      Tz = new vector<double>;
      TE = new vector<double>;
      TQGPfrac = new vector<double>;
      Tpx = new vector<double>;
      Tpy = new vector<double>;
      TdEdt = new vector<double>;
      Tdpxdt = new vector<double>;
      Tdpydt = new vector<double>;
      Tdpzdt = new vector<double>;
    }

  //One final step for easy automation of MARTINI:
  //hydroTauMax is reset for the case where writing to evolution.dat ended early
  //(due to all cells freezing out):
  hydroTauMax = hydroTau0 + hydroDtau*(int)((double)lattice->size()/( (hydroXmax/hydroDx)*(hydroXmax/hydroDx)*2.*(hydroZmax/hydroDz) ) ) ;
  cout << "hydroTauMax = " << hydroTauMax << endl;

  //Let's test getHydroValues:
  if(evolution == 1 && fixedTemperature == 0){
    double zrange = 5.;
    //double xrange = 5.;
    //double x, y, z, t;
    int numpoints = 20;
    HydroInfo hydroInfo;
    
    cout << "hydroTau0 = " << hydroTau0 << endl;
    cout << "hydroTauMax = " << hydroTauMax << endl;
    
    //First, the temperature at the center of the collision at tau = hydroTau0:
    hydroInfo = hydroSetup->getHydroValues(0., 0., 0., hydroTau0, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
					   hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
    cout << "Temperature at (0,0,0,0) = " << hydroInfo.T << endl;

    for(int inum=0; inum<numpoints; inum++){
      double hydroZ = (double)inum*zrange/(double)numpoints;
      double hydroTau = sqrt((hydroTauMax-hydroTau0)*(hydroTauMax-hydroTau0)-hydroZ*hydroZ );
      if(hydroTau > hydroTau0 && hydroTau < hydroTauMax){
	hydroInfo = hydroSetup->getHydroValues(0., 0., hydroZ, hydroTauMax-hydroTau0, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
					       hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
	cout << (double)inum*(zrange/(double)numpoints) << " " << hydroInfo.T << " " << hydroInfo.vx
	     << " " << hydroInfo.vz << endl;
      }
    }
    //for(int inum=0; inum<numpoints; inum++){
    //hydroInfo = hydroSetup->getHydroValues(double)inum*(xrange/(double)numpoints), 0., 0., 1., hydroXmax, hydroZmax, 
    //	       hydroTauMax, hydroTau0, hydroDx, hydroDz, hydroDtau, hydroWhichHydro, 
    //	       fixedDistribution, lattice, trackHistory);
    //cout << (double)inum*(xrange/(double)numpoints) << " " << hydroInfo.T << " " << hydroInfo.vx
    // << " " << hydroInfo.vz << endl;
    //}
    
    //for(int inum=0; inum<numpoints; inum++){
    //hydroInfo = hydroSetup->getHydroValues(0., 0., 0., 1.+(double)inum*(taurange/(double)numpoints), hydroXmax, hydroZmax, 
    //		       hydroTauMax, hydroTau0, hydroDx, hydroDz, hydroDtau, hydroWhichHydro, 
    //		       fixedDistribution, lattice, trackHistory);
    //cout << 1.+(double)inum*(taurange/(double)numpoints) << " " << hydroInfo.T << " " << hydroInfo.vx
    // << " " << hydroInfo.vz << endl;
    //}
  }

  return true;
}

vector<Parton>* MARTINI::next()
{
  int mt;
  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  plist = new vector<Parton> *[2];  
  plist[0] = new vector<Parton>;    // plist[0] is the main list that holds the high momentum partons that are evolved

  int counter = 0;                  // reset counter in the beginning of every event
  
  plist[0]->clear();                // clear the parton list - not necessary 
  
  generateEvent(plist[0]);
  
  if (evolution == 1)             // evolve in medium if settings allow it
    {
      mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
      for(int i=0; i<mt; i++)       // loop over all time steps 
	{
	  counter = evolve(plist, counter, i);
	  counter+=1;
	}
    }
  return plist[0];
}

// --------------------------------------------------------------------------------------------------------------------------------------

// this has been moved to main.cpp !!
// void MARTINI::pythiaEvents()
// {
//   Parton jp1;                       // used for fixed energy runs
//   Parton jp2;
  
//   vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
//   plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
//   plist[0] = new vector<Parton>;    // plist[0] is the main list that holds the high momentum partons that are evolved
//   plist[1] = new vector<Parton>;    // plist[1] is auxiliary list holding low mom. partons to be matched when fragmenting
//                                     // plist[1] is not used currently 02/24/2009 
//   Vec4 vecp, vecp2;                 // pythia four-vector object to store four-momenta
//   double p;                         // momentum |p| of the parton 
//   int counter = 0;                  // counter will provide unique color ID for produced partons
//   int numEvents = 0;                // will store number of actually accepted events (some are rejected by fragmentation)
//   int mt;                           // maximal time steps
//   int posy;                         // bin number 
//   int p_imax;                       // maximum number of particles in list
//   int id;                           // parton ID
//   int pisum = 0;                    
//   int piPlusSum = 0, piMinusSum=0;;                    
//   int KPlusSum = 0, KMinusSum=0;                    
//   const int bins = 40;              // number of bins
//   const int phiBins = 6;            // number of bins
//   double scale = 20.;               // maximum p_t in GeV in the binning (bin size = scale/bins)
//   double tscale = PI;               // maximum theta
//   double rscale = 10.;              // maximum r
//   double xscale = 20.;              // maximum x
//   int binning[bins];                // array with all the bins
//   int binning_piplus[bins];         // array with all the bins
//   int binning_Kplus[bins];          // array with all the bins
//   int binning_piminus[bins];         // array with all the bins
//   int binning_Kminus[bins];          // array with all the bins
//   int binning_tmp[bins];            // array with all the bins
//   int binning_sq[bins];             // array with all the bins for squares for error
//   int binning_hm[bins];             // array with all the bins
//   int binning_g[bins];              // array with all the bins for gluons
//   int binning_q[bins];              // array with all the bins for quarks
//   int binning_gamma[bins];          // array with all the bins for photons
//   int binning_gamma_sq[bins];       // array with all the bins for photons for squares for error
//   int binning_gamma_tmp[bins];      // array with all the bins for photons
//   int binning_r[bins];              // array with all the bins for the radius in the transverse plane
//   double binning_x[bins][bins];        // array with all the bins for the xy positions in the transverse plane
//   int binning_phi_pt[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_tmp[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_sq[phiBins][bins];// array with bins in phi and pt for squares for error
//   int binning_phi_pt_g[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_tmp_g[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_sq_g[phiBins][bins];// array with bins in phi and pt for squares for error
//   int binning_phi_pt_q[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_tmp_q[phiBins][bins];// array with bins in phi and pt
//   int binning_phi_sq_q[phiBins][bins];// array with bins in phi and pt for squares for error
//   double binning_mom[bins][bins];   
//   int binning_theta[bins];          // array with all the bins for the angle theta
//   int totalSum = 0;
//   int otherPartons = 0;                
//   int hmsum = 0;
//   int totalNNs;
//   double ymax = 1; //0.5 for one unit of rapidity around y=0.
//   double pt2=0.;
//   double pt2r=0.;
//   double Etot=0;
//   double NColTot=0.;
//   double qtTotAll=0.;
//   double Erun;
//   double Ejet=10.;
//   double r, theta;
//   int posr, postheta, posxi, posyi, posphi;
//   // make them both quarks with x GeV for testing:
  
//   if ( fixedEnergy == 1 )
//     {
//       jp1.id(1); jp2.id(1);
//       jp1.p(Ejet,0.,0.); jp2.p(Ejet,0.,0.);
//       jp1.col(101); jp1.acol(102); jp2.col(103); jp2.acol(104);
//       jp1.x(4.); jp1.y(0.); jp1.z(0.);
//       jp2.x(4.); jp2.y(0.); jp2.z(0.);
//       //jp1.splits(0); jp2.splits(0);
//     }
  
//   // init binning
//   for(int iy=0; iy<bins; iy++)
//     {
//       binning_r[iy]=0;
//       binning_theta[iy]=0;
//       binning[iy]=0;
//       binning_piplus[iy]=0;
//       binning_Kplus[iy]=0;
//       binning_piminus[iy]=0;
//       binning_Kminus[iy]=0;
//       binning_tmp[iy]=0;
//       binning_sq[iy]=0;
//       binning_hm[iy]=0;
//       binning_q[iy]=0;
//       binning_g[iy]=0;
//       binning_gamma[iy]=0;
//       binning_gamma_sq[iy]=0;
//       binning_gamma_tmp[iy]=0;
//     }

//   for(int iy=0; iy<bins; iy++)
//     for(int ix=0; ix<bins; ix++)
//       {
// 	binning_x[iy][ix]=0.;
// 	binning_mom[iy][ix]=0.;
//       }

//   for(int iy=0; iy<bins; iy++)
//     for(int ix=0; ix<phiBins; ix++)
//       {
// 	binning_phi_sq[ix][iy]=0.;
// 	binning_phi_pt[ix][iy]=0.;
// 	binning_phi_tmp[ix][iy]=0.;
// 	binning_phi_sq_q[ix][iy]=0.;
// 	binning_phi_pt_q[ix][iy]=0.;
// 	binning_phi_tmp_q[ix][iy]=0.;
// 	binning_phi_sq_g[ix][iy]=0.;
// 	binning_phi_pt_g[ix][iy]=0.;
// 	binning_phi_tmp_g[ix][iy]=0.;
//       }


//   mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  
//   ofstream foutt("test.dat",ios::app); 

//   for (int j=0; j<runs; j++)      // loop over all events
//     {
//       qtTot=0.;
//       NCol=0;
//       pt2=0.;
//       Erun=0.;
//       counter = 0;                // reset counter in the beginning of every event
     

//       plist[0]->clear();          // clear the parton list
//       plist[1]->clear();          // clear auxiliary list (obsolete at the moment)
     
//       if( fixedEnergy == 0 )
// 	{
// 	  totalNNs=generateEvent(plist[0]); // version that samples number of collisions with Glauber model
// 	}
//       else
// 	{
// 	  plist[0]->push_back( jp1 );
// 	  plist[0]->push_back( jp2 );
// 	}
            
//       if (evolution == 1)         // evolve in medium if settings allow it

// 	{
// 	  for(int i=0; i<mt; i++) // loop over all time steps 
// 	    {
// 	      counter = evolve(plist, counter, i);
// 	      counter+=1;
// 	    }
// 	}

//       // bin partons:
//       p_imax = plist[0]->size();
//       for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
// 	{
// 	  double pl,pt,theta, En, y, xini, yini;
// 	  id = plist[0]->at(p_i).id();
	  
// 	  pt = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
// 	  pl = plist[0]->at(p_i).p().pz(); //p_long
// 	  En = sqrt(pt*pt+pl*pl);
// 	  theta = atan(pt/pl);
// 	  //if (theta<0) cout << "theta=" << theta << endl;
// 	  if (theta<0) theta = PI+theta;
// 	  postheta = floor(theta*(bins/tscale));

// 	  r = sqrt(pow(plist[0]->at(p_i).xini(),2.)+pow(plist[0]->at(p_i).yini(),2.));
// 	  xini = plist[0]->at(p_i).xini();
// 	  yini = plist[0]->at(p_i).yini();
// 	  //cout << "x=" << xini << endl;
// 	  //cout << "y=" << yini << endl;
// 	  posr = floor(r*(bins/rscale));
// 	  posxi = floor((xini+10)*(bins/xscale));
// 	  posyi = floor((yini+10)*(bins/xscale));
// 	  //cout << "posxi=" << posxi << endl;
// 	  //cout << "posyi=" << posyi << endl;

// 	  //countSplits+=plist[0]->at(p_i).splits();
// 	  //avSplits+=plist[0]->at(p_i).splits()/plist[0]->size();
// 	  //if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[0]->at(p_i).splits();
// 	  // if (p>500.) countHardGluons++;
	  
// 	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
	  
// 	  if ( pt>0.) // bin initial positions of all partons that are within the intersting y range 
// 	                              //and above a certain p_t //// abs(y)<=ymax && 
// 	    {
// 	      //if( (id > 0 && id < 4) || id == 21 )
// 	      //{
// 		  if(posxi>0 && posxi<bins && posyi>0 && posyi<bins) 
// 		    {
// 		      binning_x[posxi][posyi]+=1./static_cast<double>(p_imax);
// 		      binning_mom[posxi][posyi]+=plist[0]->at(p_i).pini().px();
// 		      //cout << "adding one at " << posxi << ", " << posyi << endl;
// 		    }
// 		  if(posr>0 && posr<bins) binning_r[posr]+=1;
// 		  if(postheta>0 && postheta<bins) binning_theta[postheta]+=1;
// 		  //}
// 	    }
// 	}

//       p_imax = plist[0]->size();
//       for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
// 	{
// 	  double pl,pt,phi;
// 	  double En;
// 	  double y; // rapidity
// 	  id = plist[0]->at(p_i).id();
// 	  p = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
// 	  pt = sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
// 	  Erun += sqrt(pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.)+pow(plist[0]->at(p_i).p().pz(),2.));
// 	  pt2 += (pow(plist[0]->at(p_i).p().px(),2.)+pow(plist[0]->at(p_i).p().py(),2.));
// 	  posy = floor(pt*(bins/scale));
// 	  //cout << "pt=" << pt << ", posy=" << posy << endl;
// 	  pl = plist[0]->at(p_i).p().pz(); //p_long
// 	  En = sqrt(pt*pt+pl*pl);
// 	  y = 0.5*log((En+pl)/(En-pl)); //rapidity
// 	  phi = asin(sqrt(pow(plist[0]->at(p_i).p().py(),2.))/pt);   // azimuthal angle
// 	  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi

// 	  //cout << "p=" << p << endl;
// 	  //countSplits+=plist[0]->at(p_i).splits();
// 	  //avSplits+=plist[0]->at(p_i).splits()/plist[0]->size();
// 	  //if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[0]->at(p_i).splits();
// 	  // if (p>500.) countHardGluons++;
	  
// 	  if (fixedEnergy) // don't use a rapidity cut for the fixed energy calculation
// 	    {
// 	      if( id > 0 && id < 4 ) if(posy>=0 && posy<bins) binning_q[posy]+=1;
// 	      if( id == 21 ) if(posy>=0 && posy<bins) binning_g[posy]+=1;
// 	    }
// 	  else if ( abs(y)<=ymax )
// 	    {
// 	      if( id > 0 && id < 4 )
// 		{
// 		  if(posy>=0 && posy<bins) 
// 		    {
// 		      binning_q[posy]+=1;
// 		      if(posphi>=0 && posphi<phiBins)
// 			binning_phi_pt_q[posphi][posy] += 1;
// 		    }
// 		}
// 	      else if( id == 21 )
// 		{
// 		  if(posy>=0 && posy<bins) 
// 		    {
// 		      binning_g[posy]+=1;
// 		      if(posphi>=0 && posphi<phiBins)
// 			binning_phi_pt_g[posphi][posy] += 1;
// 		    }
// 		}
// 	    }
// 	}
      
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  for(int ix=0; ix<phiBins; ix++)
// 	    {
// 	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy]-binning_phi_tmp_q[ix][iy];
// 	      binning_phi_sq_q[ix][iy]+=binning_phi_tmp_q[ix][iy]*binning_phi_tmp_q[ix][iy];
// 	      binning_phi_tmp_q[ix][iy]=binning_phi_pt_q[ix][iy];
// 	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy]-binning_phi_tmp_g[ix][iy];
// 	      binning_phi_sq_g[ix][iy]+=binning_phi_tmp_g[ix][iy]*binning_phi_tmp_g[ix][iy];
// 	      binning_phi_tmp_g[ix][iy]=binning_phi_pt_g[ix][iy];
// 	    }
// 	}

//       pt2r+=pt2/p_imax;
//       Etot+=Erun/p_imax;
//       NColTot+=static_cast<double>(NCol)/static_cast<double>(p_imax);
//       qtTotAll+=qtTot/p_imax;

//       // fragmentation
//       if (fragmentationSwitch == 1)
// 	{

// 	  if ( fullEvent == 0 ) 
// 	    {
// 	      numEvents+=fragmentation( plist );
// 	      double pl;
// 	      double En;
// 	      double y; // rapidity
// 	      double eta; //pseudo-rapidity
// 	      double phi; // angle with respect to the reaction plane
// 	      p_imax = pythia.event.size();
// 	      for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
// 		{
// 		  if (pythia.event[p_i].isFinal()) totalSum++;
// 		  id = pythia.event[p_i].id();
// 		  // count pi_0s (111) or pi+ (211)
// 		  if ( id == 111 ) // pions (pi0)
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      phi = asin(sqrt(pow(pythia.event[p_i].py(),2.))/p);        // azimuthal angle
// 		      posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
// 		      //cout << "p_T=" << p << ", px=" << pythia.event[p_i].px() << ", py=" << pythia.event[p_i].py() 
// 		      //     << ", p_T=" << sqrt(pow(pythia.event[p_i].px(),2.)+pow(pythia.event[p_i].py(),2.)) 
// 		      //     << ", phi=" << phi << ", phi_deg=" << phi/PI*180. 
// 		      //     << ", posphi=" << posphi << endl;
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning[posy] += 1;
// 			  pisum += 1;
// 			  if(posphi>=0 && posphi<phiBins)
// 			    binning_phi_pt[posphi][posy] += 1;
// 			}
// 		    }
// 		  if ( id == 211 ) // pions (pi+)
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      phi = asin(sqrt(pow(pythia.event[p_i].py(),2.))/p);        // azimuthal angle
// 		      posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_piplus[posy] += 1;
// 			  piPlusSum += 1;
// 			}
// 		    }
// 		  if ( id == -211 ) // pions (pi-)
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_piminus[posy] += 1;
// 			  piMinusSum += 1;
// 			}
// 		    }
// 		  if ( id == 321 ) // Kaons (K+)
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_Kplus[posy] += 1;
// 			  KPlusSum += 1;
// 			}
// 		    }
// 		  if ( id == -321 ) // Kaons (K-)
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_Kminus[posy] += 1;
// 			  KMinusSum += 1;
// 			}
// 		    }
// 		  if ( id == 22 && ( pythia.event[p_i].status()<90 || pythia.event[p_i].status()>99 ) ) // photons
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_gamma[posy] += 1;
// 			}
// 		    }
// 		  if ( id < 0 && pythia.event[p_i].isHadron() ) // hadrons
// 		    {		
// 		      p = pythia.event[p_i].pT();                                // p_trans
// 		      posy = floor(p*(bins/scale));                              // bin number
// 		      pl = pythia.event[p_i].pz();                               // p_long
// 		      En = pythia.event[p_i].e();                                // energy
// 		      y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 		      //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 		      if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			{
// 			  binning_hm[posy] += 1;
// 			  hmsum += 1;
// 			}
// 		    }
// 		}
// 	    }
// 	  else 
// 	    {
// 	      for (int i=0; i<totalNNs; i++)
// 		{
// 		  numEvents+=fragmentation( plist, i );
// 		  double pl;
// 		  double En;
// 		  double y; // rapidity
// 		  double eta; //pseudo-rapidity
// 		  double phi; // angle with respect to the reaction plane
// 		  p_imax = pythia.event.size();
// 		  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
// 		    {
// 		      if (pythia.event[p_i].isFinal()) totalSum++;
// 		      id = pythia.event[p_i].id();
// 		      // count pi_0s (111) or pi+ (211)
// 		      if ( id == 111 ) // pions (pi0)
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  phi = asin(sqrt(pow(pythia.event[p_i].py(),2.))/p);        // azimuthal angle
// 			  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
// 			  //cout << "p_T=" << p << ", px=" << pythia.event[p_i].px() << ", py=" << pythia.event[p_i].py() 
// 			  //     << ", p_T=" << sqrt(pow(pythia.event[p_i].px(),2.)+pow(pythia.event[p_i].py(),2.)) 
// 			  //     << ", phi=" << phi << ", phi_deg=" << phi/PI*180. 
// 			  //     << ", posphi=" << posphi << endl;
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning[posy] += 1;
// 			      pisum += 1;
// 			      if(posphi>=0 && posphi<phiBins)
// 				binning_phi_pt[posphi][posy] += 1;
// 			    }
// 			}
// 		      if ( id == 211 ) // pions (pi+)
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  phi = asin(sqrt(pow(pythia.event[p_i].py(),2.))/p);        // azimuthal angle
// 			  posphi = floor(phi*(phiBins/(PI/2.)));                     // bin number in phi
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_piplus[posy] += 1;
// 			      piPlusSum += 1;
// 			    }
// 			}
// 		      if ( id == -211 ) // pions (pi-)
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_piminus[posy] += 1;
// 			      piMinusSum += 1;
// 			    }
// 			}
// 		      if ( id == 321 ) // Kaons (K+)
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_Kplus[posy] += 1;
// 			      KPlusSum += 1;
// 			    }
// 			}
// 		      if ( id == -321 ) // Kaons (K-)
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_Kminus[posy] += 1;
// 			      KMinusSum += 1;
// 			    }
// 			}
// 		      if ( id == 22 && ( pythia.event[p_i].status()<90 || pythia.event[p_i].status()>99 ) ) // photons
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_gamma[posy] += 1;
// 			    }
// 			}
// 		      if ( id < 0 && pythia.event[p_i].isHadron() ) // hadrons
// 			{		
// 			  p = pythia.event[p_i].pT();                                // p_trans
// 			  posy = floor(p*(bins/scale));                              // bin number
// 			  pl = pythia.event[p_i].pz();                               // p_long
// 			  En = pythia.event[p_i].e();                                // energy
// 			  y = 0.5*log((En+pl)/(En-pl));                              // rapidity
// 			  //eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
// 			  if(posy>=0 && posy<bins && abs(y)<=ymax)
// 			    {
// 			      binning_hm[posy] += 1;
// 			      hmsum += 1;
// 			    }
// 			}
// 		    }
// 		}
// 	    }
// 	}
      
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  binning_tmp[iy]=binning[iy]-binning_tmp[iy];
//        	  binning_sq[iy]+=binning_tmp[iy]*binning_tmp[iy];
// 	  binning_tmp[iy]=binning[iy];
// 	  binning_gamma_tmp[iy]=binning_gamma[iy]-binning_gamma_tmp[iy];
//        	  binning_gamma_sq[iy]+=binning_gamma_tmp[iy]*binning_gamma_tmp[iy];
// 	  binning_gamma_tmp[iy]=binning_gamma[iy];
// 	  for(int ix=0; ix<phiBins; ix++)
// 	    {
// 	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy]-binning_phi_tmp[ix][iy];
// 	      binning_phi_sq[ix][iy]+=binning_phi_tmp[ix][iy]*binning_phi_tmp[ix][iy];
// 	      binning_phi_tmp[ix][iy]=binning_phi_pt[ix][iy];
// 	    }
// 	}

//       if(j%1000==0)
// 	{
// 	  cout << "#" << j  << endl; // write every 1000th time step
// 	  foutt << "#" << j << endl;
// 	}
//     } // end loop over all events
//   foutt.close();
      
//   // output parton positions:
//   fstream foutth("theta.dat",ios::out); 
//   foutth.precision(12);  
//   foutth << "> Number_of_events = " << runs << endl; 
//   foutth << "> quarks " << endl;
//   for(int iy=0; iy<bins; iy++)
//     {
//       foutth << (iy)/(bins/tscale) << " " << bins*static_cast<double>(binning_theta[iy])/tscale << endl; //devide by runs later
//     }
//   foutth.close();

//   fstream foutx("xy.dat",ios::out); 
//   foutx.precision(12);  
//   foutx << "> Number_of_events = " << runs << endl; 
//   foutx << "> y,x " << endl;
//   for(int iy=0; iy<bins; iy++)
//     for(int ix=0; ix<bins; ix++)
//       {
// 	foutx << (iy)/(bins/xscale) - 10 << " " << (ix)/(bins/xscale) -10 << " " << bins*bins*binning_x[iy][ix]/xscale/xscale
// 	      << " " << bins*bins*static_cast<double>(binning_mom[iy][ix])/xscale/xscale
// 	      << endl; //devide by runs later
// 	if (ix == bins-1) foutx << endl;
//       }
//   foutx.close();

//   fstream foutr("r.dat",ios::out); 
//   foutr.precision(12);  
//   foutr << "> Number_of_events = " << runs << endl; 
//   foutr << "> r " << endl;
//   for(int iy=0; iy<bins; iy++)
//     {
//       foutr << (iy)/(bins/rscale) << " " << bins*static_cast<double>(binning_r[iy])/rscale  << endl; //divide by runs later
//     }
//   foutr.close();

//   fstream foutl("x.dat",ios::out); 
//   foutl.precision(12);  
//   foutl << "> Number_of_events = " << runs << endl; 
//   foutl << "> x " << endl;
//   for(int ix=0; ix<bins; ix++)
//     {
//       foutl << (ix)/(bins/xscale) - 10 << " " << bins*static_cast<double>(binning_x[ix][20])/xscale  << endl; //divide by runs later
//     }
//   foutl.close();

//   fstream fouty("y.dat",ios::out); 
//   fouty.precision(12);  
//   fouty << "> Number_of_events = " << runs << endl; 
//   fouty << "> y " << endl;
//   for(int iy=0; iy<bins; iy++)
//     {
//       fouty << (iy)/(bins/xscale) - 10 << " " << bins*static_cast<double>(binning_x[20][iy])/xscale  << endl; //divide by runs later
//     }
//   fouty.close();

  
//   pt2r/=(runs*maxTime);
//   Etot/=(runs);
//   qtTotAll/=runs;
//   NColTot/=runs;



//   // output partons:
//   fstream foutp("partons.dat",ios::out); 
//   foutp.precision(12);  
//   foutp << "> Number_of_events = " << runs << endl; 
//   foutp << "> quarks gluons" << endl;
//   for(int iy=0; iy<bins; iy++)
//     {
//       foutp <<  ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " << bins*static_cast<double>(binning_q[iy])/scale/2. //devide by runs later
// 	    << " " << bins*static_cast<double>(binning_g[iy])/scale/2. << endl; //devide by runs later
//       //divide by 2 since we have 2 initial particles.
//     }
//   foutp.close();
//   //cout << "total number of splittings per run=" << countSplits/runs 
//   //	   << ", splittings per parton=" << avSplits/runs << endl;
//   //cout << "total number of splittings of the initial partons per run=" << countSplitsOfInitialPartons/runs << endl;
//   //cout << "hard gluons with p>500 GeV = " << countHardGluons << endl;
  
//   fstream foutphiq("quarkphi.dat",ios::out); 
//   foutphiq << "> quarks. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
//   foutphiq.precision(12);  
//   foutphiq << "> Number_of_events = " << numEvents << endl; 
//   for(int iy=0; iy<bins; iy++)
//     {
//       foutphiq << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[0][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[1][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[2][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[3][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[4][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_q[5][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_q[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
// 		<< endl; //devide by numEvents later 
//       // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
//       // second output is the sum of squares to calculate the error in the end. 
//     }
//   foutphiq.close();

//   fstream foutphig("gluonphi.dat",ios::out); 
//   foutphig << "> gluons. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
//   foutphig.precision(12);  
//   foutphig << "> Number_of_events = " << numEvents << endl; 
//   for(int iy=0; iy<bins; iy++)
//     {
//       foutphig << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[0][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[1][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[2][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[3][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[4][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 		<< bins*static_cast<double>(binning_phi_pt_g[5][iy])/scale/(2.*ymax) << " " 
// 		<< bins*bins*static_cast<double>(binning_phi_sq_g[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
// 		<< endl; //devide by numEvents later 
//       // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
//       // second output is the sum of squares to calculate the error in the end. 
//     }
//   foutphig.close();

//   cout << endl;
//   if (fragmentationSwitch==1)
//     {
//       fstream fout("pi0.dat",ios::out); 
//       fout << "> pi0:" << endl;
//       fout.precision(12);  
//       fout << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  fout << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning[iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
// 	       << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       fout.close();

//       fstream foutphi("pi0phi.dat",ios::out); 
//       foutphi << "> pi0. phi (value, error) = 0-15 15-30 30-45 45-60 60-75 75-90" << endl;
//       foutphi.precision(12);  
//       foutphi << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutphi << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_phi_pt[0][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[0][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 	       << bins*static_cast<double>(binning_phi_pt[1][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[1][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 	       << bins*static_cast<double>(binning_phi_pt[2][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[2][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 	       << bins*static_cast<double>(binning_phi_pt[3][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[3][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 	       << bins*static_cast<double>(binning_phi_pt[4][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[4][iy])/scale/(2.*ymax)/scale/(2.*ymax)<< " " 
// 	       << bins*static_cast<double>(binning_phi_pt[5][iy])/scale/(2.*ymax) << " " 
// 	       << bins*bins*static_cast<double>(binning_phi_sq[5][iy])/scale/(2.*ymax)/scale/(2.*ymax)
// 	       << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutphi.close();
      
//       fstream fout2("hminus.dat",ios::out); 
//       fout2 << endl << "> h-:" << endl;
//       fout2.precision(12);  
//       fout2 << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  fout2 << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_hm[iy])/scale/(2.*ymax) << endl; 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	}
//       fout2.close();

//       fstream foutgamma("photons.dat",ios::out); 
//       foutgamma << "> gamma:" << endl;
//       foutgamma.precision(12);  
//       foutgamma << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutgamma << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 		    << bins*static_cast<double>(binning_gamma[iy])/scale/(2.*ymax) << " " 
// 		    << bins*bins*static_cast<double>(binning_gamma_sq[iy])/scale/(2.*ymax)/scale/(2.*ymax)
// 		    << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutgamma.close();

//       fstream foutpiplus("pi+.dat",ios::out); 
//       foutpiplus << "> pi+:" << endl;
//       foutpiplus.precision(12);  
//       foutpiplus << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutpiplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_piplus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutpiplus.close();
//       fstream foutpiminus("pi-.dat",ios::out); 
//       foutpiminus << "> pi-:" << endl;
//       foutpiminus.precision(12);  
//       foutpiminus << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutpiminus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_piminus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutpiminus.close();

//       fstream foutKplus("K+.dat",ios::out); 
//       foutKplus << "> K+:" << endl;
//       foutKplus.precision(12);  
//       foutKplus << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutKplus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_Kplus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutKplus.close();
//       fstream foutKminus("K-.dat",ios::out); 
//       foutKminus << "> K-:" << endl;
//       foutKminus.precision(12);  
//       foutKminus << "> Number_of_events = " << numEvents << endl; 
//       for(int iy=0; iy<bins; iy++)
// 	{
// 	  foutKminus << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " 
// 	       << bins*static_cast<double>(binning_Kminus[iy])/scale/(2.*ymax) << endl; //devide by numEvents later 
// 	  // divide by 2*ymax to get 1 unit of rapidity (\Delta y = 2*ymax). (scale/bins)=\Delta p_t = bin width. 
// 	  // second output is the sum of squares to calculate the error in the end. 
// 	}
//       foutKminus.close();
//     }  
//   cout.precision(6);

//   cout << "K+/pi+=" << static_cast<double>(KPlusSum)/static_cast<double>(piPlusSum) <<endl;
//   cout << "K-/pi-=" << static_cast<double>(KMinusSum)/static_cast<double>(piMinusSum) <<endl;
//   cout << "K+/K-=" << static_cast<double>(KPlusSum)/static_cast<double>(KMinusSum) <<endl;
//   cout << "sigmaGen=" << pythia.info.sigmaGen() << endl;
//   cout << "pisum=" << pisum << endl;
//   cout << "sum of negatively charged hadrons=" << hmsum << endl;
//   cout << "used events=" << numEvents << endl;
//   cout << "time steps=" << mt << endl;
//   cout << "total number of particles=" << totalSum << endl; 
//   cout << "total number of partons != u,d,s,g (mainly diquarks with p_t < 1 GeV):" << otherPartons << endl; 
//   cout << "qhat =" << pt2r << " GeV^2/fm" << endl;
//   cout << "<q_T^2>/mfp =" << qtTotAll/maxTime << " GeV^2/fm" << endl;
//   cout << "dE/dx =" << (Ejet-Etot)/maxTime << " GeV/fm" << endl;
//   cout << "number of collisions =" << NColTot << endl;
//   cout << "partons: " << p_imax << endl;
//   cout << "mfp=" << maxTime/NColTot << endl;
//   delete plist[0];
//   delete plist[1];
//   delete plist;
// }


//Here, we add the routines needed for heavy quark and quarkonium production:

double MARTINI::CornellPotential(double r) {

  //In GeV:                                                                                                                                                                    
  if(r < 0.2) r = 0.2 ;
  double potential = -0.0985/r+0.812*r
    -(-0.0985/hydroTau0+0.812*hydroTau0);
  //A cutoff for large separations:                                                                         
  //if(potential > 0.88){
  //potential = 0.88 ; }
  if(potential > 1.45){
    potential = 1.45 ; }

  return potential ;
}


double MARTINI::F(double r, double T){
  double F;
  if(r < 0.2){
    F = 0.;
  }
  else{
    double c_term = 0.14657/(T/T_C_HQ-0.98) ;
    double mu_term = 0.0301+0.0063*T/T_C_HQ ;
    double exp_term = exp(-mu_term*r*r/(hbarc*hbarc)) ;
    
    //The magnitude of the force, from a fit to quenched
    //QCD data:
    F = 31.841*mu_term*exp_term
      +0.6199*exp_term/(r*r)
      +1.1258*exp_term-57.823*r*r*mu_term*exp_term
      +7.528*mu_term*r*exp_term/(T/T_C_HQ-0.98) ;
    F = hbarc*F;
  }
  return F;
}

int MARTINI::generateTestEvent_HQ(vector<Parton> *plist, double L){

  Parton parton1, parton2, parton3, parton4;

  //We add 4 partons to plist, all charm quarks and anti-quarks, at the corners of a square.
  //This will test the fragmentation. Here we instantiate the plist, in some cases with properties not
  //entirely relevant, but done for completeness:
  
  totalNNs = 2;

  parton1.id(4);
  parton2.id(-4);
  parton3.id(4);
  parton4.id(-4);

  parton1.status(1);
  parton2.status(1);
  parton3.status(1);
  parton4.status(1);

  parton1.mass(CHARM_MASS);
  parton2.mass(CHARM_MASS);
  parton3.mass(CHARM_MASS);
  parton4.mass(CHARM_MASS);

  parton1.x(0.5*L);
  parton2.x(-0.5*L);
  parton3.x(-0.5*L);
  parton4.x(0.5*L);

  parton1.y(0.5*L);
  parton2.y(-0.5*L);
  parton3.y(0.5*L);
  parton4.y(-0.5*L);

  parton1.z(0.);
  parton2.z(0.);
  parton3.z(0.);
  parton4.z(0.);

  parton1.xini(0.5*L);
  parton2.xini(-0.5*L);
  parton3.xini(-0.5*L);
  parton4.xini(0.5*L);

  parton1.yini(0.5*L);
  parton2.yini(-0.5*L);
  parton3.yini(0.5*L);
  parton4.yini(-0.5*L);

  parton1.zini(0.);
  parton2.zini(0.);
  parton3.zini(0.);
  parton4.zini(0.);

  parton1.tini(0.);
  parton2.tini(0.);
  parton3.tini(0.);
  parton4.tini(0.);

  parton1.col(1);
  parton2.col(1);
  parton3.col(0);
  parton4.col(0);

  parton1.acol(0);
  parton2.acol(0);
  parton3.acol(1);
  parton4.acol(1);

  parton1.frozen(0);
  parton2.frozen(0);
  parton3.frozen(0);
  parton4.frozen(0);

  parton1.p(0., 0., 0.);
  parton2.p(0., 0., 0.);
  parton3.p(0., 0., 0.);
  parton4.p(0., 0., 0.);

  Vec4 pvec(CHARM_MASS, 0., 0., 0.);
  parton1.pini(pvec);
  parton2.pini(pvec);
  parton3.pini(pvec);
  parton4.pini(pvec);

  //Tricking the fragmentation routine into thinking the quarks evolved in the medium:
  parton1.tFinal(1.);
  parton2.tFinal(1.);
  parton3.tFinal(1.);
  parton4.tFinal(1.);

  parton1.eventNumber(0);
  parton2.eventNumber(0);
  parton3.eventNumber(0);
  parton4.eventNumber(0);

  parton1.elasticCollisions(0);
  parton2.elasticCollisions(0);
  parton3.elasticCollisions(0);
  parton4.elasticCollisions(0);

  parton1.source(0);
  parton2.source(0);
  parton3.source(0);
  parton4.source(0);

  parton1.antiI(1);
  parton2.antiI(0);
  parton3.antiI(3);
  parton4.antiI(2);

  plist->push_back(parton1);
  plist->push_back(parton2);
  plist->push_back(parton3);
  plist->push_back(parton4);

  return 1;
}

int MARTINI::generateEvent_HQ(vector<Parton> *plist)
{
  if (fullEvent==1)
    {
      // Sample Nu (which gives the number of nucleons at radial position r) for both nuclei
      // then lay them on top of each other and sample the number of charm and bottom events with 
      // area = totalHQXSec.
      
      //cout << "Number of binary collisions = " << glauber->TAB() << endl; 
      //cout << " A=" << glauber->nucleusA() << endl;

      double A;
      double Z;
      Parton parton;                               // new Parton object
      double b = glauberImpactParam;
      int posx;
      int posy;
      int n1=0;
      int n2=0;
      double r = 1.2*pow(glauber->nucleusA(),1./3.);

      nucleusA.clear();
      nucleusB.clear();
      sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons

      double cellLength = sqrt(inelasticXSec*0.1); // compute cell length in fm (1 mb = 0.1 fm^2)
      double latticeLength = 4.*r;                 // spread the lattice 2r in both (+/-) directions (total length = 4r)
      int ixmax = ceil(latticeLength/cellLength);
      ixmax*=2;
      double xmin = -latticeLength;
      
      int nucALat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus A
      int nucBLat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus B
      int collLat[ixmax][ixmax];                   // lattice with the positions of the interactions -> for information only

      int countCollisionsA[100];
      int countCollisionsB[100];

      for (int i=0; i<100; i++)
	{
	  countCollisionsA[i]=0;
	  countCollisionsB[i]=0;
	}

      for (int i = 0; i < ixmax; i++)              // initialize cells: zero nucleons in every cell
	for (int j = 0; j < ixmax; j++)
	  {
	    nucALat[i][j] = 0;
	    nucBLat[i][j] = 0;
	    collLat[i][j] = 0;
	  }

      for (int i = 0; i < nucleusA.size(); i++)    // fill nucleons into cells
	{
	  posx = floor((nucleusA.at(i).x-xmin-b/2.)/cellLength);
	  posy = floor((nucleusA.at(i).y-xmin)/cellLength);
	  nucALat[posx][posy] += 1;
	      
	  posx = floor((nucleusB.at(i).x-xmin+b/2.)/cellLength);
	  posy = floor((nucleusB.at(i).y-xmin)/cellLength);
	  nucBLat[posx][posy] += 1;

	}

      int done;
      int numberOfpp;
      int eventNumber=0;
      numberOfpp=0;
      //cout << "jetXSec=" << jetXSec << endl;
      //cout << "inelasticXSec=" << inelasticXSec << endl;
      A=glauber->nucleusA();
      if (A==197.) //gold - Au 
	Z=79.; 
      else if (A==63.) //copper - Cu
	Z=29.;
      else if (A==208.) //lead - Pb
	Z=82.;
      else
	{
	  Z=0;
	}
      for (int ix = 0; ix < ixmax; ix++) 
	for (int iy = 0; iy < ixmax; iy++)
	  {
	    for (int i = 0; i < nucALat[ix][iy]; i++)
	      for (int j = 0; j < nucBLat[ix][iy]; j++)
		{
		  //cout << "ix+b/2=" << static_cast<int>(ix+b/2.) << " # there=" << nucALat[static_cast<int>(ix+b/2.)][iy] << endl;
		  if ( random->genrand64_real1() < totalHQXSec/inelasticXSec )
		    {
		      countCollisionsA[i]+=1;
		      countCollisionsB[j]+=1;
		      double xPositionInCell = (random->genrand64_real1())*cellLength;
		      double yPositionInCell = (random->genrand64_real1())*cellLength;
		      // see if colliding nucleons are neutrons
		      if ( random->genrand64_real1() > Z/A ) n1=1;
		      else n1=0;
		      if ( random->genrand64_real1() > Z/A ) n2=1;
		      else n2=0;

		      int pythiaWorked = 0;
		      while(pythiaWorked == 0){
			pythiaWorked = pythia.next(n1,n2);
		      }                                                // generate event with pythia

		      done = 0;
		      for (int ip = 0; ip < pythia.event.size(); ++ip) 
			{
			  if (pythia.event[ip].status()>0 && (abs(pythia.event[ip].id())==4 || abs(pythia.event[ip].id()) == 5) )
			    // If the parton is final, and is a heavy quark, then put it into the list.
			    {
			      parton.id(pythia.event[ip].id());                     // set parton id
			      parton.status(pythia.event[ip].status());             // set parton status
			      parton.mass(pythia.event[ip].m());                // set mass
			      if (fixedTemperature==0)
				{
				  //Random numbers added to the initial positions:
				  double dx, dy, dz;
				  dx = dy = dz = 0.;
				  if(abs(parton.id()) == 4){
				    dx = gsl_ran_gaussian(gsl_rand, charmWidth);
				    dy = gsl_ran_gaussian(gsl_rand, charmWidth);
				  }
				  if(abs(parton.id()) == 5){
				    dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
				    dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
				  }
				  //cout << "moveBeforeTau0 = " << moveBeforeTau0 << endl;
				  if( moveBeforeTau0 == 1){
				    Vec4 partonP = pythia.event[ip].p();
				    double pzP = partonP.pz();
				    double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
						     +partonP.px()*partonP.px()
						     +partonP.py()*partonP.py()
						     +partonP.pz()*partonP.pz() );
				    double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
				    dx += partonP.px()*dt/EP;
				    dy += partonP.py()*dt/EP;
				    dz += partonP.pz()*dt/EP;
				    //cout << "pz = " << pzP << ", EP = " << EP << ", vz = " << pzP/EP << endl;
				    //cout << "dz = " << dt << ", dt = " << dt << ", vz = " << dz/dt << endl;
				    //cout << "dtau = " << sqrt(dt*dt-dz*dz) << endl;
				  }

				  parton.x(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
				  parton.y(xmin+(iy)*cellLength+yPositionInCell+dy);
				  parton.xini(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
				  parton.yini(xmin+(iy)*cellLength+yPositionInCell+dy);
				  parton.z(dz);
				  parton.zini(dz);
				  parton.tini(hydroTau0);
				  if ( done == 0 )
				    {
				      posx = floor((parton.x()-xmin)/cellLength);
				      posy = floor((parton.y()-xmin)/cellLength);
				      collLat[posx][posy] += 1;
				                  
				      //cout << numberOfpp << endl;
				      numberOfpp++;
				      done = 1;
				    }
				}
			      else
				{
				  double dx, dy, dz;
				  dx = dy = dz = 0.;
				  if(abs(parton.id()) == 4){
				    dx = gsl_ran_gaussian(gsl_rand, charmWidth);
				    dy = gsl_ran_gaussian(gsl_rand, charmWidth);
				  }
				  if(abs(parton.id()) == 5){
				    dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
				    dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
				  }

				  if( moveBeforeTau0 == 1){
				    Vec4 partonP = pythia.event[ip].p();
				    double pzP = partonP.pz();
				    double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
						     +partonP.px()*partonP.px()
						     +partonP.py()*partonP.py()
						     +partonP.pz()*partonP.pz() );
				    double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
				    dx += partonP.px()*dt/EP;
				    dy += partonP.py()*dt/EP;
				    dz += partonP.pz()*dt/EP;
				  }

				  parton.x(dx);     
				  parton.y(dy);
				  parton.xini(dx);     
				  parton.yini(dy);
				  parton.z(dz);
				  parton.zini(dz);
				  parton.tini(0.);
				}
			      parton.col(pythia.event[ip].col());                 // set color 
			      parton.acol(pythia.event[ip].acol());                 // set anti-color
			      parton.frozen(0);                                     // parton is not frozen (will evolve)
			      if(setInitialPToZero == 0){
				parton.p(pythia.event[ip].p());  // set momentum
				parton.pini(pythia.event[ip].p());                    // set initial momentum
			      }
			      else{
				parton.p(0.,0.,0.);
				Vec4 pinitial(pythia.event[ip].mass(),0.,0.,0.);
				parton.pini(pinitial);                    // set initial momentum
			      }
			      parton.tFinal(0.);                                    // Set the initial freezeout time to zero
			      parton.eventNumber(eventNumber);
			      parton.elasticCollisions(0);                          // set initial no. of el colls to zero
			      parton.source(0);                                     // All initial partons have itsSource=0
			      plist->push_back(parton);                             // add the parton to the main list
			    }
			}
		      eventNumber++;
		    }
		}
	  }

      //Now, determine the partners of each heavy quark:
      for(int i = 0; i < plist->size(); i++){
	if(i%2 == 0){
	  plist->at(i).antiI(i+1);
	}
	else
	  plist->at(i).antiI(i-1);
      }

      totalNNs = eventNumber;
      //output for plot
      ofstream fout1("./output/density.dat",ios::out); 
      ofstream fout2("./output/densityB.dat",ios::out); 
      ofstream fout3("./output/colls.dat",ios::out); 
      ofstream fout4("./output/NN.dat",ios::app); 
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
	    if ( j == ixmax-1 ) fout1 << endl;
	  }
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
	    if ( j == ixmax-1 ) fout2 << endl;
	  }
      
      for (int i = 0; i < ixmax; i++)
	for (int j = 0; j < ixmax; j++)
	  {
	    fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
	    if ( j == ixmax-1 ) fout3 << endl;
	  }
      
      fout4 << numberOfpp << endl;
      
      fout1.close();
      fout2.close();
      fout3.close();
      fout4.close();
      return eventNumber;
    }
  else if (fullEvent==0) // sample only one hard collision per heavy-ion event
    {
      Parton parton;                                                // new Parton object
      double A;
      double Z;
      int n1;
      int n2;
      A=glauber->nucleusA();
      if (A==197.) //gold - Au 
	Z=79.; 
      else if (A==63.) //copper - Cu
	Z=29.;
      else if (A==208.) //lead - Pb
	Z=82.;
      else
	{
	  Z=0;
	}
      if ( random->genrand64_real1() > Z/A ) n1=1;
      else n1=0;
      if ( random->genrand64_real1() > Z/A ) n2=1;
      else n2=0;

      int pythiaWorked = 0;
      while(pythiaWorked == 0){
	pythiaWorked = pythia.next(n1,n2);                                                // generate event with pythia
      }

      ReturnValue rv;
      //if (evolution && fixedTemperature==0) 
      if (!allFromCenter)
	{
	  if(nbinFromFile == 1){
	    //Sample a random element of the list of number of collisions. -CFY 11/2/2010
	    double randomr = random->genrand64_real2();
	    int randomint = (int)(randomr*((double)Nbin));
	    rv.x = binary[randomint][0];
	    rv.y = binary[randomint][1];
	  }
	  else if (glauberEnvelope)
	    rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
	  else
	    rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
	  cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
	}
      else
	{
	  rv.x=0.;
	  rv.y=0.;
	}
      if (trackHistory)
	{
	  rv.x=initialXjet;
	  rv.y=initialYjet;
	}
      for (int i = 0; i < pythia.event.size(); ++i) 
	{
	  if (pythia.event[i].isFinal() && (abs(pythia.event[i].id()) == 4 || abs(pythia.event[i].id()) == 5) )
	    {
	      double dx, dy, dz;
	      dx = dy = dz = 0.;
	      if(abs(pythia.event[i].id() ) == 4){
		dx = gsl_ran_gaussian(gsl_rand, charmWidth);
		dy = gsl_ran_gaussian(gsl_rand, charmWidth);
	      }
	      if(abs(pythia.event[i].id() ) == 5){
		dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
		dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
	      }

	      if( moveBeforeTau0 == 1){
		Vec4 partonP = pythia.event[i].p();
		double pzP = partonP.pz();
		double EP = sqrt(pythia.event[i].m()*pythia.event[i].m()
				 +partonP.px()*partonP.px()
				 +partonP.py()*partonP.py()
				 +partonP.pz()*partonP.pz() );
		double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
		dx += partonP.px()*dt/EP;
		dy += partonP.py()*dt/EP;
		dz += partonP.pz()*dt/EP;
	      }

	      parton.id(pythia.event[i].id());                      // set parton id
	      parton.status(pythia.event[i].status());              // set parton status
	      parton.mass(pythia.event[i].m());                // set mass
	                  
	      parton.x(rv.x+dx);                                   // set position
	      parton.y(rv.y+dy);
	      if (fixedTemperature==0)
		parton.tini(hydroTau0);
	      else
		parton.tini(0.);
	      parton.xini(rv.x+dx);                                // set initial position to remember
	      parton.yini(rv.y+dy);
	                  
	      parton.z(dz);
	      parton.zini(dz);
	      parton.col(pythia.event[i].col());                 // set color 
	      parton.acol(pythia.event[i].acol());                  // set anti-color
	      parton.frozen(0);                                     // parton is not frozen (will evolve)
	      if(setInitialPToZero == 0){
		parton.p(pythia.event[i].p());                        // set momentum
		parton.pini(pythia.event[i].p());                     // set momentum
	      }
	      else{
		parton.p(0.,0.,0.);                                   // set momentum
		Vec4 pinitial(pythia.event[i].mass(), 0.,0.,0.);
		parton.pini(pinitial);                               // set momentum
	      }
	      parton.tFinal(0.);                                    // Set the freezeout time initially to zero
	      parton.elasticCollisions(0);                          // set initial no. of el colls to zero
	      //parton.splits(0);                                   // set initial no. of splittings to zero
	      parton.source(0);                                     // All intial partons have itsSource=0
	      plist->push_back(parton);                             // add the parton to the main list
	    }
	}

      //Now, determine the partners of each heavy quark:
      for(int i = 0; i < plist->size(); i++){
	if(i%2 == 0){
	  plist->at(i).antiI(i+1);
	}
	else
	  plist->at(i).antiI(i-1);
      }

      return 1;
    }
}

// Evolves heavy quarks according to Langevin dynamics. The evolution here is 
// pretty reasonable, if you want to see something bad, look at generateEvent_HQ. -CFY 4/10/2011
int MARTINI::evolve_HQ(vector<Parton>  **plist, int counter, int it)
{

  cout.precision(6);

  HydroInfo hydroInfo, antihydroInfo;

  Vec4 vecp, antivecp;                          // momentum four-vectors
  double pmu[4], pmuRest[4], pmuRestNew[4], pmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 
  double antipmu[4], antipmuRest[4], antipmuRestNew[4], antipmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 

  int imax = plist[0]->size();        // number of partons in the list
  int id, antiI; //newID;                      // original and new parton's ID
  //int col, acol, newCol, newAcol;     // number of the color string attached to the old, and new parton
  int ix, iy, iz, itau;               // parton's position in cell coordinates
  int antiix, antiiy, antiiz, antiitau; // parton's position in cell coordinates
  int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction

  double dt = dtfm/hbarc;             // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
                                      // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
  double dtflow, antidtflow;                      // dt in the fluid's rest frame
  double T, antiT, T_ave;             // temperature
  double M, antiM;                    // Parton's mass [GeV]
  double x, y, z;                     // parton's position in [fm]
  double antix, antiy, antiz       ;  // parton's position in [fm]
  double t, tau, antitau;             // lab time and lab tau
  double vx, vy, vz;                  // flow velocity of the cell
  double antivx, antivy, antivz;      // flow velocity of the cell
  double beta;                        // absolute value of the flow velocity
  double antibeta;                        // absolute value of the flow velocity
  double gamma;                       // gamma factor 
  double antigamma;                       // gamma factor 
  double umu[4], umub[4];             // Fluid's 4-velocity
  double antiumu[4], antiumub[4];             // Fluid's 4-velocity
  const double rd=12.3;               // ratio of degrees of freedom g_QGP/g_H (for Kolb hydro only)
  
  double kappa, sigma, eta;
  double kicktriplet[3];
  double r, force, F_x, F_y, F_z;

  ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);

  if ( hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);

  //cout << "counter = " << counter << ", it = " << it << endl;

  for ( int i=0; i<imax; i++)                                   // loop over all partons
    { 
      //cout << "plist[0]->at(" << i << ").frozen() = " << plist[0]->at(i).frozen() << endl;
      if(plist[0]->at(i).frozen() == 0){
	counter++;                                                // counter that characterizes every step uniquely
	id = plist[0]->at(i).id();                                // id of parton i 
	
	if(abs(id)!=4 && abs(id)!=5) continue;                    // We immediately leave the loop of the parton is not a heavy quark!
	
        
	t = it*dtfm;        //had a +tau0 here. now let partons travel doing their vacuum shower stuff first...
	
	vecp = plist[0]->at(i).p();                               // four-vector p of parton i
	M = plist[0]->at(i).mass();                               // Mass of the parton (only used for charm and bottom quarks)
	pmu[1] = vecp.px();
	pmu[2] = vecp.py();
	pmu[3] = vecp.pz();
	pmu[0] = sqrt(M*M+pmu[1]*pmu[1]+pmu[2]*pmu[2]+pmu[3]*pmu[3]); //vecp loaded into an array for ease with LorentzBooster
	
	x = plist[0]->at(i).x();                                  // x value of position in [fm]
	y = plist[0]->at(i).y();                                  // y value of position in [fm]
	z = plist[0]->at(i).z();                                  // z value of position in [fm]
	
	// boost - using original position...
	ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
	iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
	// note that x and y run from -hydroXmax to +hydroXmax
	// and ix and iy from 0 to 2*hydroXmax
	// hence the (hydroXmax+x or y) for both
	if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) 
	  {
	    iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
	  }
	
	if (fixedTemperature == 0)
	  {
	    tau = 0.;           //just give tau some value - will change below
	    
	    if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
	    //cout << "tau = " << tau << endl;
	    if (tau<hydroTau0||tau>hydroTauMax-1) continue;       // do not evolve if tau<tau0 or tau>tauMax
	    
	    // get the temperature and flow velocity at the current position and time:
	    hydroInfo = hydroSetup->getHydroValues(x, y, z, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
						     hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
	    
	    T = hydroInfo.T;
	    //cout << "T = " << T << "for particle " << i << endl;
	    vx = hydroInfo.vx;
	    vy = hydroInfo.vy;
	    vz = hydroInfo.vz;
	    
	    //cout << "At t = " << t << " and z = " << z << ", vz = " << vz << endl;
	    //cout << "For the charm quark, vz = " << pmu[3]/pmu[0] << endl;
	    if (!trackHistory)
	      {
		// warn if out of range and move on to next parton (should not happen!)
		if ( ix < 0 || ix >= ixmax ) 
		  {
		    cout << "WARNING - x out of range" << endl;
		    plist[0]->at(i).frozen(1);
		    continue;
		  }
		if ( iy < 0 || iy >= ixmax )
		  {
		    cout << "WARNING - y out of range" << endl;
		    plist[0]->at(i).frozen(1);
		    continue;
		  }
		if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) && ( iz < 0 || iz >= izmax ) ) 
		  {
		    cout << "WARNING - z out of range" << endl;
		    plist[0]->at(i).frozen(1);
		    continue;
		  }
	      }
	    
	    if (T<hydroTfinal && hydroWhichHydro == 2)                 // for 2D hydro evolutions stop at T_final (can change that)
	      {
		plist[0]->at(i).tFinal(t);
		continue;                                             // do not evolve if T<T_c
	      }
	    
	    beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
	    gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
	    umu[0] = gamma;
	    umu[1] = gamma*vx;
	    umu[2] = gamma*vy;
	    umu[3] = gamma*vz;
	    LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
	    dtflow = dt*pmuRest[0]/pmu[0];                          // dtfm in the fluid's rest frame
	    
	  }
	else // the fixed temperature case:
	  {
	    beta = 0.;
	    gamma = 1.; 
	    T = fixedT;
	    umu[0] = 1.; umu[1] = umu[2] = umu[3] = 0.;
	    LorentzBooster(umu, pmu, pmuRest);
	    dtflow = dt*pmuRest[0]/pmu[0];
	  }  
	
	//After all of the various hydro switches, we need to evolve the heavy quarks. This is optimized either for 
	//quarkonium or for single heavy quarks. We differentiate between the two with the switch examineHQ. First, for single 
	//heavy quarks:
	if(examineHQ == 1){
	  //Finally evolve the heavy quarks' positions:
	  plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	  plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	  plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	  
	  // Sample a triplet of kicks for the heavy quark. For now, we 
	  // take the simplest case: constant 2piT*D_c, however later 
	  // both kicktriplet[3] and eta will be determined from HTL 
	  // results at LO or NLO order:
	  if(T > T_C_HQ){
	    plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
	    kappa = 4.*M_PI*T*T*T/TwoPiTD_HQ;
	    sigma = sqrt(kappa*dtflow);
	    kicktriplet[0] = gsl_ran_gaussian(gsl_rand, sigma);
	    kicktriplet[1] = gsl_ran_gaussian(gsl_rand, sigma);
	    kicktriplet[2] = gsl_ran_gaussian(gsl_rand, sigma);
	    
	    eta = 2.*M_PI*T*T/(3.*M);
	    
	    //Update the momentum in the rest frame:
	    pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
	    pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
	    pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
	    pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
				 +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
	    
	    //The flow velocity for boosting back:
	    umub[0] = umu[0];
	    umub[1] = -umu[1];
	    umub[2] = -umu[2];
	    umub[3] = -umu[3];
	    
	    //The boost back into the lab frame:
	    LorentzBooster(umub, pmuRestNew, pmuNew);
	    
	    //Update the particle's momentum:
	    plist[0]->at(i).p(pmuNew[1], pmuNew[2], pmuNew[3]);
	  }
	  else{ plist[0]->at(i).frozen(1);}
      }
	
	//Next, we examine quarkonium:
	if(examineHQ == 2){
	  if(id == 4 || id == 5){ //Only the quarks, not the anti-quarks!
	    antiI = plist[0]->at(i).antiI();
	    //We must duplicate the steps at the beginning of the program for the partner anti-quarks position, momentum, and 
	    //other properties:
	    antivecp = plist[0]->at(antiI).p();                               // four-vector p of parton i
	    antiM = plist[0]->at(antiI).mass();                               // Mass of the parton (only used for charm and bottom quarks)
	    antipmu[1] = antivecp.px();
	    antipmu[2] = antivecp.py();
	    antipmu[3] = antivecp.pz();
	    antipmu[0] = sqrt(antiM*antiM+antipmu[1]*antipmu[1]+antipmu[2]*antipmu[2]+antipmu[3]*antipmu[3]); //vecp loaded into an array for ease with LorentzBooster
	    
	    antix = plist[0]->at(antiI).x();                                  // x value of position in [fm]
	    antiy = plist[0]->at(antiI).y();                                  // y value of position in [fm]
	    antiz = plist[0]->at(antiI).z();                                  // z value of position in [fm]
	    
	    antiix = floor((hydroXmax+antix)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
	    antiiy = floor((hydroXmax+antiy)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
	    if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) 
	      {
		antiiz = floor((hydroZmax+antiz)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
	      }
	    
	    // get the temperature and flow velocity at the current position and time:
	    
	    if(fixedTemperature == 0){
	      antitau = 0.;                                                 // Initially set antitau = 0.
	      
	      //cout << "antiz = " << antiz << ", t = " << t << endl;
	      if (antiz*antiz<t*t) antitau = sqrt(t*t-antiz*antiz);         // determine tau (z*z) should always be less than (t*t)
	      //cout << "OK before the call continue statement..." << endl;
	      if (antitau<hydroTau0 || antitau>hydroTauMax-1) continue;       // do not evolve if tau<tau0 or tau>tauMax
	      //cout << "x = " << x << ", y = " << y << ", z = " << z << ", t = " << t << endl;
	      //cout << "antix = " << antix << ", antiy = " << antiy << ", antiz = " << antiz << ", t = " << t << endl;

	      antihydroInfo = hydroSetup->getHydroValues(antix, antiy, antiz, t, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
							 hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
	      
	      if (!trackHistory)
		{
		  // warn if out of range and move on to next parton (should not happen!)
		  if ( antiix < 0 || antiix >= ixmax ) 
		    {
		      cout << "WARNING - x out of range" << endl;
		      plist[0]->at(i).frozen(1);
		      plist[0]->at(antiI).frozen(1);
		      continue;
		    }
		  if ( antiiy < 0 || antiiy >= ixmax )
		    {
		      cout << "WARNING - y out of range" << endl;
		      plist[0]->at(i).frozen(1);
		      plist[0]->at(antiI).frozen(1);
		      continue;
		  }
		  if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) && ( antiiz < 0 || antiiz >= izmax ) ) 
		    {
		      cout << "WARNING - z out of range" << endl;
		      plist[0]->at(i).frozen(1);
		      plist[0]->at(antiI).frozen(1);
		      continue;
		    }
		}
	      
	      antiT = antihydroInfo.T;
	      T_ave = 0.5*(T+antiT);
	      
	      antivx = antihydroInfo.vx;
	      antivy = antihydroInfo.vy;
	      antivz = antihydroInfo.vz;
	      
	      antibeta = sqrt(antivx*antivx+antivy*antivy+antivz*antivz);                           // absolute value of flow velocity in units of c
	      antigamma = 1./sqrt(1.-antibeta*antibeta);                            // gamma factor
	      antiumu[0] = antigamma;
	      antiumu[1] = antigamma*antivx;
	      antiumu[2] = antigamma*antivy;
	      antiumu[3] = antigamma*antivz;
	      
	      LorentzBooster(antiumu, antipmu, antipmuRest);                      // Boost into the fluid's rest frame
	      antidtflow = dt*antipmuRest[0]/antipmu[0];                          // dtfm in the fluid's rest frame
	      
	    }
	    else // the fixed temperature case:
	      {
		antibeta = 0.;
		antigamma = 1.; 
		antiT = fixedT; T_ave = fixedT;
		
		antiumu[0] = 1.; antiumu[1] = antiumu[2] = antiumu[3] = 0.;
		LorentzBooster(antiumu, antipmu, antipmuRest);
		antidtflow = dt*antipmuRest[0]/antipmu[0];
		antix = plist[0]->at(antiI).x();                                  // x value of position in [fm]
		antiy = plist[0]->at(antiI).y();                                  // y value of position in [fm]
		antiz = plist[0]->at(antiI).z();                                  // z value of position in [fm]
	      }
	    // We sample two triplets of kicks now:
	    
	    //Finally, we take into account the heavy quark-antiquark interaction. This is done in the lab frame, 
	    //and is correct strictly non-relativistically. This is good for observables such as TOTAL yields, and 
	    //relative yields of excited states, but is probably very bad for something like a differential yield of 
	    //J/psi up to high p_T
	    if(T_ave > T_C_HQ){
	      plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
	      plist[0]->at(antiI).tFinal(t); // update the time - finally will be the final time
	      
	      plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	      plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	      plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
	      plist[0]->at(antiI).x(antix+antivecp.px()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
	      plist[0]->at(antiI).y(antiy+antivecp.py()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
	      plist[0]->at(antiI).z(antiz+antivecp.pz()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
	      
	      kappa = 4.*M_PI*T*T*T/TwoPiTD_HQ;
	      sigma = sqrt(kappa*dtflow);
	      eta = 0.5*kappa/(M*T);
	      
	      kicktriplet[0] = gsl_ran_gaussian(gsl_rand, sigma);
	      kicktriplet[1] = gsl_ran_gaussian(gsl_rand, sigma);
	      kicktriplet[2] = gsl_ran_gaussian(gsl_rand, sigma);
	      
	      //Update the momentum in the rest frame:
	      pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
	      pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
	      pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
	      pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
				   +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
	      
	      //The flow velocity for boosting back:
	      umub[0] = umu[0];
	      umub[1] = -umu[1];
	      umub[2] = -umu[2];
	      umub[3] = -umu[3];
	      
	      //The boost back into the lab frame:
	      LorentzBooster(umub, pmuRestNew, pmuNew);
	      
	      //Now, the anti-quark:
	      kappa = 4.*M_PI*antiT*antiT*antiT/TwoPiTD_HQ;
	      sigma = sqrt(kappa*dtflow);
	      eta = 0.5*kappa/(M*antiT);
	      
	      kicktriplet[0] = gsl_ran_gaussian(gsl_rand, sigma);
	      kicktriplet[1] = gsl_ran_gaussian(gsl_rand, sigma);
	      kicktriplet[2] = gsl_ran_gaussian(gsl_rand, sigma);
	      
	      //Update the momentum in the rest frame:
	      antipmuRestNew[1] = antipmuRest[1]*(1.-eta*antidtflow)+kicktriplet[0];
	      antipmuRestNew[2] = antipmuRest[2]*(1.-eta*antidtflow)+kicktriplet[1];
	      antipmuRestNew[3] = antipmuRest[3]*(1.-eta*antidtflow)+kicktriplet[2];
	      antipmuRestNew[0] = sqrt(M*M+antipmuRestNew[1]*antipmuRestNew[1]
				       +antipmuRestNew[2]*antipmuRestNew[2]+antipmuRestNew[3]*antipmuRestNew[3]);
	      
	      //The flow velocity for boosting back:
	      antiumub[0] = antiumu[0];
	      antiumub[1] = -antiumu[1];
	      antiumub[2] = -antiumu[2];
	      antiumub[3] = -antiumu[3];
	      
	      //The boost back into the lab frame:
	      LorentzBooster(antiumub, antipmuRestNew, antipmuNew);
	      
	      //Finally, we include the effect of the HQ potential, parametrized from
	      // lattice data:
	      r = sqrt( (x-antix)*(x-antix)+(y-antiy)*(y-antiy)+(z-antiz)*(z-antiz) );
	      force = F(r, T_ave);
	      //force = 0.; //We are debugging!
	      //if(r > 0.4 && r < 1.5){
	      //force = 0.0985/(r*r)+0.812;
	      //}//Still debugging!

	      if(force > 0.){
		F_x = force*(antix-x)/r;
		F_y = force*(antiy-y)/r;
		F_z = force*(antiz-z)/r;
	      }
	      else F_x = F_y = F_z = 0.;

	      //Update the particle's momentum:
	      plist[0]->at(i).p(pmuNew[1]+F_x*dt*hbarc, pmuNew[2]+F_y*dt*hbarc, pmuNew[3]+F_z*dt*hbarc);
	      plist[0]->at(antiI).p(antipmuNew[1]-F_x*dt*hbarc, antipmuNew[2]-F_y*dt*hbarc, antipmuNew[3]-F_z*dt*hbarc);
	      
	      //Update the anti-particle's momentum:
	      //cout << "pmu1: " << pmuNew[1] << " " << pmuNew[2] << " " << pmuNew[3] << endl;
	      //cout << "pmu2: " << antipmuNew[1] << " " << antipmuNew[2] << " " << antipmuNew[3] << endl;
	    }
	    else{
	      plist[0]->at(i).frozen(1);
	      plist[0]->at(antiI).frozen(1);
	    }
	  }
	}
      }
    }
  //foutt.close();
  return counter;
}

//How many neighbors do the heavy quarks have?
//int MARTINI::examine_nearest_neighbor_HQ(vector<Parton> ** plist){
//if(fullEvent == 0 || totalNNs <= 1){
//  cout << "Not enough heavy quarks to make an interesting analysis here..." << endl;
//}
//if(fullEvent == 1 || totalNNs > 1){
//  int np = plist[0]->size();
//  int neighbors = 0;
//  for(int ip=0; ip<np; ip++){
//    int id = plist[0]->at(ip).id();
//    int antiI;
//    if(id>0){
//if(ip%2 == 0){
//  antiI = ip+1;
//}
//else antiI = ip-1;
//
//for(int ip2=ip+1; ip2<np; ip++){
	  

int MARTINI::fragmentation_HQ(vector<Parton> ** plist){
  pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below

  //If fullEvent == 0, this is relatively simple: we only determine the invariant mass of the pair 
  //in the CM frame and use this to hadronize the pairs into mesons:
  if(fullEvent == 0 || totalNNs == 1){
    int id1 = plist[0]->at(0).id();
    int id2 = plist[0]->at(1).id();
    //The separation of the quarks in the lab frame:
    double dx = plist[0]->at(0).x()-plist[0]->at(1).x();
    double dy = plist[0]->at(0).y()-plist[0]->at(1).y();
    double dz = plist[0]->at(0).z()-plist[0]->at(1).z();
    double r = sqrt(dx*dx+dy*dy+dz*dz);
    //After kinetic freezeout, the evolution of the heavy quarks is completely adiabatic. Therefore, 
    //the Cornell potential should be used to determine the binding energy of the pair:
    double potential;
    potential = CornellPotential(r);
    //else potential = CornellPotential(0.);
    //The relative momentum of the pair:
    Vec4 p1 = plist[0]->at(0).p();
    Vec4 p2 = plist[0]->at(1).p();
    double M1 = plist[0]->at(0).mass();
    double M2 = plist[0]->at(1).mass();
    //// While the 1S definition of the charm and bottom mass is appropriate for dynamics, 
    //// the \bar{MS} definition is more appropriate for hadronization in the color evaporation model:
    //double M1, M2;
    //if(abs(id1) == 4){
    //M1 = CHARM_MASS;
    //}
    //else
    //M1 = BOTTOM_MASS;
    //if(abs(id2) == 4){
    //M2 = CHARM_MASS;
    //}
    //else
    //M2 = BOTTOM_MASS;
    double p1mu[4], p2mu[4];
    p1mu[1] = p1.px();
    p1mu[2] = p1.py();
    p1mu[3] = p1.pz();
    p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
    p2mu[1] = p2.px();
    p2mu[2] = p2.py();
    p2mu[3] = p2.pz();
    p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);

    //Let's just calculate with the invariant mass and the potential:
    double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
			-(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
			-(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
			-(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );

    //double E_CM = prel2/mred+potential;
    double E_CM = M_inv+potential;
    //Right now, generateEvent_HQ will never generate b\bar{c} states when fullEvent==0, however
    //we code this case in here in case this change is made:
    if(abs(id1)==4 && abs(id2)==4){
      //If the binding energy falls within certain ranges, either open charm or charmonium states 
      //are formed. We distinguish between three cases: pairs which we say hadronize into J/Psi particles
      //(and therefore have particle ID 443), particles which are bound yet which we will not specify at this 
      //point their state (which have the special particle ID ), and open heavy flavor mesons, whose momenta are 
      //determined with Peterson fragmentation:
      if(E_CM < E_JPSI){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MJPsi*MJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(443, 81, 0, 0, pnewvec, MJPsi);
      }
      if(E_CM > E_JPSI && E_CM < E_CCB_BOUND){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MExcitedJPsi*MExcitedJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(443000, 81, 0, 0, pnewvec, MExcitedJPsi);
      }
      if(E_CM > E_CCB_BOUND){
	double z1, z2;
	z1 = SamplePetersonFunction(M1);
	z2 = SamplePetersonFunction(M2);
	double pnew1[4], pnew2[4];
	pnew1[1] = z1*p1mu[1];
	pnew1[2] = z1*p1mu[2];
	pnew1[3] = z1*p1mu[3];
	pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
	pnew2[1] = z2*p2mu[1];
	pnew2[2] = z2*p2mu[2];
	pnew2[3] = z2*p2mu[3];
	pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
	Vec4 pnew1vec, pnew2vec;
	pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	//For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	//fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
	pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
      }
    }
    //Now, for dealing with bottom quarks:
    if(abs(id1)==5 && abs(id2)==5){
      if(E_CM < E_Y){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MUpsilon*MUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(553, 81, 0, 0, pnewvec, MUpsilon);
      }
      if(E_CM > E_Y && E_CM < E_BBB_BOUND){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MExcitedUpsilon*MExcitedUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(553000, 81, 0, 0, pnewvec, MExcitedUpsilon);
      }
      if(E_CM > E_BBB_BOUND){
	double z1, z2;
	z1 = SamplePetersonFunction(M1);
	z2 = SamplePetersonFunction(M2);
	double pnew1[4], pnew2[4];
	pnew1[1] = z1*p1mu[1];
	pnew1[2] = z1*p1mu[2];
	pnew1[3] = z1*p1mu[3];
	pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
	pnew2[1] = z2*p2mu[1];
	pnew2[2] = z2*p2mu[2];
	pnew2[3] = z2*p2mu[3];
	pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
	Vec4 pnew1vec, pnew2vec;
	pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	//For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	//fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
	pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
      }
    }
    if((abs(id1)==4 && abs(id2)==5) || (abs(id1)==5 && abs(id2)==4) ){
      if(E_CM < E_BBC){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MB_c*MB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(543, 81, 0, 0, pnewvec, MB_c);
      }
      if(E_CM > E_BBC && E_CM < E_BBC_BOUND){
	double pnewJPsi[4];
	pnewJPsi[1] = p1mu[1]+p2mu[1];
	pnewJPsi[2] = p1mu[2]+p2mu[2];
	pnewJPsi[3] = p1mu[3]+p2mu[3];
	pnewJPsi[0] = sqrt(MExcitedB_c*MExcitedB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	Vec4 pnewvec;
	pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	pythia.event.append(543000, 81, 0, 0, pnewvec, MExcitedB_c);
      }
      if(E_CM > E_BBC_BOUND){
	double z1, z2;
	z1 = SamplePetersonFunction(M1);
	z2 = SamplePetersonFunction(M2);
	double pnew1[4], pnew2[4];

	pnew1[1] = z1*p1mu[1];
	pnew1[2] = z1*p1mu[2];
	pnew1[3] = z1*p1mu[3];

	pnew2[1] = z2*p2mu[1];
	pnew2[2] = z2*p2mu[2];
	pnew2[3] = z2*p2mu[3];

	//For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	//fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	if(abs(id1)==4){
	  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
	  Vec4 pnew1vec;
	  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	  pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
	}
	else{
	  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
	  Vec4 pnew1vec;
	  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	  pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
	}
	if(abs(id2)==4){
	  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
	  Vec4 pnew2vec;
	  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	  pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
	}
	else{
	  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
	  Vec4 pnew2vec;
	  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	  pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
	}
      }
    }
  }
  if(fullEvent==1 && totalNNs > 1 ){
    //Determining recombinant production in a full event will be more difficult, but crucial. 
    //First, the element of the permutation group describing the pairing of heavy quarks. 
    //The quark from collision i, i=0, ..., totalNNs, is paired with the anti-quark in collision
    //pairing[i]:
    int pairing[totalNNs];
    for(int ipair=0; ipair<totalNNs; ipair++){
      pairing[ipair] = ipair;
    }

    //Whether or not the initial heavy quark pair is not hadronizing statistically. This can happen for two reasons:
    //either the pair forms quarkonium and does not interact strongly with the medium, or the pair never 
    //experienced in-medium evolution.
    int notstatistical[totalNNs];
    //int notevolved[totalNNs];
    for(int ipair=0; ipair<totalNNs; ipair++){
      notstatistical[totalNNs] = 0.;
      //notevolved[totalNNs] = 0.;

      //First, those pairs which do not hadronize statistically due to being bound:
      int id1, id2;
      id1 = plist[0]->at(2*ipair).id();
      id2 = plist[0]->at(2*ipair+1).id();

      double dx, dy, dz;
      dx = plist[0]->at(2*ipair).x()-plist[0]->at(2*ipair+1).x();
      dy = plist[0]->at(2*ipair).y()-plist[0]->at(2*ipair+1).y();
      dz = plist[0]->at(2*ipair).z()-plist[0]->at(2*ipair+1).z();

      double r = sqrt(dx*dx+dy*dy+dz*dz);
      double potential = CornellPotential(r);

      //The relative momentum of the pair:
      Vec4 p1 = plist[0]->at(2*ipair).p();
      Vec4 p2 = plist[0]->at(2*ipair+1).p();

      double M1 = plist[0]->at(2*ipair).mass();
      double M2 = plist[0]->at(2*ipair+1).mass();
      //double M1, M2;
      //if(abs(id1) == 4){
      //M1 = CHARM_MASS;
      //}
      //else
      //M1 = BOTTOM_MASS;
      //if(abs(id2) == 4){
      //M2 = CHARM_MASS;
      //}
      //else
      //M2 = BOTTOM_MASS;

      double p1mu[4], p2mu[4];
      p1mu[1] = p1.px();
      p1mu[2] = p1.py();
      p1mu[3] = p1.pz();
      p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
      p2mu[1] = p2.px();
      p2mu[2] = p2.py();
      p2mu[3] = p2.pz();
      p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);

      double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
			  -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
			  -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
			  -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );

      //Finally, the energy of the pair in the CM frame:
      double E_CM = M_inv+potential;
      //All that for this:
      if((abs(id1) ==4 && abs(id2) == 4) && E_CM<E_CCB_BOUND) notstatistical[ipair] = 1;
      if((abs(id1) ==5 && abs(id2) == 5) && E_CM<E_BBB_BOUND) notstatistical[ipair] = 1;
      if(( (abs(id1) ==4 && abs(id2) == 5) || (abs(id1) == 5 && abs(id2) == 4) ) && E_CM<E_BBC_BOUND) notstatistical[ipair] = 1;

      //Finally, whether or not the pair evolved in medium:
      //double z = 0.5*(plist[0]->at(2*ipair).z()+plist[0]->at(2*ipair+1).z());
      //double t_final = 0.5*(plist[0]->at(2*ipair).tFinal()+plist[0]->at(2*ipair+1).tFinal() );
      //double tau = sqrt(t_final*t_final-z*z);

      //Only need to check the final time of one particle to see if both were evolved, when we are working with
      //examineHQ==2:
      double t_final = plist[0]->at(2*ipair).tFinal();
      if(t_final == 0.){
	notstatistical[ipair] = 1;
	//notevolved[ipair] = 1;
      }
    }

    //We perform the Metropolis algorithm to determine the pairing of these heavy quarks:
    for(int iter=0; iter<N_SAMPLES; iter++){
      //Select a random 2-cycle:
      int cycle[2] ;
      cycle[0] = gsl_rng_uniform_int(gsl_rand, totalNNs);
      bool select_random = true;
      while(select_random){
	int random_int = gsl_rng_uniform_int(gsl_rand, totalNNs);
	if(random_int != cycle[0]){
	  cycle[1] = random_int;
	  select_random = false;
	}
      }
      
      //Skip the rest of the loop if the pair is not hadronizing statistically:
      if(notstatistical[cycle[0]]==1) continue;
      if(notstatistical[cycle[1]]==1) continue;
      
      int pairing_new[totalNNs];
      //The new pairing is mostly the same...
      for(int ipair=0; ipair<totalNNs; ipair++){
	pairing_new[ipair] = pairing[ipair];
      }
      //Except that it is multiplied on the left by the 2-cycle:
      pairing_new[cycle[0]]=pairing[cycle[1]];
      pairing_new[cycle[1]]=pairing[cycle[0]];
      
      //This is fine, but we need the exact locations of these particles in plist:
      int i1, i2, antiI1, antiI2;
      if(plist[0]->at(2*cycle[0]).id()>0) i1 = 2*cycle[0];
      else i1 = 2*cycle[0]+1;
      if(plist[0]->at(2*cycle[1]).id()>0) i2 = 2*cycle[1];
      else i2 = 2*cycle[1]+1;
      if(plist[0]->at(2*pairing[cycle[0]]).id()<0) antiI1 = 2*pairing[cycle[0]];
      else antiI1 = 2*pairing[cycle[0]]+1;
      if(plist[0]->at(2*pairing[cycle[1]]).id()<0) antiI2 = 2*pairing[cycle[1]];
      else antiI2 = 2*pairing[cycle[1]]+1;
      
      //Determine whether or not to accept or reject this pairing,
      //based on the Metropolis algorithm.
      //First, determine the difference in the total potential energy for the new and old pairings:
      double r1, r2, r1new, r2new;
      double dV;
      
      //We are going to reuse these dx variables for calculating the various r's and ultimately, dV:
      double dxtemp, dytemp, dztemp;
      dxtemp = plist[0]->at(i1).x()-plist[0]->at(antiI1).x();
      dytemp = plist[0]->at(i1).y()-plist[0]->at(antiI1).y();
      dztemp = plist[0]->at(i1).z()-plist[0]->at(antiI1).z();
      r1 = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
      dxtemp = plist[0]->at(i2).x()-plist[0]->at(antiI2).x();
      dytemp = plist[0]->at(i2).y()-plist[0]->at(antiI2).y();
      dztemp = plist[0]->at(i2).z()-plist[0]->at(antiI2).z();
      r2 = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
      dxtemp = plist[0]->at(i1).x()-plist[0]->at(antiI2).x();
      dytemp = plist[0]->at(i1).y()-plist[0]->at(antiI2).y();
      dztemp = plist[0]->at(i1).z()-plist[0]->at(antiI2).z();
      r1new = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
      dxtemp = plist[0]->at(i2).x()-plist[0]->at(antiI1).x();
      dytemp = plist[0]->at(i2).y()-plist[0]->at(antiI1).y();
      dztemp = plist[0]->at(i2).z()-plist[0]->at(antiI1).z();
      r2new = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);

      //cout << "r1 = " << r1 << ", r2 = " << r2 << ", r1new = " << r1new << ", r2new = " << r2new << endl;

      dV = CornellPotential(r1new)+CornellPotential(r2new)-CornellPotential(r1)-CornellPotential(r2);
      double BoltzmannFactor = exp(-dV/T_C_HQ);

      //The Metropolis algorithm: if the Boltzmann factor is greater than unity, immediately accept the permutation
      //(this is arbitrary but determines the next step):
      //cout << "dV = " << dV << endl;
      //cout << "T_C_HQ = " << T_C_HQ << endl;
      //cout << "Boltzmann factor = " << BoltzmannFactor << endl;
      if(BoltzmannFactor > 1){
	for(int ipair=0; ipair<totalNNs; ipair++){
	  pairing[ipair] = pairing_new[ipair];
	}
      }
      //If it is less than unity, accept the permutation at a rate that maintains detailed balance:
      else{
	double rdouble = gsl_rng_uniform(gsl_rand);
	//cout << "rdouble = " << rdouble << endl;
	if(rdouble < BoltzmannFactor){
	  for(int ipair=0; ipair<totalNNs; ipair++){
	    pairing[ipair] = pairing_new[ipair];
	  }
	}
      }
      ////Output the new permutation for testing:
      //for(int ipair=0; ipair<totalNNs; ipair++){
      //cout << pairing[ipair] << " ";
      //}
      //cout << endl;
    }

    ////Output the new permutation for testing:
    //cout << "The new permutation: " ;
    //for(int ipair=0; ipair<totalNNs; ipair++){
    //cout << pairing[ipair] << " ";
    //}
    //cout << endl;

    //Now that the pairings have been sampled, we loop through the full event, adding J/psi's, 
    //Upsilons, b\bar{c} states, c\bar{b} states, and unhadronized free heavy quarks to the 
    //event record:
    for(int ie=0; ie<totalNNs; ie++){
      int i1, i2;
      if(plist[0]->at(2*ie).id()>0) i1 = 2*ie;
      else i1 = 2*ie+1;
      if(plist[0]->at(2*pairing[ie]).id()<0) i2 = 2*pairing[ie];
      else i2 = 2*pairing[ie]+1;
      int id1 = plist[0]->at(i1).id();
      int id2 = plist[0]->at(i2).id();
      //Distinguish between diagonal and recombinant production:
      int status;
      if(ie == pairing[ie]) status = 81;
      else status = 82;
      
      //Determine again the binding energy in the CM frame:
      //The separation of the quarks in the lab frame:
      double dx = plist[0]->at(i1).x()-plist[0]->at(i2).x();
      double dy = plist[0]->at(i1).y()-plist[0]->at(i2).y();
      double dz = plist[0]->at(i1).z()-plist[0]->at(i2).z();
      double r = sqrt(dx*dx+dy*dy+dz*dz);
      //After kinetic freezeout, the evolution of the heavy quarks is completely adiabatic. Therefore, 
      //the Cornell potential should be used to determine the binding energy of the pair:
      double potential;
      //if(notevolved[ie] == 1){
      //potential = CornellPotential(0.);
      //}
      potential = CornellPotential(r);

      //The relative momentum of the pair:
      Vec4 p1 = plist[0]->at(i1).p();
      Vec4 p2 = plist[0]->at(i2).p();
      double M1 = plist[0]->at(i1).mass();
      double M2 = plist[0]->at(i2).mass();
      //double M1, M2;
      //if(abs(id1) == 4){
      //M1 = CHARM_MASS;
      //}
      //else
      //M1 = BOTTOM_MASS;
      //if(abs(id2) == 4){
      //M2 = CHARM_MASS;
      //}
      //else
      //M2 = BOTTOM_MASS;
      double p1mu[4], p2mu[4];
      p1mu[1] = p1.px();
      p1mu[2] = p1.py();
      p1mu[3] = p1.pz();
      p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
      p2mu[1] = p2.px();
      p2mu[2] = p2.py();
      p2mu[3] = p2.pz();
      p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);

      double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
			  -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
			  -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
			  -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );

      //Finally, the energy of the pair in the CM frame:
      double E_CM = M_inv+potential;

      if(abs(id1)==4 && abs(id2)==4){
	//If the binding energy falls within certain ranges, either open charm or charmonium states 
	//are formed. We distinguish between three cases: pairs which we say hadronize into J/Psi particles
	//(and therefore have particle ID 443), particles which are bound yet which we will not specify at this 
	//point their state (which have the special particle ID ), and open heavy flavor mesons, whose momenta are 
	//determined with Peterson fragmentation:
	//if(r < 2.0){
	//cout << "(dx, dy, dz) = (" << dx << ", " << dy << ", " << dz << ")" << endl;
	//cout << "tFinal of the charm quark is " << plist[0]->at(i1).tFinal() << endl;
	//cout << "tFinal of the charm anti-quark is " << plist[0]->at(i2).tFinal() << endl;
	//cout << "r for the pairing of charm quark " << i1 << " with charm anti-quark " << i2 << " is " << r << endl;
	//cout << "potential = " << potential << endl;
	//cout << "p for the charm quark is (" << p1.px() << ", " << p1.py() << ", " << p1.pz() << ")" << endl;
	//cout << "p for the charm anti-quark is (" << p2.px() << ", " << p2.py() << ", " << p2.pz() << ")" << endl;
	//cout << "p^mu for the charm quark is (" << p1mu[0] << ", " << p1mu[1] << ", " << p1mu[2] << ", " << p1mu[3] << ")" << endl;
	//cout << "p^mu for the charm anti-quark is (" << p2mu[0] << ", " << p2mu[1] << ", " << p2mu[2] << ", " << p2mu[3] << ")" << endl;
	//cout << "M_inv = " << M_inv << endl;
	//cout << "E_CM = " << E_CM << endl;
	//}

	if(E_CM < E_JPSI){
	  if(status == 82){
	    cout << "Recombinant J/Psi formed: tFinal(" << i1 << ") = " << plist[0]->at(i1).tFinal() << ", tFinal(" << i2 << ") = " << plist[0]->at(i2).tFinal() << endl;
	  }
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MJPsi*MJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(443, status, 0, 0, pnewvec, MJPsi);
	}
	if(E_CM > E_JPSI && E_CM < E_CCB_BOUND){
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MExcitedJPsi*MExcitedJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(443000, status, 0, 0, pnewvec, MExcitedJPsi);
	}
	if(E_CM > E_CCB_BOUND){
	  double z1, z2;
	  z1 = SamplePetersonFunction(M1);
	  z2 = SamplePetersonFunction(M2);
	  //cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
	  double pnew1[4], pnew2[4];
	  pnew1[1] = z1*p1mu[1];
	  pnew1[2] = z1*p1mu[2];
	  pnew1[3] = z1*p1mu[3];
	  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
	  pnew2[1] = z2*p2mu[1];
	  pnew2[2] = z2*p2mu[2];
	  pnew2[3] = z2*p2mu[3];
	  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
	  Vec4 pnew1vec, pnew2vec;
	  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);

	  //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	  //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	  pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
	  pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
	}
      }
      //Now, for dealing with bottom quarks:
      if(abs(id1)==5 && abs(id2)==5){
	if(E_CM < E_Y){
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MUpsilon*MUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(553, status, 0, 0, pnewvec, MUpsilon);
	}
	if(E_CM > E_Y && E_CM < E_BBB_BOUND){
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MExcitedUpsilon*MExcitedUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(553000, status, 0, 0, pnewvec, MExcitedUpsilon);
	}
	if(E_CM > E_BBB_BOUND){
	  double z1, z2;
	  z1 = SamplePetersonFunction(M1);
	  z2 = SamplePetersonFunction(M2);
	  //cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
	  double pnew1[4], pnew2[4];
	  pnew1[1] = z1*p1mu[1];
	  pnew1[2] = z1*p1mu[2];
	  pnew1[3] = z1*p1mu[3];
	  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
	  pnew2[1] = z2*p2mu[1];
	  pnew2[2] = z2*p2mu[2];
	  pnew2[3] = z2*p2mu[3];
	  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
	  Vec4 pnew1vec, pnew2vec;
	  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	  //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	  //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	  pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
	  pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
	}
      }
      if((abs(id1)==4 && abs(id2)==5) || (abs(id1)==5 && abs(id2)==4) ){
	if(E_CM < E_BBC){
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MB_c*MB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(543, status, 0, 0, pnewvec, MB_c);
	}
	if(E_CM > E_BBC && E_CM < E_BBC_BOUND){
	  double pnewJPsi[4];
	  pnewJPsi[1] = p1mu[1]+p2mu[1];
	  pnewJPsi[2] = p1mu[2]+p2mu[2];
	  pnewJPsi[3] = p1mu[3]+p2mu[3];
	  pnewJPsi[0] = sqrt(MExcitedB_c*MExcitedB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
	  Vec4 pnewvec;
	  pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
	  pythia.event.append(543000, status, 0, 0, pnewvec, MExcitedB_c);
	}
	if(E_CM > E_BBC_BOUND){
	  double z1, z2;
	  z1 = SamplePetersonFunction(M1);
	  z2 = SamplePetersonFunction(M2);	
	  //cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
	  double pnew1[4], pnew2[4];
	  
	  pnew1[1] = z1*p1mu[1];
	  pnew1[2] = z1*p1mu[2];
	  pnew1[3] = z1*p1mu[3];
	  
	  pnew2[1] = z2*p2mu[1];
	  pnew2[2] = z2*p2mu[2];
	  pnew2[3] = z2*p2mu[3];
	  
	  //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
	  //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
	  if(abs(id1)==4){
	    pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
	    Vec4 pnew1vec;
	    pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	    pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
	  }
	  else{
	    pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
	    Vec4 pnew1vec;
	    pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
	    pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
	  }
	  if(abs(id2)==4){
	    pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
	    Vec4 pnew2vec;
	    pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	    pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
	  }
	  else{
	    pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
	    Vec4 pnew2vec;
	    pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
	    pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
	  }
	}
      }
    }
  }
} 

// Uses the rejection method:
double MARTINI::SamplePetersonFunction(double M){
  double eps_Q;
  eps_Q = (CHARM_MASS/M)*(CHARM_MASS/M)*EPS_C;
  
  double z;
  bool zNotDetermined = true;
  while(zNotDetermined){

    double zSampled, fSampled;
    zSampled = gsl_rng_uniform(gsl_rand);
    fSampled = gsl_rng_uniform(gsl_rand)/eps_Q;

    double fValue = 1./(zSampled*(1.-1./zSampled-eps_Q/(1.-zSampled))*(1.-1./zSampled-eps_Q/(1.-zSampled)) );
    if(fSampled < fValue){
      z = zSampled;
      zNotDetermined = false;
    }
  }

  return z;
}

//// For the decays of open heavy flavor only:
//int MARTINI::decayOpenHeavyFlavor(){
//int eventSize = pythia.event.size();
//for(int ie=0; ie<eventSize; ie++){
