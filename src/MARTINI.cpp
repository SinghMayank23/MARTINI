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
#include <iomanip>

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
    binary_info_ptr = new binarycollision_info();
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
    delete binary_info_ptr;
}

// evolve every parton by one time step. This is the core of MARTINI.
int MARTINI::evolve(vector<Parton>  *plist, vector<Source> *slist, int counter, int it)
{
    HydroInfo hydroInfo;
    ReturnValue f;                  // will contain \Delta p, the change in momentum due to a process
    Norms n;                        // holds the integrals \int dGamma_i dp, currently three for all radiative process
    double qqRate, qgRate;          // hold the integrals \int dGamma_i domega for the elastic processes
    double gqRate, ggRate;
    double conversionqg;            // hold the total conversion rate q->g
    double conversionqgamma;        // hold the total conversion rate q->gamma
    double conversiongq;            // hold the total conversion rate g->q
    
    Vec4 vecp, newvecp;             // momentum four-vectors
    Vec4 vecpRest, newvecpRest;     // momentum four-vectors in rest frame of fluid cell
    Vec4 vecq;                      // four-vectors of momentum transfer
    Vec4 vecRecoil, vecRecoilRest;  // four-vectors of recoil momentum in fluid rest frame and lab frame
    Vec4 vecThermal, vecThermalLab;
    Parton newOne, newHole;         // parton object for the additionally created parton
    Source newSource;
    
    int imax = plist->size();       // number of partons in the list
    
    int id, newID;                  // original and new parton's ID
    int col, acol, newCol, newAcol; // number of the color string attached to the old, and new parton
    int ix, iy, iz, itau;           // parton's position in cell coordinates
    int ixmax, izmax;               // maximum cell number in x- and y-direction, and z-direction
    int radiate;                    // switches to hold what process will happen
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
    int GluonToQuark;              // holds the decision whether a gluon converts to a quark (=1) or anti-quark (=0)
    //Sangyong's addition
    int mother;
    int daughter;
    int temp_pos;
    int in_coherence;              // flag for whether mother and daughter are still together
    int ind_proc;                  // which process to take place
    
    double dt = dtfm/hbarc;        // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
                                   // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
    double T;                      // temperature
    double p, q, qt;               // parton momentum in [GeV]
    double x, y, z;                // parton's position in [fm]
    double px, py, pz;             // parton's momentum
    double t, tau, eta;            // lab time and lab tau
    double vx, vy, vz, veta;       // flow velocity of the cell
    double vetaZ, vetaZTau;        // flow velocity of the cell
    double beta;                   // absolute value of the flow velocity
    double gamma;                  // gamma factor 
    double cosPhi, cosPhiRest;     // angle between pvec and v, and angle between pvec_restframe and v
    double cosPhiRestEl;           // angle between pvec_restframe and v in case of elastic collisions with transv. mom. transfer
    double boostBack;              // boost factor to boost back from rest- to lab-frame
    double pRest;                  // rest frame abs. value of momentum
    double pxRest, pyRest, pzRest; // rest frame momentum components
    double pRecoil;                // lab frame abs. value of recoil particle
    double newpx, newpy, newpz;    // temporary momentum components
    double pxRecoil, pyRecoil, pzRecoil;
    double pxThermalLab, pyThermalLab, pzThermalLab;
    double omega;                  // energy transfered in elastic collision
    double totalQuarkProb = 0.;    // total probability that the quark undergoes some interaction
    double totalGluonProb = 0.;    // total probability that the gluon undergoes some interaction
    double randx;                  // catch the random number
    double EMProb;                 // temporary holder for photon radiation rate
    double AccProb;                // accumulated probability
    double NextProb;               // next probability to consider
    double delp, k;                // sampled momentum in splitting and momentum of new parton
    double pdummy = 1.;
    double dE, dpx, dpy, dpz;      // energy deposited to medium
    
    ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
    izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
    
    
    for ( int i=0; i<imax; i++)    // loop over all partons
    { 
        elastic_qq = 0;            // in the beginning assume nothing will happen
        elastic_qg = 0;            // then decide according to rates if processes occur 
        elastic_gq = 0;            // then decide according to rates if processes occur 
        elastic_gg = 0;
        radiate = 0;
        radiatePhoton = 0;
        radiate_ggg = 0;
        radiate_gqq = 0;
        convert_quark_to_gamma = 0;
        convert_quark = 0;
        convert_gluon = 0;
        
        counter++;   // counter that characterizes every step uniquely
                     // it is used to give new partons a unique color index
        if (plist->at(i).frozen() == 2) 
            continue; // if a particle is frozen out, don't evolve it any further
        
        id = plist->at(i).id();   // id of parton i 
        vecp = plist->at(i).p();  // four-vector p of parton i
        p = vecp.pAbs();          // |p| of parton i 
        
        x = plist->at(i).x();     // x value of position in [fm]
        y = plist->at(i).y();     // y value of position in [fm]
        z = plist->at(i).z();     // z value of position in [fm]
        
        // boost - using original position...
        ix = floor((hydroXmax+x)/hydroDx+0.0001);    // x-coordinate of the cell we are in now
        iy = floor((hydroXmax+y)/hydroDx+0.0001);    // y-coordinate of the cell we are in now
                                                     // note that x and y run from -hydroXmax to +hydroXmax
                                                     // and ix and iy from 0 to 2*hydroXmax hence the (hydroXmax+x or y) for both
        t = it*dtfm; // now let partons travel doing their vacuum shower stuff first...
        
        double pz_local = vecp.pz();
        double rapidity = 0.5*log( (p+pz_local)/(p-pz_local) );
        if(trackHistory and abs(rapidity)<1.0)
        {
            Tt->push_back(t);
            Tx->push_back(x);
            Ty->push_back(y);
            Tz->push_back(z);
            TE->push_back(p);
            Tpx->push_back(vecp.px());
            Tpy->push_back(vecp.py());
        }
        
        double eps = 1e-15;
        if (tauEtaCoordinates == 1) 
        {
            eta = 0.;
            if (fabs(t) < eps and fabs(z) < eps) 
                eta = 0.;
            else if ( t > 0 and fabs(z)<t ) 
                eta = 0.5 * log ( (t+z)/(t-z) );    //z is converted to eta in tau-eta coordinates
            else 
            {
                plist->at(i).frozen(2);
                continue;
            }
        }
        iz = floor((hydroZmax+eta)/hydroDz+0.0001); // z-coordinate of the cell we are in now
        // hydroZmax, hydroDz are actually hyDroEtamax, hydroDeta in tau-eta coordinates

        if (fixedTemperature == 0)
        {
            tau = 0.;           //just give tau some value - will change below
            plist->at(i).tFinal(t); // update the time - finally will be the final time
            plist->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
            // +++ if you want to move the partons before tau0 do the position update here.
            if ( moveBeforeTau0 == 1 )
            {
                plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
                plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            }
            
            if(z*z<t*t)
                tau = sqrt(t*t-z*z);    // determine tau (z*z) should always be less than (t*t)
            else 
                continue;
            
            if (tau < hydroTau0)
            {
                continue;       // do not evolve if tau < tau0 or tau > tauMax
            }
            if (tau > hydroTauMax - hydroDtau)
            {
                plist->at(i).frozen(2);
                continue;       // do not evolve if tau < tau0 or tau > tauMax
            }
            
            // +++ if you DO NOT want to move the partons before tau0 do the position update here.
            if ( moveBeforeTau0 == 0 )
            {
                plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
                plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            }
            
            // stop evolution if the position of the parton isout of range
            if ( ix < 0 || ix >= ixmax || iy < 0 || iy >= ixmax ||  iz < 0 || iz >= izmax ) 
            {
                plist->at(i).frozen(2);
                continue;
            }
            
            // get the temperature and flow velocity at the current position and time:
            hydroInfo = hydroSetup->getHydroValues(x, y, z, t);
            
            T = hydroInfo.T;
            vx = hydroInfo.vx;
            vy = hydroInfo.vy;
            vz = hydroInfo.vz;
            
            // parton stops interacting with medium if T < hydroTfinal
            // but still evolves
            if ( T < hydroTfinal ) continue;
            
            // the initial momentum vector
            px = vecp.px();
            py = vecp.py();
            pz = vecp.pz();
            
            // absolute value of flow velocity
            beta = sqrt(vx*vx+vy*vy+vz*vz);
            // angle between the particle and flow directions.
            cosPhi = (px*vx + py*vy + pz*vz)/(vecp.pAbs()*beta);
            gamma = 1./sqrt(1.-beta*beta);
            // boost p to fluid cell's rest frame
            pRest = p * gamma * (1.-beta*cosPhi);

            if(pCutPropToT)
                pCut = eLossCut*T; //standard: 4*T
            
            double eLossCut0 = eLossCut;
            // special treatment for negative particles
            if(plist->at(i).recoil() == -1) 
                eLossCut0 = 1.5;
            
            // parton stops evolution if energy is below threshold
            if(pRest/T < eLossCut0) 
                continue;
            
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
            vecpRest.e(sqrt(pxRest*pxRest + pyRest*pyRest + pzRest*pzRest));
            
            // angle between pRest and flow vel.
            cosPhiRest = (pxRest*vx + pyRest*vy + pzRest*vz)/(pRest*beta);
        }
        else // the fixed temperature case:
        {
            plist->at(i).tFinal(t); // update the time - finally will be the final time
            beta = 0.;
            gamma = 1.; 
            cosPhiRest=1.;
            pRest = p;
            T = fixedT;
            vecpRest=vecp;
            plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
            plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            plist->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
            x = plist->at(i).x();                              // x value of position in [fm]
            y = plist->at(i).y();                              // y value of position in [fm]
            z = plist->at(i).z();                              // z value of position in [fm]
            if (pRest < eLossCut) // momentum cut: to avoid sampling issues later on  RMY Jan 11 2022
            {
                plist->at(i).frozen(2);
                continue;
            }
        }
      
        boostBack = gamma*(1.+beta*cosPhiRest);  // boost factor back from rest- to lab-frame 
        
        // Sangyong's addition starts
        // see if parton i has a mother or a daughter
        in_coherence = 0; // default no
        if (formationTime)
        {
            if(plist->at(i).daughter() != -1) 
                in_coherence = hasDaughter(plist, i, T);
            if(plist->at(i).mother() != -1) 
                in_coherence = hasMother(plist, i, T);
        }
        // has... returns 1 if still coherent
        // has... returns 0 if separated. also unhooks them.
        // Sangyong's addition ends

        alpha_s_rad = alpha_s;
        alpha_s_elas = alpha_s;
        
        if(runningRad || runningElas)
        {
            double alpha_s_max = 0.42;
            double alpha_s_min = 0.15;
            double beta0 = 9./(4.*M_PI);
            double lambda_QCD = 0.20;
            double eps = 1e-1;
            double Ca = 0.;
            if (fabs(id) > 0 and fabs(id) < 4) Ca = 4./3.;
            else if (id == 21) Ca = 3.;
            double mD_sq = (2./3.)*M_PI*alpha_s*pow(T, 2.)*(2.*Nc+Nf);
            double qmax_sq = 6.*pRest*T;
            double q_hat = Ca*alpha_s*mD_sq*T*log(1+qmax_sq/mD_sq);
        
            if(runningRad)
            {
                double mu = renormFacRad*pow(q_hat*pRest, 0.25);
                if(mu < lambda_QCD) mu = lambda_QCD + eps;
                alpha_s_rad = 1./(beta0*log(pow(mu/lambda_QCD, 2.)));
                
                if(alpha_s_rad < 0.)
                {
                  cout << "[Warning]::alpha_s_rad is smaller than 0." << endl;
                  alpha_s_rad = 0.;
                }
                if(alpha_s_rad > alpha_s_max) alpha_s_rad = alpha_s_max;
            }
        
            if(runningElas)
            {
                double qmin_sq =0.05*T;
                double scat_rate = Ca*alpha_s*T*(log(qmax_sq/qmin_sq) + 
                                                 log((qmin_sq+mD_sq)/(qmax_sq+mD_sq)));
                double mfp = 1./scat_rate;
                double mu = renormFacElas*sqrt(q_hat*mfp);
                if(mu < lambda_QCD) mu = lambda_QCD + eps;
                alpha_s_elas = 1./(beta0*log(pow(mu/lambda_QCD, 2.)));
              
                if(alpha_s_elas < 0.)
                {
                  cout << "[Warning]::alpha_s_elas is smaller than 0." << endl;
                  alpha_s_elas = 0.;
                }
                if(alpha_s_elas > alpha_s_max)       alpha_s_elas = alpha_s_max;
                else if (alpha_s_elas < alpha_s_min) alpha_s_elas = alpha_s_min;
            }
        }

        // get total probabilities for radiative processes
        if (doRadiative == 1 and pRest/T > AMYpCut) 
        {
            n = rates->integrate(pRest/T,T,alpha_s_rad, Nf);
            // n stores the areas under the prob. distr.
            // Warn if probability for any process is larger than 1:
            if (dt*n.Gamma/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability to emit during one time step " << dt*n.Gamma/gamma 
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*n.Gamma_ggg/gamma >1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability ggg during one time step " << dt*n.Gamma_ggg/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << ", T=" << T << endl;
            if (dt*n.Gamma_gqq/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability gqq during one time step " << dt*n.Gamma_gqq/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*n.Gamma_em/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability qqgamma during one time step " << dt*n.Gamma_em/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
        }
      
        // get total probabilities for elastic processes

        if (doElastic == 1 ) 
        {
            qqRate = elastic->totalRate(pRest, T, alpha_s_elas, Nf, 1);
            gqRate = qqRate*9./4.;
            qgRate = elastic->totalRate(pRest, T, alpha_s_elas, Nf, 3);
            ggRate = qgRate*9./4.; 
            conversionqg = rates->Gammaqg( pRest, T, alpha_s_elas );
            conversiongq = rates->Gammagq( pRest, T, alpha_s_elas, Nf );
            if ( photonSwitch == 1 ) 
            {
                conversionqgamma = rates->Gammaqgamma( pRest, T, alpha_s );
                if ( abs(id) == 2 ) 
                    conversionqgamma*=(4./9.); // multiplying by (ef/e)^2
                else if ( abs(id) == 1 or abs(id) == 3 ) 
                    conversionqgamma*=(1./9.); // multiplying by (ef/e)^2
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
  

      if ( abs(id) > 0 and abs(id) < 4 ) // if parton is a quark (u, d, or s), let it evolve like a quark
      {
        // Check if we want a perturbative calculation. If so, do that here
        // RY June 11 2021
        if ( pertCalc ) 
        {
            // Commented out: Just create the photons and add them to the parton list
            // negative: the parton list can become too large, positive: more control at the 
            // analysis stage. Rouz, August 2021 
            // first ensure that both histograms for perturbative calculation
            // have been allocated. They should be.
            //if ( conversion_photons != NULL and amy_photons != NULL)
            //{
                
                //Parton part = plist->at(i);
                //double chargeSquared = abs(id) == 2 ? 4/9. : 1/9.;
                //double boostFactor =  gamma*(1.0 - beta*cosPhi);
                //double weightMartini;//, weightAngular;
          
                //if (part.status() > 0 ) //&& useRateTableforConversion)
                //{
                //    if(abs(part.p().eta()) < etaCut_pert)
                //    // if the particle is in the pseudorapidity window, apply collinear approximation
                //    {
                //        weightMartini = chargeSquared*rates->Gammaqgamma( pRest, T, alpha_s );
                //        conversion_photons->Fill(part.p().pT(), weightMartini*boostFactor*dt/gamma);
                //    }
                //    if(pRest/T > AMYpCut)
                //    {   
                //        Norms tempNorm = rates->integrate(pRest/T, T, alpha_s, Nf);
                //        ReturnValue dd = rates->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 4);
                //        Parton pp;
                //        pp.p(dd.y*T*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                //             dd.y*T*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
                //             dd.y*T*vecp.pz()/vecp.pAbs()*boostBack);
                //        if ( abs(pp.p().eta()) < etaCut_pert)
                //        {
                //             amy_photons->Fill(pp.p().pT(), chargeSquared*boostFactor*tempNorm.Gamma_em*dt/gamma);
                //        } 
                //    }
                //}
            //}
            //else
            //{
            //    cout<<"ERROR: Histograms for perturbative calculation are not allocated. EXIT."<<endl;
            //    exit(-1);
            //}
            //Deal with the bremsstrahlung photon first:
            int n_photons;
            if (pRest/T > AMYpCut)
            {            
                double brem_rate = rates->integrate( pRest, T, alpha_s_rad, Nf).Gamma_em*exaggerate_brem*dt/gamma;
                brem_rate = abs(id) == 2 ? brem_rate*(4./9) : brem_rate*(1./9);
                double brem_rate_compare = random->genrand64_real1();
                n_photons = floor(brem_rate);
                if (brem_rate - n_photons > brem_rate_compare)
                    n_photons += 1;
                ReturnValue d;
                double DELTAP;
                for (int j = 0; j < n_photons; j++)
                { 
                    d = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 4); 
                    DELTAP = d.y*T;
                    //cout << "\t brem photon number "<< j << " delp: "<< delp << endl;
                    if ( DELTAP < pRest and DELTAP > 0. and DELTAP > pCut)
                    {
                        //p = (pRest - DELTAP)*boostBack;
                        //k = DELTAP;
                        Parton photonAMY;
                        photonAMY.p(DELTAP*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                                    DELTAP*vecp.py()/vecp.pAbs()*boostBack,
                                    DELTAP*vecp.pz()/vecp.pAbs()*boostBack);
                        photonAMY.id(22);        // emitted parton is a photon
                        photonAMY.mass(0.);
                        photonAMY.frozen(2);     // photons do not interact anymore
                        photonAMY.x(x);          // set the new parton's initial position
                        photonAMY.y(y);
                        photonAMY.z(z);
                        photonAMY.tFinal(t);    // creation time of this photon
                        photonAMY.tauAtEmission(tau); //creation (proper) time of the photon
                        photonAMY.source(2);     // AMY photon
                        photonAMY.daughter_of(-1); //don't update mother and daughter information
                        photonAMY.mother_of(-1);  // this is a perturbative calculation, piggy-backing off the main MC engine
                        photonAMY.recoil(plist->at(i).recoil());//inherit the recoil code of mother
                        plist->push_back(photonAMY); // add the photon to the list of partons
                    }
                }
            }
            //cout<<"num brem photons : "<<n_photons<<endl;
            //Deal with conversion photon next:
            double convrate = rates->Gammaqgamma(pRest, T, alpha_s_elas)*exaggerate_conv*dt/gamma;
            convrate = abs(id) == 2 ? convrate*(4./9) : convrate*(1./9);
            double conv_rate_compare = random->genrand64_real1();
            n_photons = floor(convrate);
            if (convrate - n_photons > conv_rate_compare)
                n_photons += 1;
            //cout << "num conv photons : "<<n_photons<<endl;
            for ( int j = 0 ; j < n_photons; j++)
            {
                Parton photonConv;
                photonConv.col(0);
                photonConv.acol(0);
                photonConv.id(22);
                photonConv.mass(0.);
                photonConv.frozen(2); 
                photonConv.x(x);
                photonConv.y(y);
                photonConv.z(z);
                photonConv.tFinal(t); // creation time of this photon
                photonConv.tauAtEmission(tau);
                photonConv.p(px, py, pz);
                photonConv.daughter_of(-1);
                photonConv.mother_of(-1);
                photonConv.source(1);//conversion photon
                photonConv.recoil(plist->at(i).recoil()); // inherit the recoil code of mother
                plist->push_back(photonConv);
            }
            // Done with the perturbative stuff.
        }    
          // in the following decide what to do in this time step
          // Here begin the modifications that will allow elastic collisions below the 
          // AMYpCut, while not sampling radiative splittings below this scale. -CFY 9/28/2011
          if (doRadiative == 1 and doElastic == 1 ) //radiative + elastic 
          {
              if ( abs(id) == 2 ) // multiply by (e_f/e)^2
                  n.Gamma_em*=4./9.;
              else
                  n.Gamma_em*=1./9.;

              if(pRest/T>AMYpCut)
              {
                  totalQuarkProb = dt/gamma*(n.Gamma)
                                 + dt/gamma*(qqRate)+dt/gamma*(qgRate)
                                 + dt/gamma*conversionqg+dt/gamma*conversionqgamma;
              }
              else
              {
                  totalQuarkProb = dt/gamma*(qqRate)+dt/gamma*(qgRate)
                                 + dt/gamma*conversionqg+dt/gamma*conversionqgamma;
              }

              if ( photonSwitch == 1 && pRest/T>AMYpCut) 
              {
                  totalQuarkProb += dt/gamma*(n.Gamma_em);
                  EMProb = (dt/gamma*(n.Gamma_em))/totalQuarkProb;
              }
              else
              {
                  EMProb = 0.0;
              }

              if (totalQuarkProb>1)
              {
                  cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                       << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
                  cout << "Parton with pRest : "<<pRest<<" and source : "<<plist->at(i).source()<<endl;
              }
        
              if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
              {
                  
                  if(pRest/T>AMYpCut)
                  {
                      // Sangyong's addition:
                      // in_coherence is the new flag
                      // added if(in_coherence==0) to radiate* and convert*
                      // this includes coherence effects in radiations, but not in conversions
                      // aways do elastic - elastic was lumped together with the conversions
                              // I am modifying this to use one random number
                      // we either have
                              // totalQuarkProb = dt/gamma*(n.Gamma)
                      //  +dt/gamma*(qqRate)+dt/gamma*(qgRate)
                      //  +dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                      // or
                              // totalQuarkProb = dt/gamma*(n.Gamma)
                      //  +dt/gamma*(qqRate)+dt/gamma*(qgRate)
                      //  +dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                      //  +dt/gamma*(n.Gamma_em);
                      
                      randx = random->genrand64_real1();
                      
                      AccProb = 0.0;
                      NextProb = (dt/gamma*(qqRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qq = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(qgRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qg = 1;
                      
                      AccProb += NextProb;
                      NextProb = dt/gamma*(n.Gamma)/totalQuarkProb;
                      if(AccProb <= randx && randx < (AccProb + NextProb))
                      {
                          if(in_coherence==0) 
                              radiate = 1;
                      }

                      AccProb += NextProb;
                      NextProb = (dt/gamma*conversionqg)/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark = 1;
                      }

                      // +dt/gamma*conversionqgamma; 
                      AccProb += NextProb;
                      NextProb = (dt/gamma*conversionqgamma)/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                             convert_quark_to_gamma = 1;
                      }

                      // EMProb = dt/gamma*(n.Gamma_em)/totalQuarkProb
                      // this is set to zero when photonSwitch == 0
                      AccProb += NextProb;
                      NextProb = EMProb;
                      if( AccProb <= randx && randx <= (AccProb + NextProb) )
                      {
                          if(in_coherence==0) 
                              radiatePhoton = 1;
                      }
                  }// if(pRest/T>AMYpCut)
                  else
                  {
		              randx = random->genrand64_real1();

                      AccProb = 0.0;
                      NextProb = (dt/gamma*(qqRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qq = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(qgRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qg = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(conversionqg))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark = 1;
                      }

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(conversionqgamma))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark_to_gamma = 1;
                      }
                  }// pRest/T < AMYpCut
              }// something happens to the quark
              else  // nothing happens to the quark
              {
                  continue;
              }// nothing happens
          }// rad + elastic
          //Here ends the modification of the rates for quarks, on to gluons. -CFY
          else if (doRadiative == 1 and doElastic == 0 and pRest/T>AMYpCut) //radiative only
          {
              // Sangyong's addition:
              // in_coherence is the new flag defined in Parton.h
              // originally, I added if(in_coherence==0) to radiate*
              // but since this part is radiative only, that does not make sense. removed.
              totalQuarkProb = dt/gamma*(n.Gamma);
              if ( photonSwitch == 1 ) totalQuarkProb += dt/gamma*(n.Gamma_em);

              if (totalQuarkProb>1)
                  cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                       << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
              if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
              {
                  if ( random->genrand64_real1() < dt/gamma*(n.Gamma)/totalQuarkProb )
                      radiate = 1; 
                  else if ( photonSwitch == 1 )
                      radiatePhoton = 1;
              }
              else // nothing happens to the quark
              {
                  continue;
              }
          } // rad only
          else if (doRadiative == 0 and doElastic == 1) //elastic only
          // no sense in implementing coherence effect here, either
          {
                totalQuarkProb = dt/gamma*(qqRate)+dt/gamma*(qgRate)
		               + dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                // dt/gamma is dt_rest-frame
                
                if (totalQuarkProb>1)
                cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                     << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
                if ( random->genrand64_real1() < totalQuarkProb )  // check if something happens
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
                else 
                {
                    continue;
                }
          }
          else
              continue;

          // now do what has been decided before
          if( radiate == 1 ) // see if emission happens
          {
              if(pRest/T>AMYpCut) // do not evolve partons with momenta below this scale 
              {
                  // do process 1, q->qg
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 1); 
                  delp = f.y*T;
                  if(delp > pRest) continue;

                  p = (pRest - delp)*boostBack;    // quark's new momentum in the lab frame
                  newvecp=p/vecp.pAbs()*vecp;      // new quark momentum
                  plist->at(i).p(newvecp);         // change quark's momentum to the new one

                  if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }

                  if (delp > pCut) // if new parton is kept (f.y*T > threshold [in GeV])
                  {
                      newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack,  // emitted gluon's momentum
                               delp*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
                               delp*vecp.pz()/vecp.pAbs()*boostBack);
                      col = plist->at(i).col();              // color of original quark  
                      acol = plist->at(i).acol();            // anti-color
                      plist->at(i).increment_hard_radiation(); // Rouz : count hard radiations
                      if (col!=0)                               // if we had a quark
                      {
                          plist->at(i).col(100000000+counter);  // set color to new color
                          newOne.col(col);                      // set new gluon's color to quark's original color
                          newOne.acol(100000000+counter);       // set new gluon's anti-color to quark's new color
                      }
                      else if (acol!=0) // if we had an anti-quark
                      {
                          plist->at(i).acol(100000000+counter);  // set anti-color to new color
                          newOne.col(100000000+counter); // set new particle's color
                          newOne.acol(acol);            // set new particle's anti-color   
                      }

                      newOne.id(21);    // emitted parton is a gluon
                      newOne.init_counts();//Rouz : initialize the counts
                      newOne.mass(0.);
                      newOne.frozen(0);
                      newOne.x(x);      // set the new parton's initial position
                      newOne.y(y);
                      newOne.z(z);
                      newOne.source(11); // AMY parton, gluon emitted from quark or anti-quark
                      newOne.recoil(plist->at(i).recoil());
		      
                      //Sangyong's addition for splits
                      daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                      newOne.daughter_of(i); // daughter of the current parton
                      newOne.mother_of(-1);  // mother of nobody
                      newOne.p_at_split(newOne.p()); // momentum at the split point
                      
                      plist->at(i).daughter_of(-1); // daughter of nobody
                      plist->at(i).mother_of(daughter);
                      plist->at(i).p_at_split(plist->at(i).p());
                      //Sangyong's addition end
                      plist->push_back(newOne);  // add the gluon to the list of partons if (f.y>AMYpCut)
                  }
                  else
                  {
                    // delp < pCut: soft radiation.
                    plist->at(i).increment_soft_radiation();
                  }
                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          else if( radiatePhoton == 1 ) // see if photon emission happens
          {
              if(pRest/T>0.01) // do not evolve partons with momenta below this scale
              {
                  // do process 4, q->qgamma
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 4); 
                  delp = f.y*T;

                  if(delp > pRest || delp < 0.) continue;           // do nothing at this moment
                  else if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }
                  p = (pRest - delp)*boostBack;
                  k = delp;
                  newvecp = p/vecp.pAbs()*vecp;               // four vector for new quark momentum
                  plist->at(i).p(newvecp);                 // update new quark momentum
                   
                  newOne.p(k*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                           k*vecp.py()/vecp.pAbs()*boostBack,
                           k*vecp.pz()/vecp.pAbs()*boostBack);
                  newOne.id(22);                 // emitted parton is a photon
                  newOne.init_counts(); // Rouz's addition. New parton has had no scatterings so far
                  //newOne.splits(0);
                  newOne.mass(0.);
                  newOne.frozen(2);              // photons do not interact anymore
                  newOne.x(x);                   // set the new parton's initial position
                  newOne.y(y);
                  newOne.z(z);
                  newOne.tauAtEmission(tau);
                  newOne.source(2);              // AMY photon
                  newOne.recoil(plist->at(i).recoil());

                  //Sangyong's addition for splits
                  daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                  newOne.daughter_of(i); // daughter of the current parton
                  newOne.mother_of(-1);  // mother of nobody
                  newOne.p_at_split(newOne.p()); // momentum at the split point
                  
                  plist->at(i).daughter_of(-1); // daughter of nobody
                  plist->at(i).mother_of(daughter);
                  plist->at(i).p_at_split(plist->at(i).p());
                  //Sangyong's addition

                  plist->push_back(newOne); // add the photon to the list of partons
                  if (delp < pCut) // Rouz : update count of soft or hard radiations
                    plist->at(i).increment_soft_radiation();
                  else
                    plist->at(i).increment_hard_radiation();

                  if(outputSource and delp > pRest - pCut)
                  {
                      // if delp > pRest, jet is converted to photon
                      // and gains energy from thermal medium
                      // dE is negative but the direction is same as mother jet.
                      // if pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      dE = pRest-delp;
                      dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                      dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                      dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();

                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);
                      slist->push_back(newSource);
                  }
              }
          }                             
          else if ( convert_quark == 1  ) // do q->g
          {
              col = plist->at(i).col();   // color of original quark  
              acol = plist->at(i).acol(); // anti-color
              if (col!=0)                 // if we had a quark
              {
                  plist->at(i).acol(200000000+counter);  // add a new anti-color for the gluon
                  newOne.col(200000000+counter);         // create new thermal gluon
                  newOne.acol(300000000+counter);                        
              }
              else if (acol!=0)
              {
                  plist->at(i).col(200000000+counter);   // add a new color for the gluon
                  newOne.acol(200000000+counter);        // create new thermal gluon
                  newOne.col(300000000+counter);                        
              }
              plist->at(i).id(21);    // convert to gluon
              plist->at(i).source(14); // Conversion gluon
              double choice = random->genrand64_real1();
              if ( choice < 0.5 )
              {
                newOne.id(21);          // add a thermal gluon
                newOne.mass(0.); 
                newOne.frozen(2);
                newOne.x(x);            // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, -1));  // sample momentum from Bose distribution at temperature T
                newOne.recoil(plist->at(i).recoil());

                //Sangyong's addition for formation time of radiation
                newOne.daughter_of(-1); // daughter of nobody
                newOne.mother_of(-1);  // mother of nobody
                newOne.p_at_split(newOne.p()); // momentum at the conversion point

                plist->push_back(newOne); 
              }
              else 
              {
                 newOne.id(id);                        // add a thermal quark of same flavor
                 if (col!=0)                               // if we had a quark
                 {
                     newOne.col(300000000+counter);            // attach the color string to the thermal quark
                     newOne.acol(0);                        
                 }
                 else if (acol!=0)                         // if we had an anti-quark
                 {
                     newOne.acol(300000000+counter);           // attach the color string to the thermal anti-quark
                     newOne.col(0);                        
                 }
                 if(abs(id) < 3) 
                     newOne.mass(0.33);
                 else 
                     newOne.mass(0.5);

                 newOne.frozen(2);
                 newOne.x(x);                              // set the new parton's initial position
                 newOne.y(y);
                 newOne.z(z);
                 newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
                 newOne.recoil(plist->at(i).recoil());
                 
                 //Sangyong's addition for formation time of radiation
                 newOne.daughter_of(-1); // daughter of nobody
                 newOne.mother_of(-1);  // mother of nobody
                 newOne.p_at_split(newOne.p()); // momentum at the conversion point
                 
                 plist->push_back(newOne); 
              }
              continue;                                 //prevent the new gluon from interacting again in this time step 
          }
          else if ( convert_quark_to_gamma == 1 and pRest/T>0.01 ) // do q->gamma
          {
              col = plist->at(i).col();  // color of original quark  
              acol = plist->at(i).acol();// anti-color
              plist->at(i).acol(0);      // make color neutral photon
              plist->at(i).col(0);       // make color neutral photon
              if (col!=0)                // if we had a quark
              {
                  newOne.col(col);                      // attach the color string to the thermal gluon
                  newOne.acol(400000000+counter);                        
              }
              else if (acol!=0)
              {
                  newOne.acol(acol);                    // attach the color string to the thermal gluon
                  newOne.col(400000000+counter);                        
              }

              plist->at(i).id(22);    // convert to photon
              plist->at(i).source(1); // Conversion photon "source" code.
              plist->at(i).init_counts(); // Rouz: this sets all these counts to zero for the photon. they don't interact in MARTINI
              //TODO: here, why not a coin flip for thermal q/qbar vs gluon? Rouz Dec 8, 2021
              newOne.id(21);          // add a thermal gluon
              newOne.init_counts(); // Rouz: new parton, hasn't had any interactions yet
              newOne.mass(0.); 
              newOne.frozen(0);
              newOne.x(x);            // set the new parton's initial position
              newOne.y(y);
              newOne.z(z);
              newOne.tauAtEmission(tau);
              newOne.p(random->thermal(T, -1)); // sample momentum from Bose distribution at temperature T
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point

              plist->push_back(newOne); 

              newOne.id(id);  // add a thermal quark of same flavor as initial quark
              if (col!=0)     // if we had a quark
              {
                  newOne.col(400000000+counter); // attach the color string to the thermal quark
                  newOne.acol(0);                        
              }
              else if (acol!=0) // if we had an anti-quark
              {
                  newOne.acol(400000000+counter); // attach the color string to the thermal anti-quark
                  newOne.col(0);                        
              }
              if ( abs(id) < 3) 
                  newOne.mass(0.33);
              else 
                  newOne.mass(0.5);

              newOne.frozen(0);
              newOne.x(x); // set the new parton's initial position
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, 1)); // sample momentum from Fermi distribution at temperature T
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 
              continue;                                 //prevent the new gluon from interacting again in this time step 
          }
          else if ( elastic_qq == 1 ) // do qq->qq
          {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 1); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }
              
              if( transferTransverseMomentum == 1 and omega < 800. )//&& fabs(eta) < 4.
              {
                  int repeat = 1;
                  // k_min is a mimumum momenum of a thermal parton to be sample 
                  // for a recoil parton to be on-shell
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 1);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else 
              {
                  q = omega;
              }
          }
          else if ( elastic_qg == 1 ) // do qg->qg
          {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 3); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 && omega < 800. ) // && fabs(eta) < 4.)
              {
                  int repeat = 1;
                  // k_min is a mimumum momenum of a thermal parton to be sample 
                  // for a recoil parton to be on-shell
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 3);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else 
              {
                  q = omega;
              }
          }
          if ( elastic_qq == 1 || elastic_qg == 1 ) // in both cases ( qq->qq and qg->qg ) change the momentum
          {
              newvecpRest=elastic->getNewMomentum(vecpRest, omega*T, q*T, random); // takes everything in GeV
              if (beta > 1e-15) // before: beta == 0, should avoid equality check for doubles
              {
                  newpx = vx*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecpRest.px() 
                    + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.py()
                    + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.pz();
                  newpy = vy*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecpRest.py()
                    + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.px()
                    + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.pz();
                  newpz = vz*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecpRest.pz()
                    + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.px()
                    + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.py();
                  
                  // momentum four-vector in lab frame
                  newvecp.px(newpx);
                  newvecp.py(newpy);
                  newvecp.pz(newpz);
                  newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
              }
              else
                newvecp = newvecpRest;
              double newpRest = newvecpRest.pAbs();
              if(newpRest < pCut)
              {
                  if(outputSource)
                  {
                      // output energy deposition from jet particle absorption
                      dE = newvecpRest.pAbs();
                      dpx = newvecpRest.px();
                      dpy = newvecpRest.py();
                      dpz = newvecpRest.pz();

                      newSource.type(1);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
                  plist->at(i).frozen(2);
                  //newvecp *= pdummy/newvecp.pAbs();
                  // exclude soft recoil by assigning id == 0
                  if(plist->at(i).recoil() == -1) plist->at(i).id(0);
              }
    
              plist->at(i).p(newvecp);   // change quark's momentum to the new one

              vecq = vecpRest - newvecpRest; // momentum transfer q
              if (vecq.pAbs() > pCut) // Rouz : Count elastic scatterings
                plist->at(i).increment_hard_elastic();
              else
                plist->at(i).increment_soft_elastic();

              if(getRecoil and omega > 0. and plist->at(i).recoil() != -1)
              {
                  if( elastic_qq == 1 )   // quark recoil
                  {
                      // thermal quark to be recoiled
                      vecThermal = elastic->getThermalMomentum(vecq, T, 1, random);
                      vecRecoilRest = vecq + vecThermal;

                      pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                      pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                      pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                      
                      // momentum four-vector in lab frame
                      vecRecoil.px(pxRecoil);
                      vecRecoil.py(pyRecoil);
                      vecRecoil.pz(pzRecoil);
                      vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                       pyRecoil*pyRecoil +
                                       pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                          double r = random->genrand64_real1();
                          // note that when using general N_f this has to be changed
                          if (r<0.33) newID=1;
                          else if (r<0.66) newID=2;
                          else newID=3;
                          double mass;
                          if (newID<3) mass=0.33;  // set the quark's mass
                          else mass = 0.5;

                          newOne.p(vecRecoil);
                          newOne.mass(mass);

                          double matter = random->genrand64_real1();
                          if(matter < 0.5)   // thermal parton is quark
                          {
                              newOne.id(newID);
                              newOne.col(500000000+counter);
                              newOne.acol(0);
                              newHole.id(newID);
			                  newHole.col(600000000+counter);
			                  newHole.acol(0);
                          }
                          else  // thermal parton is anti-quark
                          {
                              newOne.id(-newID);
                              newOne.col(0);
                              newOne.acol(500000000+counter);
                              newHole.id(-newID);
			                  newHole.col(0);
			                  newHole.acol(600000000+counter);
                          }

                          newOne.x(x);
                          newOne.y(y);
                          newOne.z(z);
                          newOne.frozen(0);
                          newOne.recoil(1);
                          newOne.init_counts(); // Rouz: initialize the counts
                          //Sangyong's addition for formation time of radiation
                          newOne.daughter_of(-1); // daughter of nobody
                          newOne.mother_of(-1);  // mother of nobody
                          newOne.p_at_split(newOne.p()); // momentum at the conversion point

                          plist->push_back(newOne);    // add recoil parton to the parton list

                          pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                          pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                          pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();

                          // momentum four-vector in lab frame
                          vecThermalLab.px(pxThermalLab);
                          vecThermalLab.py(pyThermalLab);
                          vecThermalLab.pz(pzThermalLab);
                          vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                pyThermalLab*pyThermalLab +
                                                pzThermalLab*pzThermalLab));

                          newHole.p(vecThermalLab);
                          newHole.mass(mass);

                          newHole.x(x);
                          newHole.y(y);
                          newHole.z(z);
                          newHole.frozen(0);
                          newHole.recoil(-1);

                          //Sangyong's addition for formation time of radiation
                          newHole.daughter_of(-1); // daughter of nobody
                          newHole.mother_of(-1);  // mother of nobody
                          newHole.p_at_split(newHole.p()); // momentum at the conversion point

                          plist->push_back(newHole);    // add hole to the parton list
                      }
                      else
                      {
                          if(outputSource)
                          {
                              // momentum transfer q goes into medium
                              dE = vecq.e();
                              dpx = vecq.px();
                              dpy = vecq.py();
                              dpz = vecq.pz();

                              newSource.type(1);
                              newSource.tau(tau);
                              newSource.x(x);
                              newSource.y(y);
                              newSource.eta(eta);
                              newSource.dE(dE);
                              newSource.dpx(dpx);
                              newSource.dpy(dpy);
                              newSource.dpz(dpz);
                              newSource.vx(vx);
                              newSource.vy(vy);
                              newSource.vz(vz);

                              slist->push_back(newSource);
                          }
                      }
                  }
                  else if( elastic_qg == 1 )   // gluon recoil
                  {
                      // thermal gluon to be recoiled
                      vecThermal = elastic->getThermalMomentum(vecq, T, -1, random);
                      vecRecoilRest = vecq + vecThermal;

                      pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                      pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                      pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                      
                      // momentum four-vector in lab frame
                      vecRecoil.px(pxRecoil);
                      vecRecoil.py(pyRecoil);
                      vecRecoil.pz(pzRecoil);
                      vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                       pyRecoil*pyRecoil +
                                       pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                          newOne.p(vecRecoil);
                          newOne.id(21);
                          newOne.mass(0.);
                          newOne.col(700000000+counter);
                          newOne.acol(800000000+counter);

                          newOne.x(x);
                          newOne.y(y);
                          newOne.z(z);
                          newOne.frozen(0);
                          newOne.recoil(1);
                          newOne.init_counts();//Rouz : initialize counts
                          //Sangyong's addition for formation time of radiation
                          newOne.daughter_of(-1); // daughter of nobody
                          newOne.mother_of(-1);  // mother of nobody
                          newOne.p_at_split(newOne.p());

                          plist->push_back(newOne);    // add recoil gluon to the parton list


                          pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                          pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                          pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();

                          // momentum four-vector in lab frame
                          vecThermalLab.px(pxThermalLab);
                          vecThermalLab.py(pyThermalLab);
                          vecThermalLab.pz(pzThermalLab);
                          vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                pyThermalLab*pyThermalLab +
                                                pzThermalLab*pzThermalLab));

                          newHole.p(vecThermalLab);
                          newHole.id(21);
                          newHole.mass(0.);
                          newHole.col(900000000+counter);
                          newHole.acol(1000000000+counter);

                          newHole.x(x);
                          newHole.y(y);
                          newHole.z(z);
                          newHole.frozen(0);
                          newHole.recoil(-1);

                          //Sangyong's addition for formation time of radiation
                          newHole.daughter_of(-1); // daughter of nobody
                          newHole.mother_of(-1);  // mother of nobody
                          newHole.p_at_split(newOne.p());

                          plist->push_back(newHole);    // add hole to the parton list
                      }
                      else
                      {
                          if(outputSource)
                          {
                              // momentum transfer q goes into medium
                              dE = vecq.e();
                              dpx = vecq.px();
                              dpy = vecq.py();
                              dpz = vecq.pz();

                              newSource.type(1);
                              newSource.tau(tau);
                              newSource.x(x);
                              newSource.y(y);
                              newSource.eta(eta);
                              newSource.dE(dE);
                              newSource.dpx(dpx);
                              newSource.dpy(dpy);
                              newSource.dpz(dpz);
                              newSource.vx(vx);
                              newSource.vy(vy);
                              newSource.vz(vz);

                              slist->push_back(newSource);
                          }
                      }
                  }
              }
              else
              {
                  if(outputSource)
                  {
                      // momentum transfer q goes into medium
                      dE = vecq.e();
                      dpx = vecq.px();
                      dpy = vecq.py();
                      dpz = vecq.pz();

                      newSource.type(1);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
      }
      if ( id == 21 ) // if parton is a gluon, let it evolve like a gluon
      {
          // in the following decide what to do in this time step
          // Again, modified for allowing elastic collisions at momenta lower than AMYpCut. -CFY
          // Sangyong's addition:
          // in_coherence is the new flag
          // added if(in_coherence==0) to radiate*
          if (doRadiative == 1 and doElastic == 1 )  //radiative + elastic
              {
                  if(pRest/T>AMYpCut)
                  {
                      totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq)
                                     + dt/gamma*(gqRate)+dt/gamma*(ggRate)
                                     + dt/gamma*conversiongq; 
                  }
                  else
                  {
                      totalGluonProb = dt/gamma*(gqRate)+dt/gamma*(ggRate)
			             + dt/gamma*conversiongq;
                  }
                  if (totalGluonProb>1)
                  {
                      cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
                          << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
                      cout << " pRest: "<< pRest<< " source : "<<plist->at(i).source()<<endl;
                  }
                  if ( random->genrand64_real1() < totalGluonProb ) // check if something happens with the quark
                  {
                      if(pRest/T>AMYpCut)
                      {
                          // Sangyong's addition:
                          // in_coherence is the new flag
                          // added if(in_coherence==0) to radiate* and convert*
                          // this includes coherence effects in radiations, but not in conversions
                          // aways do elastic - elastic was lumped together with the conversions
                          // I am modifying this to use one random number
                          // we have
                          // totalGluonProb = 
                          // dt/gamma*(n.Gamma_ggg)
                          // +dt/gamma*(n.Gamma_gqq)
                          // +dt/gamma*(gqRate)
                          // +dt/gamma*(ggRate)
                          // +dt/gamma*conversiongq; 
                          
                          randx = random->genrand64_real1();
                          
                          AccProb = 0.0;
                          NextProb = (dt/gamma*(gqRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gq = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(ggRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gg = 1;

                          AccProb += NextProb;
                          NextProb = dt/gamma*(n.Gamma_ggg)/totalGluonProb;
                          if(AccProb <= randx && randx < (AccProb + NextProb))
                          {
                              if(in_coherence==0) 
                                  radiate_ggg = 1;
                          }
                          
                          AccProb += NextProb;
                          NextProb = dt/gamma*(n.Gamma_gqq)/totalGluonProb;
                          if(AccProb <= randx && randx < (AccProb + NextProb))
                          {
                              if(in_coherence==0) 
                                  radiate_gqq = 1;
                          }
                          
                          AccProb += NextProb;
                          NextProb = (dt/gamma*conversiongq)/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                          {
                              if(in_coherence==0 or plist->at(i).daughter() != -1)
                                  convert_gluon = 1;
                          }
                      }// pRest/T > AMYpCut
                      else
                      {
		                  randx = random->genrand64_real1();

                          AccProb = 0.;
                          NextProb = (dt/gamma*(gqRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gq = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(ggRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gg = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(conversiongq))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                          {
                              if(in_coherence==0 or plist->at(i).daughter() != -1)
                                  convert_gluon = 1;
                          }
                      }
                  }
                  else 
                  {
                      continue;
                  }
              }
              else if (doRadiative == 1 and doElastic == 0 and pRest/T>AMYpCut) //radiative
            // no sense implementing coherence here
              {
                  totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq); 
                  if (totalGluonProb>1) 
                      cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
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
                      continue;
                  }
              }// rad only
              else if (doRadiative == 0 and doElastic == 1 ) //elastic
              // no sense implementing coherence here
              {
                  totalGluonProb = dt/gamma*(gqRate)+dt/gamma*(ggRate)
		                 + dt/gamma*conversiongq; // dt/gamma is dt_rest-frame

                  if (totalGluonProb>1) 
                      cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
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
                              elastic_gg = 0;
                          }
                          else
                          {
                              elastic_gg = 1;
                              elastic_gq = 0;
                          }
                      }
                  }
                  else 
                  {
                      continue;
                  }
              }
            else 
                  continue;
          
          // now do what has been decided before
          if(radiate_ggg == 1) // g -> gg
          {
              if(pRest/T>AMYpCut) 
              {
                  f = rates->findValuesRejection
                    (pRest/T, T, alpha_s_rad, Nf, random, import, 3);   // do process 3

                  delp = f.y*T;
                  if(delp > pRest) continue;

                  p = (pRest - delp)*boostBack;
                  newvecp=p/vecp.pAbs()*vecp; // new gluon momentum
                  plist->at(i).p(newvecp);    // change quark's momentum to the new one

                  if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }
                  if(delp > pCut)
                  {
                      newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack, // emitted gluon's momentum (coll) 
              	               delp*vecp.py()/vecp.pAbs()*boostBack,
              	               delp*vecp.pz()/vecp.pAbs()*boostBack); 
                      acol = plist->at(i).acol();            // the gluon's anti-color
                      plist->at(i).increment_hard_radiation(); // Rouz: increment hard radiation
                      plist->at(i).acol(1100000000+counter); // set the first new gluon's anti-color to a new one 
                      newOne.id(21);                         // second new particle is a gluon too
                      newOne.col(1100000000+counter);        // set the second gluon's color to the new color 
                      newOne.acol(acol);                     // set second gluon's anti-color to the original one
                      newOne.x(x);
                      newOne.y(y);
                      newOne.z(z);
                      newOne.mass(0.);  // set second gluon's mass to zero
                      newOne.frozen(0);
                      newOne.source(12); //  RY: gluon emitted from a gluon. All photons produced in evolve(...) have itsSource=1
                      newOne.recoil(plist->at(i).recoil());
                      newOne.init_counts();//Rouz:: start the counts for this parton
                      
                      //Sangyong's addition for formation time of radiation
                      daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                      newOne.daughter_of(i); // daughter of the current parton
                      newOne.mother_of(-1);  // mother of nobody
                      newOne.p_at_split(newOne.p()); // momentum at the split point
              
                      plist->at(i).daughter_of(-1); // daughter of nobody
                      plist->at(i).mother_of(daughter);
                      plist->at(i).p_at_split(plist->at(i).p());
                      
                      plist->push_back(newOne);                  // add the second gluon to the list of partons
                  }
                  if ( delp < pCut)
                    plist->at(i).increment_soft_radiation();

                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          //pRest = p * gamma * (1.-beta*cosPhi); // boost p to fluid cell's rest frame
          if( radiate_gqq == 1 ) // g -> qq
          {
	          if(pRest/T>AMYpCut) 
              {
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 2); // do process 2
                  delp = f.y*T;
                  if(delp > pRest or delp < 0.)
                      continue;

                  p = (pRest - delp)*boostBack;
                  // choose if it is a u-ubar, d-dbar or s-sbar pair:
                  double r = random->genrand64_real1();
                  // note that when using general N_f this has to be changed
                  if ( r < 0.33 ) newID=1;
                  else if  ( r < 0.66 ) newID=2;
                  else newID=3;
                  double mass;
                  if ( newID < 3 ) mass=0.33;  // set the quark's mass
                  else mass = 0.5;         // (pythia needs that for fragmentation) 
                  
                  col = plist->at(i).col();   // gluon's color
                  acol = plist->at(i).acol(); // gluon's anti-color 
                  int col_m, acol_m, col_d, acol_d;
                  double matter = random->genrand64_real1();
                  if (matter<0.5) //flip matter and antimatter
                  {
                      newID *= -1;
                      col_m  = 0;
                      acol_m = acol;
                      col_d  = col;
                      acol_d = 0;
                  }
                  else
                  {
                      col_m  = col;
                      acol_m = 0;
                      col_d  = 0;
                      acol_d = acol;
                  }
                  
                  newvecp = p/vecp.pAbs()*vecp; // new quark momentum
                  plist->at(i).id(newID);     // turn gluon into quarks
                  plist->at(i).col(col_m);    // set quark's color
                  plist->at(i).acol(acol_m);  // set quark's anti-color to zero
                  plist->at(i).mass(mass);    // set quark's mass
                  plist->at(i).p(newvecp);    // change quark's momentum to the new one
                  plist->at(i).source(13);    // this fermion now is the result of gluon splitting -RY.
                  if(p < pCut)
                  {
                       plist->at(i).frozen(2);
                       // exclude soft recoil by assigning id == 0
                       if(plist->at(i).recoil() == -1)
                          plist->at(i).id(0);
                  }
                  plist->at(i).init_counts(); // Rouz: initialize counts for both since 
                  newOne.init_counts();       // this is sort of a conversion process
                  newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack,   // second quark's momentum (collinear)
                	       delp*vecp.py()/vecp.pAbs()*boostBack,
                	       delp*vecp.pz()/vecp.pAbs()*boostBack); 
                  newOne.id(-newID);                                // second new particle is an anti-quark (minus-sign)
                  newOne.mass(mass);                                // set anti-quark's mass  
                  newOne.col(col_d);                                    // set anti-quark's color to zero
                  newOne.acol(acol_d);                                // set anti-quark's anti-color to gluon's anti-color
                  newOne.x(x);
                  newOne.y(y);
                  newOne.z(z);
                  if(delp > pCut)
                      newOne.frozen(0);
                  else
                      newOne.frozen(2);
                  newOne.source(13); // RY: this fermion is also the result of gluon splitting 
                  newOne.recoil(plist->at(i).recoil());

                  //Sangyong's addition for formation time of radiation
                  daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                  newOne.daughter_of(i);    // daughter of the current parton
                  newOne.mother_of(-1);     // mother of nobody
                  newOne.p_at_split(newOne.p()); // momentum at the split point

                  plist->at(i).daughter_of(-1); // daughter of nobody
                  plist->at(i).mother_of(daughter);
                  plist->at(i).p_at_split(plist->at(i).p());

                  plist->push_back(newOne); // add the second quark to the list of partons
                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          if ( convert_gluon == 1 ) // do g->q
          {
              if ( random->genrand64_real1() < 0.5 ) 
                GluonToQuark = 1;
              else 
                GluonToQuark = 0;
         
              double r = random->genrand64_real1();
              if (r<0.33) newID=1;
              else if  (r<0.66) newID=2;
              else newID=3;
              
              col = plist->at(i).col();   // color of original gluon  
              acol = plist->at(i).acol(); // anti-color

              if (GluonToQuark == 1)
              {
                plist->at(i).acol(0);     // make the gluon a quark
                newOne.acol(acol);        // let the thermal gluon have the gluon's previous anti-color
                newOne.col(1200000000+counter);
                plist->at(i).id(newID);   // convert to quark
              }
              else
              {
                plist->at(i).col(0); // make the gluon an anti-quark
                newOne.col(col);     // let the thermal gluon have the gluon's previous color
                newOne.acol(1200000000+counter);
                plist->at(i).id(-newID);  // convert to anti-quark
              }
              if (newID == 1 or newID == 2) 
                plist->at(i).mass(0.33);           
              else
                plist->at(i).mass(0.55);
              plist->at(i).source(15);
              plist->at(i).init_counts();// Rouz: initialize counts 
              newOne.init_counts();      // Rouz: initialize counts
              newOne.id(21); // add a thermal gluon
              newOne.mass(0.); 
              newOne.frozen(2);
              newOne.x(x);
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, -1)); // sample momentum from Bose distribution at temperature T
              newOne.source(1);                 // 
              newOne.recoil(plist->at(i).recoil());
              
              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 

              double mass;
              if (newID<3) mass=0.33;  // set the quark's mass
              else mass = 0.5;         // (pythia needs that for fragmentation) 

              if (GluonToQuark == 1)
              {
                newOne.acol(1200000000+counter); // let the thermal anti-quark have the right anti-color
                newOne.col(0);
                newOne.id(-newID);              // add a thermal anti-quark
              }
              else
              {
                newOne.col(1200000000+counter);  // let the thermal quark have the right color
                newOne.acol(0);
                newOne.id(newID);  // add a thermal quark
              }

              newOne.mass(mass); 
              newOne.frozen(0);
              newOne.x(x);
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, 1)); // sample momentum from Fermi distribution at temperature T
              newOne.source(1);               // RY: gluon to fermion conversion. 
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 

              continue;  //prevent the new gluon from interacting again in this time step 
            }
            else if ( elastic_gq == 1 ) // do gq->gq
            {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 2); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 and omega < 800. )//&& fabs(eta) < 4.)
              {
                  int repeat = 1;
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 2);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else
              {
                  q = omega;
              }
            }
            else if ( elastic_gg == 1 ) // do gg->gg
            {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 4); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 && omega < 800. )//&& fabs(eta) < 4.)
                {
                    int repeat = 1;
                    double k_min;
                    do
                    {
                        q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 4);
                        k_min = (q-omega)/2.;
                        if(k_min < 20.) repeat = 0;
                    } while (repeat);
                }
                else
                {
                    q = omega;
                }
            }
            if ( elastic_gq == 1 || elastic_gg == 1 ) // in both cases ( gq->gq and gg->gg ) change the momentum 
            {
                newvecpRest=elastic->getNewMomentum(vecpRest, omega*T, q*T, random); // takes everything in GeV

                if (beta > 1e-15) // was beta == 0 but must avoid equality operations with doubles -RY
                {
                    newpx = vx*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecpRest.px() 
                      + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.py()
                      + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.pz();
                    newpy = vy*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecpRest.py()
                      + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.px()
                      + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.pz();
                    newpz = vz*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecpRest.pz()
                      + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.px()
                      + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.py();
                    

                    newvecp.px(newpx);
                    newvecp.py(newpy);
                    newvecp.pz(newpz);
                    newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
                }
                else
                    newvecp = newvecpRest;
               
                double newpRest = newvecpRest.pAbs();
                if(newpRest < pCut)
                {
                    if(outputSource)
                    {
                        // output energy deposition from jet particle absorption
                        dE = newvecpRest.pAbs();
                        dpx = newvecpRest.px();
                        dpy = newvecpRest.py();
                        dpz = newvecpRest.pz();

                        newSource.type(1);
                        newSource.tau(tau);
                        newSource.x(x);
                        newSource.y(y);
                        newSource.eta(eta);
                        newSource.dE(dE);
                        newSource.dpx(dpx);
                        newSource.dpy(dpy);
                        newSource.dpz(dpz);
                        newSource.vx(vx);
                        newSource.vy(vy);
                        newSource.vz(vz);

                        slist->push_back(newSource);
                    }
                    plist->at(i).frozen(2);
                    //newvecp *= pdummy/newvecp.pAbs();
                    // exclude soft recoil by assigning id == 0
                    if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                }

                plist->at(i).p(newvecp);           // change quark's momentum to the new one

                vecq = vecpRest - newvecpRest;             // momentum transfer q
                if ( vecq.pAbs() < pCut)//Rouz: elastic scattering counts for the gluon
                    plist->at(i).increment_soft_elastic();
                if ( vecq.pAbs() > pCut)
                    plist->at(i).increment_hard_elastic();

                if(getRecoil && omega > 0. && plist->at(i).recoil() != -1)
                {
                   if ( elastic_gq == 1 )  // quark recoil
                   {
                       // thermal quark to be recoiled
                       vecThermal = elastic->getThermalMomentum(vecq, T, 1, random);
                       vecRecoilRest = vecq + vecThermal;

                       pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                         + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                         + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                       pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                         + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                         + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                       pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                         + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                         + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                       
                       // momentum four-vector in lab frame
                       vecRecoil.px(pxRecoil);
                       vecRecoil.py(pyRecoil);
                       vecRecoil.pz(pzRecoil);
                       vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                        pyRecoil*pyRecoil +
                                        pzRecoil*pzRecoil));

                       if( vecRecoilRest.pAbs() > recoilCut )
                       {
                           col = plist->at(i).col();            // color of original gluon
                           acol = plist->at(i).acol();          // anti-color of original gluol

                           double r = random->genrand64_real1();
                           // note that when using general N_f this has to be changed
                           if (r<0.33) newID=1;
                           else if (r<0.66) newID=2;
                           else newID=3;
                           double mass;
                           if (newID<3) mass=0.33;    // set the quark's mass
                           else mass = 0.5;           // (pythia needs that for fragmentation) 

                            newOne.p(vecRecoil);

                            newOne.mass(mass);
                            double matter = random->genrand64_real1();

                            if(matter < 0.5)   // gq -> gq
                            {
                                newOne.id(newID);
                                newOne.col(1300000000+counter);
                                newOne.acol(0);

                                newHole.id(newID);
                                newHole.col(1400000000+counter);
                                newHole.acol(0);
                            }
		                    else  // g qbar -> g qbar
		                    {
                                newOne.id(-newID);
                                newOne.col(0);
                                newOne.acol(1300000000+counter);

                                newHole.id(-newID);
                                newHole.col(0);
                                newHole.acol(1400000000+counter);
		                    }

                            newOne.x(x);
                            newOne.y(y);
                            newOne.z(z);
                            newOne.frozen(0);
                            newOne.recoil(1);
                            newOne.init_counts(); // Rouz: initialize counts for the recoil particle off gluon 
                            //Sangyong's addition for formation time of radiation
                            newOne.daughter_of(-1); // daughter of nobody
                            newOne.mother_of(-1);  // mother of nobody
                            newOne.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newOne);  // add the gluon to the list of partons if (f.y>AMYpCut)

                            pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                            pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                            pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();


                            // momentum four-vector in lab frame
                            vecThermalLab.px(pxThermalLab);
                            vecThermalLab.py(pyThermalLab);
                            vecThermalLab.pz(pzThermalLab);
                            vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                  pyThermalLab*pyThermalLab +
                                                  pzThermalLab*pzThermalLab));

                            newHole.p(vecThermalLab);
                            newHole.mass(mass);

                            newHole.x(x);
                            newHole.y(y);
                            newHole.z(z);
                            newHole.frozen(0);
                            newHole.recoil(-1);
                            
                            //Sangyong's addition for formation time of radiation
                            newHole.daughter_of(-1); // daughter of nobody
                            newHole.mother_of(-1);  // mother of nobody
                            newHole.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newHole);    // add hole to the parton list
                       }
                       else
                       {
                           if(outputSource)
                           {
                               // momentum transfer q goes into medium
                               dE = vecq.e();
                               dpx = vecq.px();
                               dpy = vecq.py();
                               dpz = vecq.pz();

                               newSource.type(1);
                               newSource.tau(tau);
                               newSource.x(x);
                               newSource.y(y);
                               newSource.eta(eta);
                               newSource.dE(dE);
                               newSource.dpx(dpx);
                               newSource.dpy(dpy);
                               newSource.dpz(dpz);
                               newSource.vx(vx);
                               newSource.vy(vy);
                               newSource.vz(vz);

                               slist->push_back(newSource);
                           }
                       }
                   }
                   else if ( elastic_gg == 1 )  // gluon recoil
                   {
                        // thermal gluon to be recoiled
                        vecThermal = elastic->getThermalMomentum(vecq, T, -1, random);
                        vecRecoilRest = vecq + vecThermal;

                        pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                          + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                          + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                        pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                          + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                          + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                        pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                          + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                          + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                        
                        // momentum four-vector in lab frame
                        vecRecoil.px(pxRecoil);
                        vecRecoil.py(pyRecoil);
                        vecRecoil.pz(pzRecoil);
                        vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                         pyRecoil*pyRecoil +
                                         pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                            newOne.p(vecRecoil);
                            newOne.id(21);
                            newOne.mass(0.);
                            newOne.col(1500000000+counter);
                            newOne.acol(1600000000+counter);
                            newOne.x(x);
                            newOne.y(y);
                            newOne.z(z);
                            newOne.frozen(0);
                            newOne.recoil(1);
                            newOne.init_counts(); // Rouz: recoil initialize counts
                            //Sangyong's addition for formation time of radiation
                            newOne.daughter_of(-1); // daughter of nobody
                            newOne.mother_of(-1);  // mother of nobody
                            newOne.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newOne);    // add recoil parton to the parton list

                            pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                            pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                            pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();


                            // momentum four-vector in lab frame
                            vecThermalLab.px(pxThermalLab);
                            vecThermalLab.py(pyThermalLab);
                            vecThermalLab.pz(pzThermalLab);
                            vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                  pyThermalLab*pyThermalLab +
                                                  pzThermalLab*pzThermalLab));

                            newHole.p(vecThermalLab);

                            newHole.id(21);
                            newHole.mass(0.);
                            newOne.col(1700000000+counter);
                            newOne.acol(1800000000+counter);

                            newHole.x(x);
                            newHole.y(y);
                            newHole.z(z);
                            newHole.frozen(0);
                            newHole.recoil(-1);

                            //Sangyong's addition for formation time of radiation
                            newHole.daughter_of(-1); // daughter of nobody
                            newHole.mother_of(-1);  // mother of nobody
                            newHole.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newHole);    // add hole to the parton list
                        }
                        else
                        {
                            if(outputSource)
                            {
                                // momentum transfer q goes into medium
                                dE = vecq.e();
                                dpx = vecq.px();
                                dpy = vecq.py();
                                dpz = vecq.pz();

                                newSource.type(1);
                                newSource.tau(tau);
                                newSource.x(x);
                                newSource.y(y);
                                newSource.eta(eta);
                                newSource.dE(dE);
                                newSource.dpx(dpx);
                                newSource.dpy(dpy);
                                newSource.dpz(dpz);
                                newSource.vx(vx);
                                newSource.vy(vy);
                                newSource.vz(vz);

                                slist->push_back(newSource);
                            }
                        }
                    }
                }
                else
                {
                    if(outputSource)
                    {
                        // momentum transfer q goes into medium
                        dE = vecq.e();
                        dpx = vecq.px();
                        dpy = vecq.py();
                        dpz = vecq.pz();

                        newSource.type(1);
                        newSource.tau(tau);
                        newSource.x(x);
                        newSource.y(y);
                        newSource.eta(eta);
                        newSource.dE(dE);
                        newSource.dpx(dpx);
                        newSource.dpy(dpy);
                        newSource.dpz(dpz);
                        newSource.vx(vx);
                        newSource.vy(vy);
                        newSource.vz(vz);

                        slist->push_back(newSource);
                    }
                }
            }
        }
    } // loop over partons i

    return counter;
}

void MARTINI::sampleTA()
{
    int A = static_cast<int>(glauber->get_Projetile_A());
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
   if(fullEvent==1)
   {
      // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
      // then lay them on top of each other and sample the number of jet events in a tube with area=inelasticXSec
      
      //cout << " number of binary collisions = " << glauber->TAB() << endl;
      //cout << " A=" << glauber->nucleusA() << endl;
      Parton parton;                               // new Parton object
      double b = glauberImpactParam;
      int posx;
      int posy;
      int n1=0;
      int n2=0;
      double r = 1.2*pow(glauber->get_Projetile_A(), 1./3.);

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
      int Projectile_A=glauber->get_Projetile_A();
      int Target_A=glauber->get_Projetile_A();
      int Projectile_Z = glauber->get_Projetile_Z();
      int Target_Z = glauber->get_Projetile_Z();
      //cout << " probability to generate jet event =" << jetXSec/inelasticXSec << endl;
      for (int ix = 0; ix < ixmax; ix++) 
      {
        for (int iy = 0; iy < ixmax; iy++)
        {
          for (int i = 0; i < nucALat[ix][iy]; i++)
          {
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
                if(random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A) 
                  n1=2112;
                else 
                  n1=2212;
                if(random->genrand64_real1() > (double)Target_Z/(double)Target_A) 
                  n2=2112;
                else 
                  n2=2212;
                pythia.next(n1,n2);                                                // generate event with pythia
                done = 0;
                //pythia.event.list();
                for (int ip = 0; ip < pythia.event.size(); ++ip) 
                {
                  // if the parton is final, i.e., present after the showers, then put it in the list
                  // currently heavy quarks are taken but will not lose energy! 
                  // only g and u,d,s lose energy while the 3 quarks are assumed to be massless
                  if (pythia.event[ip].isFinal()
                      // && (pythia.event[i].id()<4 || pythia.event[i].id()==21)
                     )
                  {
                    parton.id(pythia.event[ip].id());                     // set parton id
                    parton.status(pythia.event[ip].status());             // set parton status
                    parton.mass(pythia.event[ip].m());                    // set mass
                    if (fixedTemperature==0)
                    {
                      //parton.x(xmin+(ix)*cellLength+cellLength/2.);       // set position
                      //parton.y(xmin+(iy)*cellLength+cellLength/2.);
                      parton.x(xmin+(ix)*cellLength+xPositionInCell);       // set position
                      parton.y(xmin+(iy)*cellLength+yPositionInCell);
                      //parton.xini(xmin+(ix)*cellLength+xPositionInCell);       // set position
                      //parton.yini(xmin+(iy)*cellLength+yPositionInCell);
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
                    }
                    parton.z(0.);
                    parton.col(pythia.event[ip].col());                     // set color 
                    parton.acol(pythia.event[ip].acol());                 // set anti-color
                    parton.frozen(0);                                     // parton is not frozen (will evolve)
                    parton.p(pythia.event[ip].p());                       // set momentum
                    parton.source(0);                                     // All initial partons have itsSource=0
                  
                    // Sangyong's addition
                    parton.mother_of(-1); // mother of nobody
                    parton.daughter_of(-1); // daughter of nobody
                    parton.p_at_split(parton.p()); // original momentum, don't need it but might as well populate it
                    // Sangyong's addition

                    plist->push_back(parton);                             // add the parton to the main list
                  }
                } 
                eventNumber++;
              }
            }
          }
        }
      }

      //for (int i=0; i<100; i++)
      //{
      //  cout << "collisions for parton i=" << i << " =" << countCollisionsA[i] << endl;
      //  cout << "collisions for parton j=" << i << " =" << countCollisionsB[i] << endl;
      //}

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
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
          if ( j == ixmax-1 ) fout1 << endl;
        }
      }
      
      for (int i = 0; i < ixmax; i++)
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
          if ( j == ixmax-1 ) fout2 << endl;
        }
      }
      
      for (int i = 0; i < ixmax; i++)
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
          if ( j == ixmax-1 ) fout3 << endl;
        }
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
        Parton parton;   // new Parton object
        int n1;
        int n2;
        int Projectile_A = glauber->get_Projetile_A();
        int Target_A = glauber->get_Target_A();
        int Projectile_Z = glauber->get_Projetile_Z();
        int Target_Z = glauber->get_Target_Z();
        if(random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A)
            n1=2112;
        else 
            n1=2212;
        if(random->genrand64_real1() > (double)Target_Z/(double)Target_A)
            n2=2112;
        else 
            n2=2212;
        pythia.next(n1,n2);    // generate event with pythia
        //pythia.event.list(); 
        
        //ReturnValue rv;
        //if (evolution && fixedTemperature==0) 
        if (!allFromCenter)
        {
            if(nbinFromFile == 1)
            {
                //Sample a random element of the list of number of collisions. -CFY 11/2/2010
                double randomr = random->genrand64_real2();
                int randomint = (int)(randomr*((double)Nbin));
                rv.x = binary[randomint][0];
                rv.y = binary[randomint][1];
                cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
            }
            else if(nbinFromFile == 2 or nbinFromFile == 3 or nbinFromFile == 4)
            {
               double x_sample, y_sample;
               binary_info_ptr->get_sample(&x_sample, &y_sample);
               rv.x = x_sample;
               rv.y = y_sample;
            }
            else if (glauberEnvelope)
                rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
            else 
                rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
            //cout << rv.x << " " << rv.y << endl;
        }
        else
        {
            rv.x=0.;
            rv.y=0.;
        }

        for (int i = 0; i < pythia.event.size(); ++i) 
        {
            // if the parton is final, i.e., present after the showers, then put it in the list
            // currently heavy quarks are taken but will not lose energy! 
            // only g and u,d,s lose energy and the 3 quarks are assumed to be massless
            if (pythia.event[i].isFinal()
                // && (pythia.event[i].id()<4 || pythia.event[i].id()==21)
               )
            {
                parton.id(pythia.event[i].id());         // set parton id
                parton.status(pythia.event[i].status()); // set parton status
                parton.mass(pythia.event[i].m());        // set mass
                
                parton.x(rv.x);                          // set position
                parton.y(rv.y);
                parton.z(0.);
                parton.col(pythia.event[i].col());       // set color 
                parton.acol(pythia.event[i].acol());     // set anti-color
                parton.frozen(0);                        // parton is not frozen (will evolve)
                
                parton.p(pythia.event[i].p());           // set momentum
                parton.source(0);                        // All intial partons have itsSource=0
                parton.recoil(0);                        // Hard partons have itsRecoil=0
                        
                // Sangyong's addition
                parton.mother_of(-1);                   // mother of nobody
                parton.daughter_of(-1);                 // daughter of nobody
                parton.p_at_split(parton.p());          // original momentum, don't need it but might as well populate it
                // Sangyong's addition
                
                plist->push_back(parton);               // add the parton to the main list
                // check jettiness of di-quarks:
                // if (abs(parton.id())>2000) 
                //   cout << "id=" << parton.id() << " px=" << parton.p().px() 
                //        << " py=" << parton.p().py() << " pz=" <<  parton.p().pz() << endl;
            }
        } 
      return 1;
    }
}


int MARTINI::fragmentation(vector<Parton> * plist, int recoil, int currentEvent )
{
    if (fragmentationMethod == 1)
    {
        if (fullEvent == 0 )
        {
            // clear the event record to fill in the partons resulting from evolution below
            pythia.event.reset(); 
            int id, col, acol, status;
            int numEvents = 0;
            double mass;
            Vec4 pvec; 
            vector<Parton> * flist = new vector<Parton>;
            
            for(int p=0; p<plist->size(); p++)
                if( plist->at(p).recoil() == recoil and plist->at(p).id() != 0 )
                    flist->push_back(plist->at(p));
            // reconnect color for recoils/holes
            if(recoil)
            {

                //vector<Parton> *templist;
                //templist = new vector<Parton>;
                //templist->insert(templist->begin(), flist->begin(), flist->begin()+10);
                ////templist->insert(templist->begin(), flist->begin(), flist->end());
                //flist->clear();
                //flist = templist;
                
                //Vec4 pvec;
                //int col, acol;
                //int f_imax = flist->size();
                //for ( int f_i=0; f_i<f_imax; f_i++)
                //{
                //    pvec = flist->at(f_i).p();
                //    col = flist->at(f_i).col();
                //    acol = flist->at(f_i).acol();
                //    cout << "hadronization::col/acol = " << col << " " << acol << " p = " << pvec;
                //}
                
                
                //flist->at(0).id(1);
                //flist->at(1).id(21);
                //flist->at(2).id(21);
                //flist->at(3).id(1);
                //flist->at(4).id(21);
                //flist->at(5).id(21);


                cout << setprecision(3);
                // loop over all partons in the sublist and reset color links
                for (int k=0; k<flist->size(); k++)  
                {
                    flist->at(k).col(0);
                    flist->at(k).acol(0);
                }
                
                unsigned it = 0;
                while(it < flist->size())
                {
                    double T = hydroTfinal;
                    Parton newOne;
                    
                    double r, rn;
                    int newID;
                    double mass;

                    bool quarkLoop = true;

                    // search upcoming q/qbar
                    if (flist->at(it).id() > 0 and flist->at(it).id() < 4)
                    {
                        quarkLoop = true;
                    }
                    else if (flist->at(it).id() < 0)
                    {
                        quarkLoop = false;
                    }
                    // if next parton is g, search upcoming q/qbar
                    else if (flist->at(it).id() == 21)
                    {  
                        //cout << "loop begins with gluon" << endl;
                        int pos = -1;
                        unsigned ii = it+1;
                        while(ii < flist->size())
                        {
                            if(flist->at(ii).id() != 21)
                            {
                                pos = ii;
                                if(flist->at(ii).id() > 0 and flist->at(ii).id() < 4)
                                    quarkLoop = true;
                                else
                                    quarkLoop = false;

                                break;
                            }
                            ii++;
                            if(ii == flist->size()) 
                            {
                                //cout << "no more quark to start loop! Add one at the front." << endl;
                                //cout << "[test]::it = " << it << endl;

                                r = random->genrand64_real1();
                                // note that when using general N_f this has to be changed
                                if (r<0.33) newID=1;
                                else if (r<0.66) newID=2;
                                else newID=3;
                                rn=random->genrand64_real1();
                                if(rn < 0.5) newID *= -1;
                                if (fabs(newID<3)) mass=0.33;
                                else mass = 0.5;

                                newOne.id(newID);
                                newOne.mass(mass);
                                newOne.col(0);
                                newOne.acol(0);

                                newOne.p(random->thermal(T, 1));

                                //cout << "put " << newOne.id() << " at " << it+1 << endl;
                                flist->insert(flist->begin()+it, newOne);
                                //for(unsigned i = 0; i != flist->size(); i++) {
                                //   cout << flist->at(i).id() << " ";
                                //}
                                //cout << endl;
                                break;
                            }
                        }
                        //cout << "pos = " << pos << endl;
                        if(pos > 0)
                        {
                            // rotate current gluon and q/qbar in pos
                            //cout << "[1]rotate::from " << pos+1 << " to " << it+1 << endl;
                            rotate( flist->begin()+pos, flist->begin()+pos+1, flist->begin()+it+1 );
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;  
                        }
                    }
                    //cout << "[Interm]::Loop is always started with q/qbar " << flist->at(it).id() << " at " << it+1 << endl;
                    //for(unsigned i = 0; i != flist->size(); i++) {
                    //cout << flist->at(i).id() << " ";
                    //}
                    //cout << endl;
                    // Now go to the position where there is non gluon
                    while(it < flist->size())
                    {
                        it++;
                        // If this q/qbar is the last, add qbar/q to close the loop
                        if(it == flist->size()) 
                        {
                            //cout << "gluon ending or this q/qbar is end! Add one" << endl;
                            r = random->genrand64_real1();
                            // note that when using general N_f this has to be changed
                            if (r<0.33) newID=1;
                            else if (r<0.66) newID=2;
                            else newID=3;
                            if(quarkLoop) newID *= -1;
                            if (fabs(newID<3)) mass=0.33;
                            else mass = 0.5;

                            newOne.id(newID);
                            newOne.mass(mass);
                            newOne.col(0);
                            newOne.acol(0);

                            newOne.p(random->thermal(T, 1));

                            //cout << "put " << newOne.id() << " at " << it+1 << endl;
                            flist->insert(flist->begin()+it, newOne);
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //   cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                        // Position where there is non-gluon
                        if(flist->at(it).id() != 21)
                            break;

                    }
                    //cout << "[Interm]::next q/bar position : " << it+1 << endl;

                    // 1: close the loop if next one is its counterpart (do nothing)
                    // 2: replace next q/qbar with upcoming counter part with next one
                    // 3: create its counter part if none is found
                    unsigned ii = it;  // ii : position of q/qbar for closing the loop
                    while(ii < flist->size())
                    {
                        // If we find q/qbar to close the loop.
                        if( (quarkLoop && flist->at(ii).id() < 0) ||
                            (!quarkLoop && flist->at(ii).id() > 0 && flist->at(ii).id() < 4) )
                        {
                            //cout << "[2]rotate::from " << ii+1 << " to " << it+1 << endl;
                            rotate( flist->begin()+ii, flist->begin()+ii+1, flist->begin()+it+1 );
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //    cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                        ii++;
                        // If we don't find any q/qbar to close the loop.
                        if(ii == flist->size()) 
                        {
                            //cout << "no more q/qbar to close loop! Add one" << endl;
                            // add q/qbar if nothing is found
                            unsigned iq = it;  // iq : position of q/qbar to add
                            while(iq < flist->size())
                            {
                                if(flist->at(iq).id() != 21 || iq == flist->size())
                                    break;
                                iq++;
                            }

                            r = random->genrand64_real1();
                            // note that when using general N_f this has to be changed
                            if (r<0.33) newID=1;
                            else if (r<0.66) newID=2;
                            else newID=3;
                            if(quarkLoop) newID *= -1;
                            if (fabs(newID<3)) mass=0.33;
                            else mass = 0.5;

                            newOne.id(newID);
                            newOne.mass(mass);
                            newOne.col(0);
                            newOne.acol(0);

                            newOne.p(random->thermal(T, 1));

                            //cout << "put " << newOne.id() << " at " << it+1 << endl;
                            flist->insert(flist->begin()+iq, newOne);

                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //    cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                    }

                    it++;
                    if(it == flist->size()) 
                    {
                        //cout << "process finished!" << endl;
                        break;
                    }
                    
                    //cout << "[finish]::it = " << it+1 << " new flist->size() = " << flist->size() << endl;
                    //if(it > 5000) exit(3);
                }
                // Now assign col/acol accordingly.
                int stringCounter = 100;
                unsigned ip = 0;
                while(ip < flist->size())
                {
                
                    bool quarkLoop = true;
                    if(flist->at(ip).id() < 0) quarkLoop = false;
                    //cout << "quarkLoop = " << quarkLoop << endl;
                
                    if(quarkLoop)
                        flist->at(ip).acol(0);
                    else
                        flist->at(ip).col(0);
                
                    // closeLoop = true : ending a loop
                    bool closeLoop = false;
                    unsigned ii = ip+1;
                    while(ii < flist->size())
                    {
                        double px1, py1, pz1, p1;
                        double px2, py2, pz2, p2;
                        double cosP1P2;
                
                        px1 = flist->at(ii-1).p().px();
                        py1 = flist->at(ii-1).p().py();
                        pz1 = flist->at(ii-1).p().pz();
                        p1 = flist->at(ii-1).p().pAbs();
                
                        px2 = flist->at(ii).p().px();
                        py2 = flist->at(ii).p().py();
                        pz2 = flist->at(ii).p().pz();
                        p2 = flist->at(ii).p().pAbs();
                        cosP1P2 = (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2);
                
                        if(quarkLoop)
                        {
                            flist->at(ii-1).col(stringCounter);
                            flist->at(ii).acol(stringCounter);
                            if(flist->at(ii).id() != 21) closeLoop = true;
                        }
                        else
                        {
                            flist->at(ii).col(stringCounter);
                            flist->at(ii-1).acol(stringCounter);
                            if(flist->at(ii).id() != 21) closeLoop = true;
                        }
                
                        stringCounter++;
                        if(closeLoop)
                        {
                            ip = ii;
                            break;
                        }
                        ii++;
                    }
                    if(quarkLoop)
                        flist->at(ip).col(0);
                    else
                        flist->at(ip).acol(0);
                
                    ip++;
                }


            }
            // end of color reconnection

            int f_imax = flist->size();

            // add all the higher momentum partons:
            // loop over all partons in the main list
            for ( int f_i=0; f_i<f_imax; f_i++)
            {
                id = flist->at(f_i).id();

                status = flist->at(f_i).status();
                if ( status<1 or status>100)   // if status was not set, cure this here
                      status = 62;
                

                pvec = flist->at(f_i).p();
                col = flist->at(f_i).col();
                acol = flist->at(f_i).acol();
                mass = flist->at(f_i).mass();

                pythia.event.append(id,status,col,acol,pvec,mass);
            }
            // does the fragmentation and the rest (if success add 1 to total events)
            if(pythia.forceHadronLevel()) 
                numEvents+=1; 

            delete flist;
            return numEvents;
        }
        else if (fullEvent == 1)
        {
            //pythia.event.list(); 
            pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
            int j = 0;
            int id, col, acol, sid, scol, sacol, status, eventNumber;
            int p_imax = plist->size();
            int p_jmax = plist->size();
            int numEvents = 0;
            //int check = 0;
            double mass, smass;
            Vec4 pvec, spvec; 
            // add all the higher momentum partons:
            for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons in the main list (those above momentum threshold)
            {
                id = plist->at(p_i).id();
                if (id == 22)   // do not feed photons into fragmentation 
                    continue;   // in this way only fragmentation and decay photons appare at the end
                                // photons from initial state and jet-medium interaction will be collected 
                                // before fragmenetation in the main.cpp

                status = plist->at(p_i).status();
                //if (status<1 || status>100) status = 1; // if status was not set, cure this here
                if (status<1 || status>100)   // if status was not set, cure this here
                      status = 62;

                pvec = plist->at(p_i).p();
                col = plist->at(p_i).col();
                acol = plist->at(p_i).acol();
                mass = plist->at(p_i).mass();
                // append( id, status, col, acol, p, mass )
                if (eventNumber==currentEvent)
                {
                    //if(id == 22)
                    //{
                    //  int source = plist->at(p_i).source();
                    //  cout << scientific << setw(16) << setprecision(6) 
                    //       << "photon: status = " << status << " source = " << source
                    //       << ", px = " << pvec.px() << ", py = " << pvec.py() << ", pz = " << pvec.pz() 
                    //       << endl;
                    //}
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
        int p_imax = plist->size();
        int p_jmax = plist->size();
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
        
        double eps = 1e-15;
        hydroZmax = 40;
        //cout << "doing fragmentation, method 2" << endl;
        if (hydroWhichHydro != 3 && hydroWhichHydro != 4 && hydroWhichHydro != 7 && hydroWhichHydro !=8) 
        {    
            hydroZmax=hydroTauMax;
            hydroDz=hydroDx;
        }
      
        DX = 1.; //fm
        DZ = 1.; //fm
        // these values are for internal consumption of this function
        
        ixmax = floor((2.*hydroXmax)/DX+0.0001);
        izmax = floor((2.*hydroZmax)/DZ+0.0001);//z is really eta

        Ncells = ixmax*ixmax*izmax;

        vector<Parton> ** psublist;               // pointer to array of vector<Parton> objects
        vector<double> ** cvec;

        psublist = new vector<Parton> *[Ncells];  // pointer to array of vector<Parton> objects
        cvec = new vector<double> *[Ncells]; 

        //cout << "*****************************"<<endl;
        //cout << "\t In MARTINI::fragmentation mode 2: "<<endl;
        //cout << "\t\tixmax: "<<ixmax << ", izmax: "<<izmax<<endl;
        //cout << "\t\thydroZmax: "<< hydroZmax << ", hydroXmax: "<<hydroXmax <<endl;
        //cout << "\t\tNCells : "<<Ncells<< endl;
        //cout << "*****************************"<<endl;
        
        for(i=0; i < Ncells; i++)
        {
            psublist[i] = new vector<Parton>;// psublist[i] is the ith sublist that holds 
                                             // the high momentum partons that were evolved
            cvec[i] = new vector<double>;
        }
      
        countQuarkJets = 0;
        countAntiQuarkJets = 0;
      
        for ( int p_i=0; p_i<p_imax; p_i++ ) // loop over all partons in the main list 
                                             // (those above momentum threshold) 
        {                                    // to sort them into sub-lists
            id = plist->at(p_i).id();
            if (abs(id) > 22)//throw out all except q, qbar and g
                continue; 

            plist->at(p_i).status(1);
            x = plist->at(p_i).x();          // x value of position in [fm]
            y = plist->at(p_i).y();          // y value of position in [fm]
            z = plist->at(p_i).z();          // z value of position in [fm]
            t = plist->at(p_i).tFinal();     // t value of lab time in [fm/c]

            if (fabs(t) < eps)
                continue;
    
            ix = floor((hydroXmax+x)/DX+0.0001);  // x-coordinate of the cell we are in now
            iy = floor((hydroXmax+y)/DX+0.0001);  // y-coordinate of the cell we are in now
                                                  // note that x and y run from -hydroXmax to +hydroXmax
                                                  // and ix and iy from 0 to 2*hydroXmax
                                                  // hence the (hydroXmax+x or y) for both
            iz = floor((hydroZmax+z)/DZ+0.0001);// z-coordinate of the cell we are in now
  
    
            i=ix+ixmax*(iy+iz*(ixmax)); //determine sublist number by position
            if (i > Ncells-1)
            {
                cout << "*****************************"<<endl;
                cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;
                cout << "ix=" << ix << ", iy=" << iy << ", iz=" << iz << ", i=" << i << " out of " << Ncells << endl;
                cout << "\t In MARTINI::fragmentation mode 2: "<<endl;
                cout << "\t\tixmax: "<<ixmax << ", izmax: "<<izmax<<endl;
                cout << "\t\thydroZmax: "<< hydroZmax << ", hydroXmax: "<<hydroXmax <<endl;
                cout << "\t\tNCells : "<<Ncells<< endl;
                cout << "*****************************"<<endl;
                exit(-1);
            }
            if (plist->at(p_i).p().pAbs()<2.) // dont add parton under 2 GeV
                continue; 
            
            if (abs(id)<21)
            {
                if (id > 0) 
                    countQuarkJets++;
                else
                    countAntiQuarkJets++;
            }

            
            psublist[i]->push_back(plist->at(p_i)); // add the parton to the correct sublist.
            
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

        for(i=0; i < Ncells; i++) //loop over all sublists (cells)
        {
            if (psublist[i]->size() == 0)
                continue;
            //cout << "There is/are " << psublist[i]->size() << " hard parton(s) in list " << i << endl;
            x = cvec[i]->at(0);
            y = cvec[i]->at(1);
            z = cvec[i]->at(2);
            t = cvec[i]->at(3);
            
            //cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;
            
            //get the temperature and flow velocity at the current position and time:
            if(fixedTemperature==0)
            {
                hydroInfo = hydroSetup->getHydroValues(x, y, z, t);
                T = hydroInfo.T;
            }
            else
                T = fixedT;
            
            if (T==0.) 
                T=0.1; // give it a minimum T
            //cout << " T in this cell =" << T << " GeV" << endl;
            
            V = DX*DX*DZ;
            
            Ngluons = (V*16./(PI*PI)*T*T*T*1.202056903/pow(hbarc,3.)); // 1.202056903 is Riemann zeta(3).
            Nquarks = (V*9./(2.*PI*PI)*T*T*T*1.202056903/pow(hbarc,3.));
            
            //cout << "N_gluons in this cell = " << Ngluons << endl;
            //cout << "N_quarks and anti-quarks in this cell = " << 6*Nquarks << endl;
            
            Igluons = round(Ngluons);
            Iquarks = round(6*Nquarks); // three flavors
            
            for (k=0; k < Ngluons; k++)
            {
                newOne.id(21);                   // add a thermal gluon
                newOne.mass(0.); 
                newOne.frozen(2);
                newOne.x(x);                     // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, -1));// sample momentum from Bose distribution at temperature T
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
                newOne.frozen(2);
                newOne.x(x);                    // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, 1));// sample momentum from Fermi distribution at temperature T
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
         
            //is = 0; 
            //while (psublist[i]->begin()+is<psublist[i]->end())  // remove heavy quarks (keep photons)
            //{
            //  if (psublist[i]->at(is).id()>22) 
            //  {
            //      psublist[i]->erase(psublist[i]->begin()+is);
            //      is--;
            //  }
            //  is++;
            //}

            for (k=0; k < psublist[i]->size(); k++)  // loop over all partons in the sublist and reset color links
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
                if (psublist[i]->at(k).id()==21 and 
                    ((psublist[i]->at(k).col()!=0 and psublist[i]->at(k).acol()==0) or 
                    (psublist[i]->at(k).acol()!=0 and psublist[i]->at(k).col()==0)) )
                {
                    if (psublist[i]->at(k).status()==1)
                      //cout << " gluon jet needs extra string attached!" << endl;
                    rn=random->genrand64_real1();         // add a thermal quark or anti-quark
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
                    newOne.frozen(2);
                    newOne.x(x);                              // set the new parton's initial position
                    newOne.y(y);
                    newOne.z(z);
                    newOne.p(random->thermal(T, 1)); //T         // sample momentum from Fermi distribution at temperature T
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
                if (psublist[i]->at(is).id()!=22 and psublist[i]->at(is).col()==0 and psublist[i]->at(is).acol()==0)
                {
                    //cout << "erasing position " << is << " with id=" << psublist[i]->at(is).id() << " and col/acol=" 
                    //   << psublist[i]->at(is).col() << "/" << psublist[i]->at(is).acol() << endl;
                    psublist[i]->erase(psublist[i]->begin()+is);
                    is--;
                }
                is++;
            }

            //for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist and print content
            //{
            //  if (psublist[i]->at(k).status() == 1)
            //      cout << " cell " << i << ": jet " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //  else if (psublist[i]->at(k).status() == 0)
            //      cout << " cell " << i << ": thermal parton " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //  else
            //      cout << " cell " << i << ": weird parton " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //}

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
    pythia.setJetpTmin(jetpTmin);
    pythia.setJetpTmax(jetpTmax);
      
    // initialize PYTHIA in the CM frame.
    pythia.init( 2212, 2212, cmEnergy ); // 2 protons, (if 14000 each 7000 GeV) /2212 is proton
    cout << endl << "PYTHIA initialized by MARTINI with sqrt(s)=" << cmEnergy << " GeV." << endl;
    
    // set the relevant cross sections that PYTHIA was so kind to compute
    elasticXSec = pythia.elasticCrossSection();
    totalXSec = pythia.totalCrossSection();
    jetXSec = pythia.jetCrossSection();
    inelasticXSec = totalXSec-elasticXSec;
    
    // output the cross sections
    cout << endl << "[MARTINI::initPythia]:" << endl;
    cout << "Jet cross section = " << scientific << setprecision(6) << setw(12) << jetXSec << " mb, with a pTmin = " << jetpTmin << " GeV." << endl;
    cout << "Total cross section = " << totalXSec << " mb." << endl;
    cout << "Total inelastic cross section = " << inelasticXSec << " mb." << endl;
  
}

bool MARTINI::init(int path_number)
{
    /// set all parameters as they were initialized
    int seed = 0;
    if (path_number > 0)
    {
        seed += path_number;
    }
    // For FIC:
    file_number = seed%5 +1;
    cout << "file_number = " << file_number << endl;
    seed *= 10000;
    pCut = settings.parm("General:pCut");
    eLossCut = settings.parm("General:eLossCut");
    fixedT = settings.parm("General:Temperature");
    maxTime = settings.parm("General:MaxTime");
    dtfm = settings.parm("General:TimeStep");
    alpha_s = settings.parm("General:AlphaS");
    runs = settings.mode("General:Events");
    Nf = settings.mode("General:NumberOfFlavors");
    jetpTmin = settings.parm("General:JetPTMin");
    jetpTmax = settings.parm("General:JetPTMax");
    fullEvent = settings.mode("General:FullEvent");
    moveBeforeTau0 = settings.mode("General:MoveBeforeTau0");
    cmEnergy = settings.parm("General:cmEnergy");
    fullVacuumShower = settings.flag("General:FullVacuumShower");
    Ldependence = settings.mode("General:Ldependence");
    fragmentationMethod = settings.mode("General:FragmentationMethod");
    rateSelector = settings.mode("General:RadiativeRateSet");
    examineHQ = settings.mode("General:examineHQ");
    cout << "examineHQ = " << examineHQ << endl;
    setInitialPToZero = settings.mode("General:setInitialPToZero");
    cout << "setInitialPToZero = " << setInitialPToZero << endl;
    charmWidth = settings.parm("General:charmWidth");
    bottomWidth = settings.parm("General:bottomWidth");
    T_C_HQ = settings.parm("General:T_C_HQ");
    TwoPiTD_HQ = settings.parm("General:TwoPiTD_HQ");
    totalHQXSec = settings.parm("General:totalHQXSec");
    tauEtaCoordinates = settings.mode("General:tauEtaCoordinates");
    cout << "tauEtaCoordinates = " << tauEtaCoordinates << endl;

    /***** Jet Conversion Additions: Rouz *****/
    useRateTableforConversion = settings.flag("General:UseConvPhotonsTable");
    cout << "General:UseConvPhotonsTable = "<<useRateTableforConversion<<endl;
    if (useRateTableforConversion && (alpha_s > 0.32 || alpha_s < 0.28) )
    {
      cout<<"Jet Conversion using interpolation of the rate table is only available for alpha_s=0.3. Exiting." << endl;
      exit(-1);
    }
    if (useRateTableforConversion && alpha_s > 0.29 && alpha_s < 0.31)
    {
      /* I put alpha_s in (0.29, 0.31) range because we usually take alpha_s at RHIC or LHC energies to be
      * one of the following values: 0.2, 0.27, 0.3, 0.4 and I didn't want to perform an exact equality 
      * test with a double.
      * Now add the codes to read the files here. There are two files to be read
      *   1. grid information file, containing grid start and end points in 3 dimensions & grid spacing
      *   2. grid file itself.
      * Start by reading in the grid information:
      */
     cout<<"Reading conversion table data: "<<endl;
     string table_settings = settings.word("General:conversionTableSettings");
     string grid_of_rates = settings.word("General:conversionTableFile");
     cout<<" Table Settings at: "<<table_settings<<endl;
     cout<<" grid at: "<<grid_of_rates<<endl;
     angularConvGrid = new ConversionAngularGrid(table_settings, grid_of_rates);
     cout<<"Read in the jet conversion grid file."<<endl;
    }
    // *************** Done with conversion photon table stuff **********************// 
    // Check if we want a perturbative calculation:
    pertCalc = settings.flag("General:PerturbativeCalculation");
    exaggerate_brem = settings.parm("General:IncreaseRateBrem");
    exaggerate_conv = settings.parm("General:IncreaseRateConv");
    cout<<"General:PerturbativeCalculation: "<<pertCalc<<endl;
    //if ( pertCalc )
    //{
        //// perturbative calculation
        //pTmin_pert  = settings.parm("General:pTmin");
        //pTmax_pert  = settings.parm("General:pTmax"); 
        //nbins_pert  = settings.parm("General:NumBins"); 
        //etaCut_pert = settings.parm("General:etaCut");
        //conversion_photons = new TH1D("conv_gamma_hist", "pT hist of conv gamma", nbins_pert, pTmin_pert, pTmax_pert);
        //amy_photons        = new TH1D("amy_gamma_hist", "pT hist of amy gamma", nbins_pert, pTmin_pert, pTmax_pert);

        //TH1::SetDefaultSumw2(); // so that ROOT would keep stat errors
    //}
    //-----------------------------------------------------------------------------//
    if (rateSelector < 0 || rateSelector > 4)
    {
        cout << "The chosen radiative rates are not available at the moment. Exiting." << endl;
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
    getRecoil = settings.flag("General:GetRecoil"); 
    cout << "getRecoil = " << getRecoil << endl;
    recoilCut = settings.parm("General:RecoilCut");
    cout << "recoilCut = " << recoilCut << endl;
    outputSource = settings.flag("General:OutputSource"); 
    cout << "outputSource = " << outputSource << endl;
    trackHistory = settings.flag("General:TrackHistory"); 
    initialXjet = settings.parm("General:InitialXjet");
    initialYjet = settings.parm("General:InitialYjet");
    initialPXjet = settings.parm("General:InitialPXjet");
    initialPYjet = settings.parm("General:InitialPYjet");
    allFromCenter = settings.flag("General:allFromCenter");
    pCutPropToT = settings.flag("General:pCutPropToT");
    formationTime = settings.flag("General:FormationTime");
    cout << "formationTime = " << formationTime << endl;
    runningElas = settings.flag("General:RunningElas");
    cout << "runningElas = " << runningElas << endl;
    runningRad = settings.flag("General:RunningRad");
    cout << "runningRad = " << runningRad << endl;
    renormFacRad = settings.parm("General:RenormFacRad");
    cout << "renormFacRad = " << renormFacRad << endl;
    renormFacElas = settings.parm("General:RenormFacElas");
    cout << "renormFacElas = " << renormFacElas << endl;
    
    hydroTau0 = settings.parm("Hydro:tau0");
    hydroTauMax = settings.parm("Hydro:taumax");
    hydroDtau = settings.parm("Hydro:dtau");
    hydroDx = settings.parm("Hydro:dx");
    //Hydrodynamical output saved in tau-eta coordinates will simply use hydroDz and hydroZmax:
    hydroDz = settings.parm("Hydro:dz");
    hydroXmax = settings.parm("Hydro:xmax");
    hydroZmax = settings.parm("Hydro:Zmax");
    hydro_nskip_tau = settings.parm("Hydro:nskip_tau");
    hydro_nskip_x = settings.parm("Hydro:nskip_x");   // for furture
    hydro_nskip_z = settings.parm("Hydro:nskip_z");   // for furture
    hydroTfinal = settings.parm("Hydro:Tfinal");
    hydroWhichHydro = settings.mode("Hydro:WhichHydro");  
    hydroSubset = settings.mode("Hydro:Subset");
    hydroViscous = settings.flag("Hydro:Viscous");  
    fixedDistribution = settings.flag("Hydro:fixedDistribution");
    
    hydroDtau = hydroDtau*hydro_nskip_tau;
    hydroDx = hydroDx*hydro_nskip_x;
    hydroDz = hydroDz*hydro_nskip_z;
    
    glauberTarget = settings.word("Glauber:Target");
    glauberProjectile = settings.word("Glauber:Projectile");
    glauberImpactParam = settings.parm("Glauber:b");
    glauberIMax = settings.mode("Glauber:MaxInterpolationPoints");
    glauberEnvelope = settings.flag("Glauber:Envelope");
    
    PDFname = settings.word("PDF:LHAPDFset");
    PDFmember = settings.mode("PDF:LHAPDFmember");
    nuclearPDFname = settings.word("PDF:nuclearPDF");
    
    nbinFromFile = settings.mode("General:Nbin_from_File");
    string binary_filename = settings.word("General:Nbin_File_Name");
    //Either the name of the evolution file, or a beginning tag for the 
    //file for the case nbinFromFile == 1. -CFY
    evolution_name = settings.word("General:evolution_name");
    //The number of collisions in the file, counted up from zero
    cout << "evolution_name = " << evolution_name << endl;
    cout << "binary_filename = " << binary_filename << endl;
    background_file_name = settings.word("General:Background_file_name");
    cout << "background_file_name = " << background_file_name << endl;
    Nbin = 0;
  
    if(nbinFromFile == 2 or nbinFromFile == 3 or nbinFromFile == 4)
    {
        binary_info_ptr->init(binary_filename, nbinFromFile);
        //binary_info_ptr->print_info();
        //binary_info_ptr->check_samples();
    }
        
    ///read in data files with transition rates
    import->init(rateSelector);
    rates   = new Rates(rateSelector);
    if ( fixedTemperature == 0 and fixedEnergy == 0)
    {
          // initialize nuclei information
          glauber->init(inelasticXSec,glauberTarget,glauberProjectile,glauberImpactParam,glauberIMax,glauberEnvelope);
    }  

    //importLRates->init();
    
    //cout << "rate=" << importLRates->getRate(10., 5., 2.) << endl;
    
    // output of the rates into a file:
    //import->show_dGamma();
    
    /// initialize random number generator with current time as seed
    long long rnum;
    
    rnum=time(0)+seed;
    
    random->init_genrand64(rnum);
    
    // set random seed as current time
    srand(rnum);
    int now=rand()/100;
    stringstream nowstr;
    nowstr << "Random:seed = " << now;
    string timestring = nowstr.str();
    // put most usual PYTHIA setups in this if-statement
    // we only want this to be done in the production run
    // and not during a QGP brick + parton gun mode. 
    // RY June 9 2021
    if ( fixedTemperature == 0 and fixedEnergy == 0)
    {
        // PYTHIA settings:
        pythia.readString("Random:setSeed = on");

        pythia.readString(timestring);    
        //pythia.readString("Random:seed = 100");

        // finish after doing hard process
        pythia.readString("PartonLevel:all = on"); // off to only get hard process partons
        pythia.readString("HadronLevel:all = off");

        stringstream ptminstr;
        ptminstr << "PhaseSpace:pTHatMin =  " << jetpTmin;
        string ptminstring = ptminstr.str();
        pythia.readString(ptminstring);
        stringstream ptmaxstr;
        ptmaxstr << "PhaseSpace:pTHatMax =  " << jetpTmax;
        string ptmaxstring = ptmaxstr.str();
        pythia.readString(ptmaxstring);

        // very important to have this option!!
        pythia.readString("Check:event = off");
    }
    if(examineHQ == 1)
    {
        pythia.readString("HardQCD:gg2ccbar = on");
        pythia.readString("HardQCD:qqbar2ccbar = on");
        pythia.readString("HardQCD:gg2bbbar = on");
        pythia.readString("HardQCD:qqbar2bbbar = on");
    }

    // set the PDF and possible nuclear effects. note: if nuclear effects are chosen, the use of LHAPDF is enforced!
    if ( nuclearEffects!=0 )
    {
        if (!useLHAPDF)
            cout << "[MARTINI]:WARNING: using nuclear effects - turned on LHAPDF against initial settings." << endl;
        pythia.readString("PDF:nuclearEffects = 1");
        stringstream PDFnameStr;
        PDFnameStr << "PDF:pSet=LHAPDF6:" << PDFname << "/" << PDFmember;
        pythia.readString(PDFnameStr.str());
        stringstream nuclearPDFnameStr;
        nuclearPDFnameStr << "PDF:nuclearPDF=" << nuclearPDFname;
        pythia.readString(nuclearPDFnameStr.str());
        stringstream projAstr, targAstr;
        projAstr << "PDF:proj_atomicNumber = " << glauber->get_Projetile_A();
        targAstr << "PDF:targ_atomicNumber = " << glauber->get_Target_A();
        pythia.readString(projAstr.str());
        pythia.readString(targAstr.str());
    }
    else
    {
        if(useLHAPDF) 
        {
            stringstream PDFnameStr;
            PDFnameStr << "PDF:pSet=LHAPDF6:" << PDFname << "/" << PDFmember;
            pythia.readString(PDFnameStr.str());
        }
        pythia.readString("PDF:nuclearEffects = 0");
    }
  
    if ( fixedTemperature == 0 and fixedEnergy == 0)
        // init PYTHIA with given center of mass energy. will be changed for full AA collision.
        initPythia();

    // read the hydro data from the data file(s)
    if (evolution and fixedTemperature==0) 
    {
        hydroSetup->readHydroData(hydroTau0, hydroTauMax, hydroDtau, 
                                  hydroXmax, hydroZmax, hydroDx, hydroDz, 
                                  hydro_nskip_tau, hydro_nskip_x, hydro_nskip_z,
                                  hydroWhichHydro, hydroTfinal, tauEtaCoordinates, 
                                  evolution_name); 
        cout << "OK after readHydroData" << endl;
        // update hydroTauMax according to the actual size of the hydro medium
        hydroTauMax = hydroSetup->get_hydro_tau_max();
    }
    if (trackHistory)
    {
        Tt = new vector<double>; // stores the history of the trajectory: t, x(t), y(t), z(t), dE/dt(t), dpx/dt(t), dpy/dt(t), dpz/dt(t) *for one particle*
        Tx = new vector<double>;
        Ty = new vector<double>;
        Tz = new vector<double>;
        TE = new vector<double>;
        Tpx = new vector<double>;
        Tpy = new vector<double>;
        TdEdt = new vector<double>;
        Tdpxdt = new vector<double>;
        Tdpydt = new vector<double>;
        Tdpzdt = new vector<double>;
        hydroSetup->output_temperature_evolution("temp_evo");
    }

    return true;
}

//Sangyong's addition
// things to do when there is a daughter
int MARTINI::hasDaughter(vector<Parton> *plist, int i, double T)
{
    int daughter;
    double x, y, z, rperp, kx, ky, kz, kox, koy, koz, ko;
    double kcrx, kcry, kcrz, kperp;
    
    daughter = plist->at(i).daughter();
    
    kox = plist->at(i).p_at_split().px();
    koy = plist->at(i).p_at_split().py();
    koz = plist->at(i).p_at_split().pz();
    ko = sqrt(kox*kox + koy*koy + koz*koz);
    if (ko == 0.0)
    {
        kox = plist->at(daughter).p_at_split().px();
        koy = plist->at(daughter).p_at_split().py();
        koz = plist->at(daughter).p_at_split().pz();
        ko = sqrt(kox*kox + koy*koy + koz*koz);
        if (ko == 0.0)
        {
             cout << "hasDaughter: this can'happen.\n" << endl;
             exit(0);
        }
    }
 
     x = plist->at(i).x() - plist->at(daughter).x();
     y = plist->at(i).y() - plist->at(daughter).y();
     z = plist->at(i).z() - plist->at(daughter).z();

     // k_orig cross (delta x) = k_orig cross (delta xperp)
     // | k_orig cross (delta x)| = k_orig rperp

     kcrx = y*koz - z*koy;
     kcry = z*kox - x*koz;
     kcrz = x*koy - y*kox;
     rperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
     kx = plist->at(i).p().px() - plist->at(daughter).p().px();
     ky = plist->at(i).p().py() - plist->at(daughter).p().py();
     kz = plist->at(i).p().pz() - plist->at(daughter).p().pz();
    
     // k_orig cross (delta k) = k_orig cross (delta kperp)
     // | k_orig cross (delta k)| = k_orig kperp

     kcrx = ky*koz - kz*koy;
     kcry = kz*kox - kx*koz;
     kcrz = kx*koy - ky*kox;
     kperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
     double k = plist->at(daughter).p().pAbs();
     double Crperp = 0.25*pow(k/T, 0.11);
    // two lows : 0.05 and 0.01
    // two highs: 0.08 and 2000000(if I rememeber correctly)
     if ( rperp*kperp < Crperp*hbarc ) 
     {
        return 1; // do only elastic
     }
     else  // decouple
     {
        plist->at(i).mother_of(-1); 
        plist->at(i).daughter_of(-1); 
        plist->at(daughter).mother_of(-1); 
        plist->at(daughter).daughter_of(-1); 
    
        return 0; // do rad plus elastic
     }
}// hasDaughter;


//Sangyong's addition
// things to do when there is a mother
int MARTINI::hasMother(vector<Parton> *plist, int i, double T)
{
    int mother;
    double x, y, z, rperp, kx, ky, kz, kox, koy, koz, ko;
    double kcrx, kcry, kcrz, kperp;
    
    mother= plist->at(i).mother();
    
    kox = plist->at(mother).p_at_split().px();
    koy = plist->at(mother).p_at_split().py();
    koz = plist->at(mother).p_at_split().pz();
    ko = sqrt(kox*kox + koy*koy + koz*koz);
    if( ko == 0.0 )
    {
        kox = plist->at(i).p_at_split().px();
        koy = plist->at(i).p_at_split().py();
        koz = plist->at(i).p_at_split().pz();
        ko = sqrt(kox*kox + koy*koy + koz*koz);
        if( ko == 0.0 )
        {
            cout << "hasMother: this can'happen.\n" << endl;
            exit(0);
        }
    }// if ko == 0

    x = plist->at(i).x() - plist->at(mother).x();
    y = plist->at(i).y() - plist->at(mother).y();
    z = plist->at(i).z() - plist->at(mother).z();
    
    // k_orig cross (delta x) = k_orig cross (delta xperp)
    // | k_orig cross (delta x)| = k_orig rperp
    
    kcrx = y*koz - z*koy;
    kcry = z*kox - x*koz;
    kcrz = x*koy - y*kox;
    rperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
    kx = plist->at(i).p().px() - plist->at(mother).p().px();
    ky = plist->at(i).p().py() - plist->at(mother).p().py();
    kz = plist->at(i).p().pz() - plist->at(mother).p().pz();
    
    // k_orig cross (delta k) = k_orig cross (delta kperp)
    // | k_orig cross (delta k)| = k_orig kperp
    
    kcrx = ky*koz - kz*koy;
    kcry = kz*kox - kx*koz;
    kcrz = kx*koy - ky*kox;
    kperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
    double k = plist->at(i).p().pAbs();
    double Crperp = 0.25*pow(k/T, 0.11);
    if ( rperp*kperp < Crperp*hbarc ) 
    {
        return 1; // do only elastic
    }
    else  // decouple
    {
        plist->at(i).mother_of(-1); 
        plist->at(i).daughter_of(-1); 
        plist->at(mother).mother_of(-1); 
        plist->at(mother).daughter_of(-1); 
        
        return 0; // do rad plus elastic
    }

}// hasMother;

void MARTINI::output_tracked_history(string filename)
{
    cout << "MARTINI:: output tracked histroy into file: " << filename << endl;
    ofstream history(filename.c_str());
    for(int i = 0; i < Tt->size(); i++)
    {
        history << scientific << setw(16) << setprecision(6)
                << (*Tt)[i] << "  " << (*Tx)[i] << "  " << (*Ty)[i] << "  " << (*Tz)[i] << "  "
                << (*TE)[i] << "  " << (*Tpx)[i] << "  " << (*Tpy)[i] << endl;
    }
    history.close();
}
int MARTINI::generate_pGun(vector<Parton> *plist, int gunID)
{
    // RY June 9 2021
    //cout<<"MARTINI: generate a parton gun of type "<< gunID<<endl; 
    //cout<<"Ensure fixed temperature and energy:"<<endl;
    if ( fixedTemperature == 0 or fixedEnergy == 0)
    {
        cout<<"[ERROR]: parton gun mode must be at fixed E and T. Exit."<<endl;
        exit(-1);
    }
    Parton parton;
    double mass;
    if ( gunID >= 1 and gunID <= 3)//quark
    {
        parton.col(1000000);
        parton.acol(0);
    }
    else if ( gunID <= -1 and gunID >= -3 )//anti-quark
    {
        parton.col(0);
        parton.acol(1000000);
    }
    else if ( gunID == 21 ) //gluon
    {
        parton.col(1000000);
        parton.acol(1000001);
    }
    else
    {
        string test = abs(gunID) < 1 ? "True" : "False";
        cout<<"gunID: "<<gunID<<" abs(gunID) < 1 ? "<< test<<endl;
        //it's not a u, ubar, d, dbar, s, sbar or g so
        // no: say it's not allowed and terminate.
        cout<<"[ERROR]: parton ID is not recognized. Exit."<<endl;
        exit(-1);
    }
    if ( abs(gunID) == 1 or abs(gunID) == 2) // d or u quark
        mass = 0.33;
    else if ( abs(gunID) == 3) // s quark
        mass = 0.55;
    else // gluon
        mass = 0.0;

    Vec4 pvec;
    double energy = sqrt(initialPXjet*initialPXjet + initialPYjet*initialPYjet + mass*mass);
    pvec.px(initialPXjet);
    pvec.py(initialPYjet);
    pvec.pz(0.0);
    pvec.e(energy);
    parton.mass(mass);
    parton.p(pvec);
    parton.id(gunID); 
    parton.source(10);
    parton.status(16);
    parton.x(initialXjet);
    parton.y(initialYjet);
    parton.z(0.0);
    parton.tFinal(0.0);
    parton.mother_of(-1); // mother of nobody
    parton.daughter_of(-1); // daughter of nobody
    parton.p_at_split(parton.p()); // original momentum, don't need it but might as well populate it
    parton.init_counts();
    plist->push_back(parton); 
    return 1;
}
