// Rates.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains sampling routines for the radiative transition rates

#include "Rates.h"

Rates::Rates(int rateSelectorIn)
{
  rateSelector = rateSelectorIn;
}

double Rates::function(double x, double y, Import *import, int process)
{
  double rate = 0.;
  if ( process == 1 ) rate = import->getRate(x, y);
  if ( process == 2 ) rate = import->getRate_gqq(x, y);
  if ( process == 3 ) rate = import->getRate_ggg(x, y);
  if ( process == 4 ) rate = import->getRate_em(x, y);
  return rate;
}

// the rate is multiplied by g^4 T in the three functions (integrate*) below - the rate in the data file is missing this.
// approximation of the area under the transition rate from k=-12 T to k=p*T (where k<-0.05)
// fator 2 in "Gamma_gqq" comes from the fact that the fitting parameters are based on integrals for k \in (-inf, 0.5p].
// same for integratePos and integrate.
Norms Rates::integrateNeg(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  double fac = 4*PI*alpha_s*4*PI*alpha_s*T;

  if (rateSelector == 1){ //LO = old-Peter
    norms.Gamma     = (0.328271 + 0.0486746/p + 0.140053/p/p - 0.782095/p/p/p)*fac;
    norms.Gamma_gqq = (0.00148121/p + 0.00400277/p/p + 0.012097/p/p/p)*fac*Nf*2;
    norms.Gamma_ggg = (0.738553 + 0.134064/p - 0.147011/p/p + 1.85453/p/p/p)*fac;
  }else if (rateSelector == 2){ //NLO
    norms.Gamma     = (0.866713 + 0.119781/p + 0.376045/p/p - 2.08173/p/p/p)*fac;
    norms.Gamma_gqq = (0.00337617/p + 0.00804572/p/p + 0.0250539/p/p/p)*fac*Nf*2;
    norms.Gamma_ggg = (1.94996 + 0.331612/p - 0.510103/p/p + 4.98656/p/p/p)*fac;
  }else if (rateSelector == 3){ //NP
    norms.Gamma     = (0.734152 + 0.0991851/p + 0.321493/p/p - 1.76695/p/p/p)*fac;
    norms.Gamma_gqq = (0.00289406/p + 0.00671999/p/p + 0.0194533/p/p/p)*fac*Nf*2;
    norms.Gamma_ggg = (1.65172 + 0.276364/p - 0.468439/p/p + 4.2474/p/p/p)*fac;
  }else if (rateSelector == 4){ //LO-smallq = old-AMY
    norms.Gamma     = (0.375912 + 0.0555188/p + 0.160354/p/p - 0.895641/p/p/p)*fac;
    norms.Gamma_gqq = (0.00166884/p + 0.00447542/p/p + 0.0136848/p/p/p)*fac*Nf*2;
    norms.Gamma_ggg = (0.845735 + 0.152988/p - 0.169329/p/p + 2.1233/p/p/p)*fac;
  }
  return norms;
}

// approximation of the area under the transition rate from k=0.05 T to k=p*T (p is just dimensionless variable here)
Norms Rates::integratePos(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  double fac = 4*PI*alpha_s*4*PI*alpha_s*T;
  if (rateSelector == 1){ //LO = old-Peter
    norms.Gamma     = (0.43537 - 0.760946/p + 4.96699/p/p - 12.2663/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0262806/p - 0.103007/p/p + 0.0820234/p/p/p + (0.00545606 + 0.0762096/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (0.979113 - 1.71037/p + 10.8086/p/p - 23.1217/p/p/p)*fac;
  }else if (rateSelector == 2){ //NLO
    norms.Gamma     = (1.09691 - 1.23324/p + 7.58513/p/p - 19.4038/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0288709/p - 0.0898382/p/p + 0.0685581/p/p/p + (0.00808837 + 0.0704805/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (2.46714 - 2.79209/p + 15.9446/p/p - 31.6451/p/p/p)*fac;
  }else if (rateSelector == 3){ //NP
    norms.Gamma     = (0.920411 - 0.939677/p + 5.72427/p/p - 14.8005/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0199054/p - 0.0563209/p/p + 0.042419/p/p/p + (0.00615394 + 0.0452676/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (2.07022 - 2.14342/p + 11.936/p/p - 23.2134/p/p/p)*fac;
  }else if (rateSelector == 4){ //LO-smallq = old-AMY
    norms.Gamma     = (0.496507 - 0.836112/p + 5.42924/p/p - 13.4397/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0281335/p - 0.109512/p/p + 0.0872221/p/p/p + (0.00594387 + 0.0810782/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (1.11663 - 1.87851/p + 11.7854/p/p - 25.0874/p/p/p)*fac;
  }
  return norms;
}

// p is passed in the unit of T to this function
Norms Rates::integrate(double p, double T, double alpha_s, int Nf)
{
  Norms norms;

  double fac = pow(4*PI*alpha_s, 2.)*T;
  if (rateSelector == 1){ //LO = old-Peter
    norms.Gamma     = (0.763641 - 0.712272/p + 5.10704/p/p - 13.0484/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0247993/p - 0.099004/p/p + 0.0941204/p/p/p + (0.00545606 + 0.0762096/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (1.71767 - 1.57631/p + 10.6615/p/p - 21.2672/p/p/p)*fac;
    norms.Gamma_em  = (-0.150335/p - 0.589552/p/p + 0.465951/p/p/p + (0.0540207 + 0.42001/p)/sqrt(p))*fac*((1./137.)/alpha_s);
  }else if (rateSelector == 2){ //NLO
    norms.Gamma     = (1.96362 - 1.11346/p + 7.96118/p/p - 21.4856/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0254948/p - 0.0817925/p/p + 0.093612/p/p/p + (0.00808837 + 0.0704805/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (4.41709 - 2.46048/p + 15.4345/p/p - 26.6585/p/p/p)*fac;
    norms.Gamma_em  = (-0.139714/p - 0.496187/p/p + 0.379805/p/p/p + (0.0748676 + 0.352413/p)/sqrt(p))*fac*((1./137.)/alpha_s);
  }else if (rateSelector == 3){ //NP
    norms.Gamma     = (1.65456 - 0.840492/p + 6.04577/p/p - 16.5675/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0170113/p - 0.0496009/p/p + 0.0618723/p/p/p + (0.00615394 + 0.0452676/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (3.72194 - 1.86706/p + 11.4675/p/p - 18.966/p/p/p)*fac;
    norms.Gamma_em  = (-0.109725/p - 0.449801/p/p + 0.359923/p/p/p + (0.0573986 + 0.306699/p)/sqrt(p))*fac*((1./137.)/alpha_s);
  }else if (rateSelector == 4){ //LO-smallq = old-AMY
    norms.Gamma     = (0.872419 - 0.780593/p + 5.58959/p/p - 14.3353/p/p/p)*fac;
    norms.Gamma_gqq = (-0.0264646/p - 0.105037/p/p + 0.100907/p/p/p + (0.00594387 + 0.0810782/p)/sqrt(p))*fac*Nf*2;
    norms.Gamma_ggg = (1.96236 - 1.72552/p + 11.6161/p/p - 22.9641/p/p/p)*fac;
    norms.Gamma_em  = (-0.162315/p - 0.640589/p/p + 0.507009/p/p/p + (0.0587912 + 0.455284/p)/sqrt(p))*fac*((1./137.)/alpha_s);
  }
  return norms;
}


/// q->g conversion rate                                                                                                       
double Rates::Gammaqg(double En, double T, double alpha_s)
{
  return 4./3.*2.*PI*alpha_s*alpha_s*T*T/(3.*En)*(0.5*log(En*T/mqs(T,alpha_s))-0.36149);
}

/// q->gamma conversion rate    
// difference is one less C_f=4/3, and an alpha_e=1/137 instead of alpha_s            
// has to be multiplied later with the ratio (e_f/e)^2                                                  
double Rates::Gammaqgamma(double En, double T, double alpha_s)
{   
    // This is the collinear conversion approximation. Original
    // return value was rate in the limit En/T --> Infty
    //return 2.*PI*alpha_s*(1./137.)*T*T/(3.*En)*(0.5*log(En*T/mqs(T,alpha_s))-0.36149);
    // RY June 9 2021
    double c = 0.041*(T/En) - 0.3615 + 1.01*exp(-1.35*En/T);
    return 2*PI*alpha_s*(1./137)*T*T*(0.5*log(En*T/mqs(T, alpha_s)) + c)/(3*En); 
}

/// g->q conversion rate                                                                                                       
double Rates::Gammagq(double En, double T, double alpha_s, int Nf)
{
  return Nf*3./8.*4./3.*2.*PI*alpha_s*alpha_s*T*T/(3.*En)*(0.5*log(En*T/mqs(T,alpha_s))-0.36149);
}

// calculates the area under the envelope function when using the rejection method
// (integrals had been solved analytically before)
// no-longler useful, but let's keep this function
double Rates::area (double y, double u, int posNegSwitch, int process, Import * import)
{
    if (process==1)
    {
        if (posNegSwitch==1)
	      return(0.5299 - 0.025/y + 0.01*log(y));
        else 
	      return(-0.0020833-0.025/y);
    }
    else if (process==2)
    {
        if (posNegSwitch==1)
	      return(1.2*function(u,0.05,import,2)*(y-0.05));
        else 
	      return( (6.877797*pow(10.,-8.)*exp(14.5/u)
		        -0.008805468*exp((2.5-y+0.98*u*y)/u))/(1.0204082-u));
    }      
    else if (process==3)
    {
        if (posNegSwitch==1)
	      return(2.05991 - 0.1/y + 0.02*log(y));
        else 
	      return(-0.00833333 - 0.1/y);
    }      
    else if (process==4) // radiating gamma
    {
        return((0.0333333*pow(y,0.3))/pow(u,0.5));
    }
    return(0.0);
}
ReturnValue Rates::findValuesRejection(double u, double T, double alpha_s, int Nf,
                       Random *random, Import *import, int process)
{
    ReturnValue f;
    Norms Pos, Neg;
    Pos = integratePos(u,T,alpha_s, Nf);
    Neg = integrateNeg(u,T,alpha_s, Nf);
    
    /* 
        process == 1 : quark radiating gluon
       process == 2 : gluon split into quark-antiquark pair 
       process == 3 : gluon radiating gluon
       process == 4 : quark radiating photon */
    
    // this switch will hold the decision whether k is positive or negative:
    // 0 : negative, 1 : positive
    int posNegSwitch = 1;
    double neg_ratio = 0;
    if      (process == 1) {neg_ratio = Neg.Gamma/(Neg.Gamma+Pos.Gamma);}
    else if (process == 2) {neg_ratio = Neg.Gamma_gqq/(Neg.Gamma_gqq+Pos.Gamma_gqq);}
    else if (process == 3) {neg_ratio = Neg.Gamma_ggg/(Neg.Gamma_ggg+Pos.Gamma_ggg);}
    if (random->genrand64_real1() < neg_ratio){posNegSwitch=0;}
    double y_low_cut = 0.05;
    double y_high_cut = u+12;
    
    double rndmax = 1.5;
    //if (process==4){rndmax=10;}
     double randA, y, fy, fyAct;
    if (process==2)
    { // g -> q + qbar is flat
        y_high_cut = 0.5*u;
        if (posNegSwitch == 1) 
        {// if k > 0
            do
            {
                y = y_low_cut + (y_high_cut-y_low_cut)*random->genrand64_real1();
                fyAct = function(u, y, import, process);
                if (fyAct > rndmax) 
                {
                    cerr << "envelope not enough:\t" <<fyAct << endl;
                }
            } while (rndmax*random->genrand64_real1() > fyAct); 
            f.y = y;
        }
        else
        {// if k < 0
            do
            {
                y = -y_low_cut - (12.-y_low_cut)*random->genrand64_real1();
                fy = 1.;
                fyAct = function(u, y, import, process);
                if (fyAct > rndmax) {cerr << "envelope not enough:\t" <<fyAct << endl;}
            } while (rndmax*random->genrand64_real1() > fyAct); 
            f.y = y; 
        }
    }
    else if (process==4)
    {// q -> q + gamma
        //y_high_cut = 0.5*u;
        // The following is commented out for envelope function issues.
        // Fix done by Shuzhe Shi, August 6 2021
        //do
        //{
        //    randA = random->genrand64_real1() * (log(y_high_cut/y_low_cut));
        //    y = exp(log(y_low_cut) + randA);
        //    fyAct = function(u, y, import, process);
        //    fy = 1./y;
        //    if (fyAct/fy > rndmax) {cerr << "envelope not enough:\t" <<fyAct << endl;}
        //    count++;
        //    if (count%10000)
        //    {
        //        cout<<"count within findValuesRejection : "<< count<<endl;
        //        cout<<" u : "<<u<<" T : "<<T<<endl;
        //        cout<<" fyAct/fy: "<<fyAct/fy<<endl;
        //        exit(-1);
        //    }
        //} while (rndmax*random->genrand64_real1() > fyAct/fy); 
           //y_high_cut = 0.5*u;
        int count = 0;
        do
        {
            randA = random->genrand64_real1() * 2*(sqrt(y_high_cut)-sqrt(y_low_cut));
            y = pow(sqrt(y_low_cut) + 0.5*randA, 2);
            fyAct = function(u, y, import, process);
            fy = 1./sqrt(y);
            count++;
            if (count > 2000)
            {
                cout << "Tried to find a k more than 2000 times, abort"<<endl;
                f.y = 0;
                return f;
            }
            //if (fyAct/fy > rndmax) {cerr << "envelope not enough:\t" <<fyAct << endl;} this might happen for the rare cases that p>>T and k/p~1.
        } while (rndmax*random->genrand64_real1() > fyAct/fy); 
        f.y = y;
    } 
    else if (process<5)
    {
        if (process==3){y_high_cut = 0.5*u;}// g -> gg
        if (posNegSwitch == 1)
        {// if k > 0
            do
            {
                randA = random->genrand64_real1() * (log(y_high_cut/y_low_cut)+1./y_low_cut - 1./y_high_cut);
                y = 1./gsl_sf_lambert_W0(exp(1./y_low_cut-log(y_low_cut)-randA));
                fy = 1./y/y + 1./y;
                fyAct = function(u, y, import, process);
                if (fyAct/fy > rndmax) {cerr << "envelope not enough:\t" <<fyAct/fy << endl;}
            } while (rndmax*random->genrand64_real1() > fyAct/fy); 
            f.y = y;
        }
        else
        {// if k < 0
            do
            {
               randA = random->genrand64_real1() * (1./y_low_cut - 1./12.);
               y = -1./(1./y_low_cut-randA);
               fy = 1./y/y;
               fyAct = function(u, y, import, process);
               if (fyAct/fy > rndmax) {cerr << "envelope not enough:\t" <<fyAct/fy << endl;}
            } while (rndmax*random->genrand64_real1() > fyAct/fy); 
            f.y = y; 
        }
    }
    else
    {
        cerr << "Invalid process number (" << process << ")" << endl;
    }
    return f;
    }

