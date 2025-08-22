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
  double rate;
  if ( process == 1 ) rate = import->getRate(x, y);
  if ( process == 2 ) rate = import->getRate_gqq(x, y);
  if ( process == 3 ) rate = import->getRate_ggg(x, y);
  if ( process == 4 ) rate = import->getRate_em(x, y);
  return rate;
}

// the rate is multiplied by g^4 T in the three functions (integrate*) below - the rate in the data file is missing this.
// approximation of the area under the transition rate from k=-12 T to k=p*T (where k<-0.05)
Norms Rates::integrateNeg(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  if (rateSelector == 1)
    {
      norms.Gamma = (0.376329 - 0.774118/(p*p) + 0.558981/pow(p,1.5) - 0.0639996/p + 0.00903106/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.0365634/(p*p*p) + 0.0119574/(p*p) + 0.0044576/p)*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_ggg = (0.846909 + 2.16166/(p*p*p) - 0.182009/(p*p) + 0.154982/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 2 )
    {
      norms.Gamma = (0.00503211 - 0.00212099/(p*p) + 0.00171696/p - 0.0000904348/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.00211772/(p*p*p)) + 0.00574703/(p*p) + 0.00088597/p + 0.0000173461/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_ggg = (0.0113179 - 0.0122324/(p*p*p) + 0.0255511/(p*p) + 0.00293152/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 3)
    {
      norms.Gamma = (0.00430578 - 0.001845/(p*p) + 0.00142608/p - 0.0000795888/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.0643767/pow(p,5.)) - 0.00303838/(p*p*p*p) + 0.00896399/(p*p*p) + 0.00298915/(p*p) 
			    + 0.000960497/p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.00968427 - 0.0132583/(p*p*p) + 0.0227694/(p*p) + 0.00237067/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 4)
    {
      norms.Gamma = (0.00441148 - 0.00188564/(p*p) + 0.00146807/p - 0.0000811793/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.0654783/pow(p,5.)) - 0.0029149/pow(p,4.) + 0.0091532/(p*p*p) + 0.00304244/(p*p) + 0.000975992/p)
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.00992202 - 0.0131448/(p*p*p) + 0.0231876/(p*p) + 0.00245127/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 5)
    {
      norms.Gamma = (0.00274605 - 0.00121989/(p*p) + 0.000826588/p - 0.000055282/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.00326374/(p*p*p)) + 0.00382808/(p*p) + 0.000577548/p + 0.0000122508/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.00617612 - 0.0131497/(p*p*p) + 0.0159918/(p*p) + 0.0012429/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 6)
    {
      norms.Gamma = (0.00211528 - 0.000950252/(p*p) + 0.000600056/p - 0.0000447104/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.186678/pow(p,6.)) + 0.0114748/pow(p,4.) + 0.00203147/(p*p) + 0.000571453/p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.00475738 - 0.0119168/(p*p*p) + 0.0128732/(p*p) + 0.000836708/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 7)
    {
      norms.Gamma = (0.00221685 - 0.000994417/(p*p) + 0.000635659/p - 0.0000464563/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.195702/pow(p,6.)) + 0.0122506/pow(p,4.) + 0.00211997/(p*p) + 0.000592027/p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.00498584 - 0.0121715/(p*p*p) + 0.0133924/(p*p) + 0.000899338/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 8 )
    {
      norms.Gamma = (0.3292194797883 - 0.6758562562413/(p*p) + 0.4871465960062/pow(p,1.5) - 0.0539320881756/p + 0.0078775876607/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*((0.0321479979527/(p*p*p)) + 0.0141918292508/(p*p) + 0.0043376733045/p - 0.0000124591578/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_ggg = (0.7408894524117 + 1.860773318657/(p*p*p) - 0.1352673301036 /(p*p) + 0.1400730913508/p)*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  return norms;
}

// approximation of the area under the transition rate from k=0.05 T to k=p*T (p is just dimensionless variable here)
Norms Rates::integratePos(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  if (rateSelector == 1)
    {
      norms.Gamma = (0.501144 - 0.455219/(p*p) + 0.147232/p - 0.171267/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-0.0000219797 - 0.0293002/(p*p) + 0.0404452/pow(p,1.8) - 0.016909/p + 0.00754396/pow(p,0.4))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (1.12664 + 1.91729/(p*p*p) + 0.119407/p - 0.354566/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 2)
    {
      norms.Gamma = (0.0494571 - 0.166125/(p*p) + 0.192479/p - 0.153128/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-0.000428053 + 0.0162373/(p*p) - 0.0114166/p + 0.00041804/pow(p,0.8) + 0.00343915/pow(p,0.2) 
			    + 2.59233*pow(10.,-9.)*p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.11091 - 0.606136/(p*p*p) + 0.385639/p - 0.340793/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 3)
    {
      norms.Gamma = (0.0470191 - 0.168738/(p*p) + 0.196352/p - 0.152618/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.0916944/pow(p,1.1) - 0.0961707/p + 0.00814928/pow(p,0.4) - 0.000212637/pow(p,0.2))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.105425 - 0.646291/(p*p*p) + 0.393123/p - 0.339446/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 4)
    {
      norms.Gamma = (0.0473874 - 0.16828/(p*p) + 0.195736/p - 0.152697/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.09135/pow(p,1.1) - 0.0958814/p + 0.00817075/pow(p,0.4) - 0.000216607/pow(p,0.2))*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.106254 - 0.639909/(p*p*p) + 0.391934/p - 0.339657/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 5)
    {
      norms.Gamma = (0.00274605 - 0.00121989/(p*p) + 0.000826588/p - 0.000055282/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-0.000451316 + 0.012578/(p*p) - 0.000829902/p - 0.00765827/pow(p,0.8) + 0.00358473/pow(p,0.2)
			    + 2.76734*pow(10.,-9.)*p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.0914979 - 0.777253/(p*p*p) + 0.416761/p - 0.335444/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 6)
    {
      norms.Gamma = (0.0377103 - 0.191536/(p*p) + 0.216734/p - 0.150216/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-0.0000605083 + 0.00211754/(p*p) + 0.0166933/p - 0.0601028/pow(p,0.6) + 0.0443532/sqrt(p) 
			    + 7.54605*pow(10.,-10.)*p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.0844763 - 0.867009/(p*p*p) + 0.431646/p - 0.332926/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 7)
    {
      norms.Gamma = (0.038249 - 0.18944/(p*p) + 0.21526/p - 0.150384/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-0.0000594666 + 0.00252652/(p*p) + 0.0160471/p - 0.0593556/pow(p,0.6) + 0.0439917/sqrt(p) 
			    + 7.42279*pow(10.,-10.)*p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.085689 - 0.850046/(p*p*p) + 0.428924/p - 0.333395/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
    }
  else if (rateSelector == 8)
    {
      norms.Gamma = (0.5321935162467 - 3.103687009651/(p*p) + 2.013878510392/p - 0.9416445905935/sqrt(p))
											*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.00046561236634 - 0.046213050541/(p*p) + 0.0999017925532/p - 0.0817080199182/pow(p,0.8) 
											+ 0.0080898702818/pow(p,0.2) - 1.252535698671*pow(10.,-8.)*p)*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (1.192331601118 - 11.52502882017/(p*p*p) + 3.30100321574/p - 1.90491167553/sqrt(p))
											*4*PI*alpha_s*4*PI*alpha_s*T;
    }

  return norms;
}

// p is passed in units of T to this function
Norms Rates::integrate(double p, double T, double alpha_s, int Nf, int Ldependence, double L)
{
  Norms norms;
  double L_fac;
  double fac = 4*PI*alpha_s*4*PI*alpha_s*T;
  if (rateSelector == 1)
    {
      norms.Gamma = (0.8777175507615367 - 0.669043368107278/(p*p) + 0.2556218941655182/p - 0.1790856284715404/sqrt(p))*fac;
      norms.Gamma_gqq = (Nf*(0.076727/(p*p*p) - 0.0101059/(p*p) + 0.0413716/p - 0.0495917/pow(p,0.8) + 0.0189409/sqrt(p)))*fac;
      norms.Gamma_ggg = (1.97443 + 8.58545/(p*p*p) - 2.12458/(p*p) + 0.569455/p - 0.401423/sqrt(p))*fac;
      // total rate for q->qgamma
      // has to be multiplied later with the ratio (e_f/e)^2 (4/9) for up, (1/9) for down or strange
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000148478 + 0.0629722/(p*p*p) - 0.0479871/p + 0.0489135/sqrt(p)); 
    }
  else if (rateSelector == 2)
    {
      norms.Gamma = (0.0544893 - 0.168246/(p*p) + 0.194196/p - 0.153219/sqrt(p))*fac;
      norms.Gamma_gqq = (Nf*(0.00866093/(p*p*p) - 0.0000364348/(p*p) + 0.0476482/p - 0.0631694/pow(p,0.8) + 0.0196501/sqrt(p)))*fac;
      norms.Gamma_ggg = (0.12277 + 2.34586/(p*p*p) - 1.24584/(p*p) + 0.57925/p - 0.370621/sqrt(p))*fac;
      // total rate for q->qgamma
      // has to be multiplied later with the ratio (e_f/e)^2 (4/9) for up, (1/9) for down or strange
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000130102 + 0.189901/(p*p*p) - 0.0903018/p + 0.0506507/sqrt(p)); 
    }
  else if (rateSelector == 3)
    {
      norms.Gamma = (0.0513249 - 0.170583/(p*p) + 0.197778/p - 0.152697/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.00579244/(p*p*p) + 0.0000208018/(p*p) + 0.0494433/p - 0.0648592/pow(p,0.8) + 0.0197331/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.115651 + 2.30562/(p*p*p) - 1.249/(p*p) + 0.586224/p - 0.369281/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000134386 + 0.213977/(p*p*p) - 0.093559/p + 0.0505387/sqrt(p)); 
    }
  else if (rateSelector == 4)
    {
      norms.Gamma = (0.0517989 - 0.170166/(p*p) + 0.197204/p - 0.152778/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(0.00622666/pow(p,3.) + 0.0000197786/(p*p) + 0.0491497/p - 0.0645873/pow(p,0.8) + 0.0197198/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.116717 + 2.31195/(p*p*p) - 1.24852/(p*p) + 0.585106/p - 0.369489/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000133608 + 0.21013/(p*p*p) - 0.0930606/p + 0.0505616/sqrt(p)); 
    }
  else if (rateSelector == 5)
    {
      norms.Gamma = (0.0435759 - 0.182181/(p*p) + 0.209574/p - 0.151182/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.00135664/(p*p*p)) - 0.000668446/(p*p) + 0.0559352/p - 0.0705371/pow(p,0.8) + 0.0200057/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.0982183 + 2.18907/(p*p*p) - 1.26188/(p*p) + 0.609636/p - 0.365416/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000157371 + 0.289558/(p*p*p) - 0.101865/p + 0.0496785/sqrt(p)); 
    }
  else if (rateSelector == 6)
    {
      norms.Gamma = (0.0398256 - 0.192486/(p*p) + 0.217334/p - 0.150261/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.00443762/(p*p*p)) - 0.00197762/(p*p) + 0.0607314/p - 0.0744365/pow(p,0.8) 
			    + 0.0201883/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.0844763 - 0.867009/(p*p*p) + 0.431646/p - 0.332926/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000056766 - 0.761379/pow(p,4.) + 0.243627/(p*p) - 0.146482/p + 0.0553315/sqrt(p)); 
    }
  else if (rateSelector == 7)
    {
      norms.Gamma = (0.0404659 - 0.190434/(p*p) + 0.215896/p - 0.15043/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(-(0.00396437/(p*p*p)) - 0.00167442/(p*p) + 0.0598072/p - 0.0737001/pow(p,0.8) + 0.0201541/sqrt(p))
	*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (0.0912236 + 2.1422/(p*p*p) - 1.27517/(p*p) + 0.623055/p - 0.363618/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.000057087 - 0.720786/pow(p,4.) + 0.234903/(p*p) - 0.144368/p + 0.0552769/sqrt(p)); 
    }
  else if (rateSelector == 8)
    {
      norms.Gamma = (0.8616261061951 - 3.291252289008/(p*p) + 2.110181881549/p - 0.9484514504401/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_gqq = Nf*(2.582964577842 /(p*p*p) - 1.701017553943/(p*p) + 1.497695542213/p - 1.196064197449/pow(p,0.8) 
												+ 0.1806660194562/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T; 
      norms.Gamma_ggg = (1.94632330923 + 61.78560172928/(p*p*p) - 30.78773357571/(p*p) + 8.040896793632/p 
												- 2.624946207817/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
      norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T*(0.0053055996096 + 2.327880611991/pow(p,3.) - 0.6675641452702/p 
												+ 0.3222533793572/sqrt(p)); 
    }
  
  if (Ldependence == 1) // fit to Simon's formation time dependent rates (beware: may be quite inaccurate)
    {
      if (pow(L,0.61)*4.6>0.2/sqrt(p)-pow(p*T,0.26)/pow(T,0.8))
	L_fac = pow(L,0.61)*4.6/pow(p*T,0.26)*pow(T,0.8);
      else
	L_fac = 1.;
      if (L_fac>1.) L_fac = 1.;
      //cout << L << endl;
      norms.Gamma     *= L_fac; 
      norms.Gamma_ggg *= L_fac; 
      norms.Gamma_gqq *= L_fac; 
      norms.Gamma_em  *= L_fac; 
    }
  else if (Ldependence == 2) // take the formation time dependent rates and introduce double the L dependence by hand
    {
      //if (pow((L-0.4),0.61*2.)*4.6>0.2/sqrt(p)-pow(p*T,0.26)/pow(T,0.8))
      L_fac = pow((L-0.4),0.61)*4.6/pow(p*T,0.26)*pow(T,0.8);
      if (L_fac>2.) L_fac = 2.;
      //cout << L << endl;
      norms.Gamma     *= L_fac; 
      norms.Gamma_ggg *= L_fac; 
      norms.Gamma_gqq *= L_fac; 
      norms.Gamma_em  *= L_fac; 
    }
  return norms;
}

// approximation of the area under the transition rate from k=-12 T to k=p*T (where k<-0.05)
Norms Rates::integrateNegOld(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  norms.Gamma = (0.3825932+0.03217262/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_gqq = Nf*(0.05749489/(p*p)-0.03112226/pow(p,1.8)+0.00445603/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_ggg = (0.85739544+0.2125156/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  return norms;
}

// approximation of the area under the transition rate from k=0.05 T to k=p*T (p is just dimensionless variable here)
Norms Rates::integratePosOld(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  norms.Gamma = (0.49384-0.463155/(p*p)+0.15763/p-0.16124/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_gqq = Nf*(0.04392766/(p*p)-0.05353506/p+0.03259168/pow(p,0.8)+0.00191753/pow(p,0.2))
    *4*PI*alpha_s*4*PI*alpha_s*T; 
  norms.Gamma_ggg = (1.1104368+1.433124/(p*p*p)+0.16666408/p-0.33858792/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
  return norms;
}

Norms Rates::integrateOld(double p, double T, double alpha_s, int Nf, int Ldependence, double L)
{
  Norms norms;
  double L_fac;
  double fac = 4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma = ((0.49384-0.463155/(p*p)+0.15763/p-0.16124/sqrt(p)) + (0.3825932+0.03217262/p))*fac;
  norms.Gamma_gqq = (Nf*(0.04392766/(p*p)-0.05353506/p+0.03259168/pow(p,0.8)+0.00191753/pow(p,0.2))
    + Nf*(0.05749489/(p*p)-0.03112226/pow(p,1.8)+0.00445603/p))*fac;
  norms.Gamma_ggg = ((1.1104368+1.433124/(p*p*p)+0.16666408/p-0.33858792/sqrt(p)) + (0.85739544+0.2125156/p))*fac;
  // total rate for q->qgamma
  // has to be multiplied later with the ratio (e_f/e)^2 (4/9) for up, (1/9) for down or strange
  norms.Gamma_em = 4*PI*alpha_s*4*PI*(1./137.)*T
    *(0.0000608043 - 0.0571554/(p*p) + 0.154845/pow(p,1.2) - 0.160349/p + 0.0558144/pow(p,0.5)); 
  if (Ldependence == 1)
    {
      if (pow(L,0.61)*4.6>0.2/sqrt(p)-pow(p*T,0.26)/pow(T,0.8))
	L_fac = pow(L,0.61)*4.6/pow(p*T,0.26)*pow(T,0.8);
      else
	L_fac = 1.;
      if (L_fac>1.) L_fac = 1.;
      //cout << L << endl;
      norms.Gamma     *= L_fac; 
      norms.Gamma_ggg *= L_fac; 
      norms.Gamma_gqq *= L_fac; 
      norms.Gamma_em  *= L_fac; 
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
  return 2.*PI*alpha_s*(1./137.)*T*T/(3.*En)*(0.5*log(En*T/mqs(T,alpha_s))-0.36149);
}

/// g->q conversion rate                                                                                                       
double Rates::Gammagq(double En, double T, double alpha_s, int Nf)
{
  return Nf*3./8.*4./3.*2.*PI*alpha_s*alpha_s*T*T/(3.*En)*(0.5*log(En*T/mqs(T,alpha_s))-0.36149);
}

// calculates the area under the envelope function when using the rejection method
// (integrals had been solved analytically before)
double Rates::area (double y, double u, int posNegSwitch, int process, Import * import)
{
  if (process==1)
    {
      if (posNegSwitch==1)
	return 0.5299 - 0.025/y + 0.01*log(y);
      else 
	return -0.0020833-0.025/y;
    }
  else if (process==2)
    {
      if (posNegSwitch==1)
	return 1.2*function(u,0.05,import,2)*(y-0.05);
      else 
	return (6.877797*pow(10.,-8.)*exp(14.5/u)
		-0.008805468*exp((2.5-y+0.98*u*y)/u))/(1.0204082-u);
    }      
  else if (process==3)
    {
      if (posNegSwitch==1)
	return 2.05991 - 0.1/y + 0.02*log(y);
      else 
	return -0.00833333 - 0.1/y;
    }      
  else if (process==4) // radiating gamma
    {
      return (0.0333333*pow(y,0.3))/pow(u,0.5);
    }      
}

ReturnValue Rates::findValuesRejection(double u, double T, double alpha_s, int Nf,
				       Random *random, Import *import, int process)
{
  ReturnValue f;
  Norms Pos, Neg;
  double x;
  double y;
  Pos = integratePos(u,T,alpha_s, Nf);
  Neg = integrateNeg(u,T,alpha_s, Nf);
  
  // this switch will hold the decision whether k is positive or negative: 0=negative, 1=positive
  int posNegSwitch = 1; 
  
  // set counters for rejections
  //int i = -1;
  //int j = 0;

  if (process==1) // for quark radiating gluon
    {
      // decide whether k shall be positive or negative (if x (uniform on [0,1]) < area(k<0)/area(all k) then k<0)
      if (random->genrand64_real1()<Neg.Gamma/(Neg.Gamma+Pos.Gamma)) posNegSwitch = 0;

      // generate random number according to envelope distribution f(y)
      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      y = 2.5/(gsl_sf_lambert_W0(2.59235*pow(10.,23.)
					 *exp(-100*random->genrand64_real1()*area(u+12.,u,posNegSwitch,1,import))));
	      // here random->genrand64_real1()*area(u+12.,posNegSwitch) 
	      // is a uniform random number on [0, area under f(x)]
	      x = random->genrand64_real1();
	      // x is uniform on [0,1]
	      //i++;
	    } while( x > (function(u,y,import,1))/((0.025/(y*y))+0.01/y)); 
	  // reject if x is larger than the ratio p(y)/f(y), f(y)=0.025/(y*y)+0.01/y
	  f.y=y;
	}
      else // if k < 0
	{
	  do
	    {
	      y = -12./(1+480.*random->genrand64_real1()*area(-0.05,u,posNegSwitch,1,import));
	      // here random->genrand64_real1()*area(u+12.,posNegSwitch) 
	      // is a uniform random number on [0, area under f(x)]
	      x = random->genrand64_real1();
	      // x is uniform on [0,1]
	      //i++;
	    } while( x > (function(u,y,import,1))/((0.025/(y*y)))); 
          // reject if x is larger than the ratio p(y)/f(y), f(y)=0.025/(y*y)
	  f.y=y;
	}
    }
  else if ( process == 2 ) // for gluon radiating quark - anti-quark pair
    {
      if (random->genrand64_real1()<Neg.Gamma_gqq/(Neg.Gamma_gqq+Pos.Gamma_gqq)) posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      y = 0.83333*(0.06*function(u,0.05,import,2)
			   +area(u/2.,u,posNegSwitch,2,import)*random->genrand64_real1())/function(u,0.05,import,2);
	      x = random->genrand64_real1();
	      //i++;
	    }while( x > (function(u,y,import,2))/(1.2*function(u,0.05,import,2))); 
	  // f(y)=1.2*function(u,0.05,import,2) (constant)
	  f.y=y;
	}
      else // if k < 0
	{
	  do
	    {
	      y = (2.5-u*log(7.81082*pow(10.,-6.)
			     *exp(14.5/u)+(-115.883+113.566*u) 
			     *area(-0.05,u,posNegSwitch,2,import)
			     *random->genrand64_real1()))/(1.-0.98*u);
	      x = random->genrand64_real1();
	      //if ((function(u,y,import,2))/(0.1*exp((0.98-1./u)*(-2.5+y))/u)>1.) cout << "WARNING..y_neg=" << y << endl;
	      //i++;
	    } while( x > (function(u,y,import,2))/(0.98*exp((1.-1./u)*(-2.5+y))/u) ); 
	  // f(y)=0.1*exp((0.98-1./u)*(-2.5+y))/u 
	  f.y=y;
	}
    }
  else   if (process==3) // for gluon radiating gluons
    {
      if (random->genrand64_real1()<Neg.Gamma_ggg/(Neg.Gamma_ggg+Pos.Gamma_ggg)) posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      y = 5./(gsl_sf_lambert_W0(2.68812*pow(10.,45.)
					 *exp(-50*random->genrand64_real1()*area(u/2.,u,posNegSwitch,3,import))));
	      x = random->genrand64_real1();
	      //if (y>u/2.) cout << " large y=" << y << endl;
	      //i++;
	    } while( x > (function(u,y,import,3))/((0.1/(y*y))+0.02/y)); 
	  // f(y)=0.1/(y*y)+0.02/y
	  f.y=y;
	}
      else // if k < 0
	{
	  do
	    {
	      y = -12./(1. + 120.*random->genrand64_real1()*area(-0.05,u,posNegSwitch,3,import));
	      x = random->genrand64_real1();
	      //i++;
	    } while( x > (function(u,y,import,3))/((0.1/(y*y)))); 
	  // f(y)=0.1/(y*y)
	  f.y=y;
	}
    }
  else if (process==4) // for quark radiating gamma
    {
      // generate random number according to envelope distribution f(y)
      do
	{
	  y = 83895.3*pow(pow(u,0.5)*random->genrand64_real1()*area(1.15*u,u,posNegSwitch,4,import),3.333333333333333);
	  // here random->genrand64_real1()*area(1.15*u,posNegSwitch) 
	  // is a uniform random number on [0, area under f(x)]
	  x = random->genrand64_real1();
	  // x is uniform on [0,1]
	  //i++;
	} while( x > (function(u,y,import,4))/((0.01/(pow(y,0.7)))/pow(u,0.5))); 
      // reject if x is larger than the ratio p(y)/f(y), f(y)=(0.01/pow(y,0.7))/pow(u,0.5)
      f.y=y;
    }
  //f.rejections=i;
  //f.acceptances=j;
  return f;
}
