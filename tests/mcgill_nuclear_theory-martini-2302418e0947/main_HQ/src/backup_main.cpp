#include<stdio.h>
#include<string>
#include<time.h>
#include<sys/timeb.h>
#include<cmath>
#include<list>
#include<gsl/gsl_sf_lambert.h>
#include "Constants.h"
#include "Random.h"
#include "Import.h"
#include "Pythia.h"
#include "Parton.h"

using namespace Pythia8;

struct ReturnValue 
{
  double x;
  double y;
  int rejections;
  int acceptances;
};

struct Norms 
{
  double Gamma;
  double Gamma_gqq;
  double Gamma_ggg;
};

struct HardPartons 
{
  Parton initParton[2];
};

double function(double x, double y, Import *import, int process)
{
  double rate;
  if ( process == 1 ) rate = import->getRate(x, y);
  if ( process == 2 ) rate = import->getRate_gqq(x, y);
  if ( process == 3 ) rate = import->getRate_ggg(x, y);
  return rate;
}

// approximation of the area under the transition rate from k=0.05 T to k=p*T (p is just dimensionless variable here)
Norms integratePos(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  norms.Gamma = (0.49384-0.463155/(p*p)+0.15763/p-0.16124/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_gqq = Nf*(0.04392766/(p*p)-0.05353506/p+0.03259168/pow(p,0.8)+0.00191753/pow(p,0.2))
    *4*PI*alpha_s*4*PI*alpha_s*T; 
  norms.Gamma_ggg = (1.1104368+1.433124/(p*p*p)+0.16666408/p-0.33858792/sqrt(p))*4*PI*alpha_s*4*PI*alpha_s*T;
  return norms;
}

// approximation of the area under the transition rate from k=-12 T to k=p*T (where k<-0.05)
Norms integrateNeg(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  norms.Gamma = (0.3825932+0.03217262/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_gqq = Nf*(0.05749489/(p*p)-0.03112226/pow(p,1.8)+0.00445603/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  norms.Gamma_ggg = (0.85739544+0.2125156/p)*4*PI*alpha_s*4*PI*alpha_s*T;
  return norms;
}

Norms integrate(double p, double T, double alpha_s, int Nf)
{
  Norms norms;
  norms.Gamma = integratePos(p,T,alpha_s, Nf).Gamma + integrateNeg(p,T,alpha_s, Nf).Gamma;
  norms.Gamma_gqq = integratePos(p,T,alpha_s, Nf).Gamma_gqq + integrateNeg(p,T,alpha_s, Nf).Gamma_gqq;
  norms.Gamma_ggg = integratePos(p,T,alpha_s, Nf).Gamma_ggg + integrateNeg(p,T,alpha_s, Nf).Gamma_ggg;
  return norms;
}

// calculates the area under the envelope function when using the rejection method
// (integrals had been solved analytically before)
double area (double y, double u, int posNegSwitch, int process, Import * import)
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
}

ReturnValue findValuesRejection(double u, double T, double alpha_s, int Nf, Random *random, Import *import, int process)
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
  //f.rejections=i;
  //f.acceptances=j;
  return f;
}

void generateSpectrum(double p, double T, double alpha_s, int Nf, Random *random, Import *import, int process)
{
  const int bins = 800;
  double scale = 100;
  int binning[bins];
  ReturnValue f;
  struct timeb start, end;
  for(int iy=0; iy<bins; iy++)
    {
      binning[iy]=0;
    }
  int posy;
  int runs=200000;
  int totrej=0;
  int totacc=runs;
  
  ftime( &start );
  for(int a=0; a<runs; a++)
    {
      f = findValuesRejection(p, T, alpha_s, Nf, random, import, process);
      posy = floor((f.y+12.)*(bins/scale));
      if(posy<bins) binning[posy]+=1;
      totrej+=f.rejections;
    }
  ftime( &end );
  
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale)-12. << " " << 0.078*bins*static_cast<double>(binning[iy])/runs << endl;
    }
  double sum=0;
  for(int iy=0; iy<bins; iy++)
    {
      sum+=binning[iy];
    }
  cout << endl;
  cout << "sum=" << sum/runs << endl;
  cout << "acceptance ratio=" << static_cast<double>(totacc)/(totacc+totrej) << endl;
  cout << "time needed for " 
       << runs << " samplings: " << (end.time*1000+end.millitm-start.time*1000-start.millitm) << " ms." <<  endl;
  cout << "time per sampling: " 
       << static_cast<double>(end.time*1000+end.millitm-start.time*1000-start.millitm)/runs << " ms." <<  endl;
}


// evolve every parton by one time step.
void evolve(vector<Parton>  *plist, double T, double alpha_s, int Nf, double dtfm, 
	    Random * random, Import * import, int counter)
{
  int imax = plist->size();
  int id, newID;
  ReturnValue f;
  double dt = dtfm/hbarc;
  double p;
  Norms n;
  Vec4 vecp, newvecp;
  Parton newOne;
  double r;
  int col, acol, newCol, newAcol;

  for ( int i=0; i<imax; i++) // loop over all partons
    { counter++;
      id = plist->at(i).id();
      vecp = plist->at(i).p();
      p = vecp.pAbs();
      n = integrate(p/T,T,alpha_s, Nf);

      // cout << "particle " << i << " has ID=" << plist->at(i).id() << " p=" << plist->at(i).p();

      if ( abs(id) > 0 && id < 4 ) // if parton is a quark (u, d, or s), let it evolve like a quark
	{
	  if(random->genrand64_real1() < dt*(n.Gamma)) // see if emission happens
	    {
	      if(p/T>4.01) // do not evolve partons with momenta below this scale 
		{
		  f = findValuesRejection(p/T, T, alpha_s, Nf, random, import, 1); // do process 1, q->qg
		  p = p - f.y*T;
		  n = integrate(p/T,T,alpha_s, Nf);   // n stores the areas
		  newvecp=p*vecp/vecp.pAbs();         // new quark momentum
		  plist->at(i).p(newvecp);            // change quark's momentum to the new one
		  col = plist->at(i).col();           // color  
		  acol = plist->at(i).acol();         // anti-color
		  if (col>0)
		    {
		      plist->at(i).col(500+counter);  // set color
		      plist->at(i).acol(0);           // set anti-color
		    }
		  else if (acol>0)
		    {
		      plist->at(i).col(0);            // set color
		      plist->at(i).acol(500+counter); // set anti-color
		    }
		  plist->at(i).splits(plist->at(i).splits()+1);
		  //if (f.y>4.01) - removed to keep all color strings
		  newOne.p(f.y*T*vecp.px()/vecp.pAbs(),f.y*T*vecp.py()/vecp.pAbs(),f.y*T*vecp.pz()/vecp.pAbs()); 
		  // emitted gluon's momentum (collinear)
		  if ( col>0 )                        // if we had a quark
		    {
		      newOne.col(col);                // set new particle's color
		      newOne.acol(500+counter);       // set new particle's anti-color 
		    }
		  else if ( acol>0 )                  // if we had an anti-quark
		    {
		      newOne.col(500+counter);        // set new particle's color
		      newOne.acol(col);               // set new particle's anti-color   
		    }			  
		  // cout << "parton " << i << " radiated gluon with p=" << newOne.p();
		  newOne.id(21);                      // emitted parton is a gluon
		  newOne.splits(0);
		  plist->push_back(newOne);           // add the gluon to the list of partons
		}
	    }
	}
      else if ( id == 21 ) // if parton is a gluon, let it evolve like a gluon
	{
	  if(random->genrand64_real1() < dt*(n.Gamma_gqq)) // g -> qq
	    {
	      if(p/T>4.01) 
		{
		  f = findValuesRejection(p/T, T, alpha_s, Nf, random, import, 2); // do process 2
		  p = p - f.y*T;
		  n = integrate(p/T,T,alpha_s, Nf);
		  // choose if it is a u-ubar, d-dbar or s-sbar pair:
		  r = random->genrand64_real1();
		  if (r<0.33) newID=1;
		  else if  (r<0.66) newID=2;
		  else newID=3;
		  newvecp=p*vecp/vecp.pAbs();       // new quark momentum
		  plist->at(i).id(newID);           // turn gluon into quark 
		  // ########## now anti-quark has always smaller momentum - take care of that
		  col = plist->at(i).col();
		  acol = plist->at(i).acol();
		  plist->at(i).col(col);
		  plist->at(i).acol(0);
		  plist->at(i).splits(plist->at(i).splits()+1);
		  plist->at(i).p(newvecp);          // change quark's momentum to the new one
		  //if (f.y>4.01) - removed to keep all color strings
		  newOne.p(f.y*T*vecp.px()/vecp.pAbs(),f.y*T*vecp.py()/vecp.pAbs(),f.y*T*vecp.pz()/vecp.pAbs()); 
		  // second quark's momentum (collinear)
		  newOne.id(-newID);             // second new particle is an anti-quark
		  newOne.col(0);
		  newOne.acol(acol);
		  newOne.splits(0);
		  plist->push_back(newOne);     // add the second quark to the list of partons
		}
	    }
	  if(random->genrand64_real1() < dt*(n.Gamma_ggg)) // g -> gg
	    {
	      if(p/T>4.01) 
		{
		  f = findValuesRejection(p/T, T, alpha_s, Nf, random, import, 3); // do process 3
		  p = p - f.y*T;
		  n = integrate(p/T,T,alpha_s, Nf);
		  newvecp=p*vecp/vecp.pAbs();       // new gluon momentum
		  plist->at(i).p(newvecp);          // change gluon's momentum to the new one
		  col = plist->at(i).col();
		  acol = plist->at(i).acol();
		  plist->at(i).acol(1000+counter);
		  plist->at(i).splits(plist->at(i).splits()+1);
		  //if (f.y>4.01) - removed to keep all color strings
		  newOne.p(f.y*T*vecp.px()/vecp.pAbs(),f.y*T*vecp.py()/vecp.pAbs(),f.y*T*vecp.pz()/vecp.pAbs()); 
		  // emitted gluon's momentum (collinear)
		  newOne.id(21);                 // second new particle is a gluon too
		  newOne.col(1000+counter);
		  newOne.acol(acol);
		  newOne.splits(0);
		  plist->push_back(newOne);     // add the second gluon to the list of partons
		}
	    }
	}
    }
}

void quarkBrick(double initp, double T, double alpha_s, int Nf, double dtfm, double maxTime, Random * random, Import * import, int runs)
{
  ReturnValue f;
  struct timeb start, end;
  double p = initp;
  double dt = dtfm/hbarc;
  double mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  cout << "mt=" << mt << endl;
  Norms n = integrate(p/T,T,alpha_s, Nf);
  if (dt*n.Gamma>1.) 
    cout << "WARNING: probability to emit during one time step " << dt*n.Gamma 
	 << " > 1. Decrease time step to improve result." << endl;
  
  const int bins = 100;
  double scale = 12;
  int binning[bins];
  for(int iy=0; iy<bins; iy++)
    {
      binning[iy]=0;
    }
  
  int posy;
  ftime( &start );
  for(int j=1; j<runs; j++)
    {
      p=initp;
      n = integrate(p/T,T,alpha_s, Nf);
      for(int i=0; i<mt; i++)
	{
	  if(random->genrand64_real1() < dt*(n.Gamma))
	    {
	      if(p/T>4.01) 
		{
		  f = findValuesRejection(p/T, T, alpha_s, Nf, random, import, 1);
		  p = p - f.y*T;
		  n = integrate(p/T,T,alpha_s, Nf);
		}
	    }
	}
      posy = floor(p*(bins/scale));
      if(posy>0 && posy<bins) binning[posy]+=1;
    }
  ftime( &end );
  
  for(int iy=0; iy<bins; iy++)
    {
      cout << (iy)/(bins/scale) << " " << bins*static_cast<double>(binning[iy])/runs/scale << endl;
    }
  
  double sum=0;
  for(int iy=0; iy<bins; iy++)
    {
      sum+= bins*static_cast<double>(binning[iy])/runs/scale*(scale/bins);
    }
  
  cout << "sum=" << sum << endl;
  cout << "time needed for " 
       << runs << " events: " << (end.time*1000+end.millitm-start.time*1000-start.millitm) << " ms." <<  endl;

}

void initPythia(Pythia *pythia)
{
  pythia->readString("Random:setSeed = on");

  // set random seed as current time
  srand(time(0));
  int now=rand()/100;
  stringstream nowstr;
  nowstr << "Random:seed = " << now;
  string timestring = nowstr.str();
  pythia->readString(timestring);    

  // turn on all hard processes
  pythia->readString("HardQCD:all = on");
  
  // finish after doing hard process
  pythia->readString("PartonLevel:all = on"); // off to only get hard process partons
  pythia->readString("HadronLevel:all = off");

  // tweak pythia 8
  pythia->readString("SigmaProcess:factorMultFac = 1.");
  pythia->readString("SigmaProcess:renormMultFac = 1.");
  pythia->readString("SigmaProcess:factorScale2 = 3");
  pythia->readString("SigmaProcess:renormScale2 = 3");
  pythia->readString("SigmaProcess:factorScale3 = 4");
  pythia->readString("SigmaProcess:renormScale3 = 4");
  
  pythia->readString("StringPT:enhancedFraction = 0.05"); //default=0.01
  pythia->readString("StringPT:sigma = 0.6"); // default=0.36

  // stop shower at scale p_T given here:
  pythia->readString("SpaceShower:pTmin = 0.33"); // default = 0.2
  pythia->readString("TimeShower:pTmin = 0.5"); // default = 0.5

  //pythia->readString("ProcessLevel:all = on");

  pythia->readString("Check:event = off");

  // initialize kind of collision
  pythia->init( 2212, 2212, 200.); // 2 protons, (if 14000 each 7000 GeV)
}

void generateEventWithShower(Pythia *pythia, vector<Parton> *plist)
{
  Parton parton;
  pythia->next();
  //pythia->event.list(); 
  //cout << "particle entries=" << pythia->process.size() << endl;
  for (int i = 0; i < pythia->event.size(); ++i) 
    {
      if (pythia->event[i].isFinal())// && (pythia->event[i].id()<4 || pythia->event[i].id()==21))
	{
	  // set parton id
	  parton.id(pythia->event[i].id());
	  // set mass
	  parton.mass(pythia->event[i].m());
	  // set color 
	  parton.col(pythia->event[i].col());
	  parton.acol(pythia->event[i].acol());
	  // set momentum
	  parton.p(pythia->event[i].p());
	  parton.splits(0);
	  //cout << "partons.initParton[" << j << "].p()=" << parton[j].p() 
	  //     << " col=" << parton[j].col() << " acol=" << parton[j].acol() << endl;
	  plist->push_back(parton);
	}
    } 
}

void generateEvent(Pythia *pythia, vector<Parton> *plist)
{
  int j = 0;
  Parton parton[2];
  pythia->next();
  
  //pythia->info.list(); pythia->event.list(); 
  //pythia->process.list(); 
  
  //cout << "particle entries=" << pythia->process.size() << endl;
  for (int i = 0; i < pythia->process.size(); ++i) 
    {
      if (abs(pythia->process[i].status()) == 23) 
	{
	  // set parton id
	  if ( pythia->process[i].isGluon() ) parton[j].id(21);
	  if ( pythia->process[i].isQuark() ) parton[j].id(pythia->process[i].id());
	  // set color 
	  parton[j].col(pythia->process[i].col());
	  parton[j].acol(pythia->process[i].acol());
	  // set momentum
	  parton[j].p(pythia->process[i].p());
	  parton[j].splits(0);
	  //cout << "partons.initParton[" << j << "].p()=" << parton[j].p() 
	  //     << " col=" << parton[j].col() << " acol=" << parton[j].acol() << endl;
	  plist->push_back(parton[j]);
	  j++;
	}
    } 
}


int fragmentation(vector<Parton> * plist, Pythia *pythia, double T, Random * random, int HIC)
{
  pythia->event.reset(); // clear the event record to fill in the partons resulting from evolution below
  int j = 0;
  int id, col, acol;
  int p_imax = plist->size();
  int numEvents = 0;
  double mass;
  Vec4 pvec; 
  //cout << "after evolution:" << endl;
  
  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
    {
      id = plist->at(p_i).id();
      pvec = plist->at(p_i).p();
      col = plist->at(p_i).col();
      acol = plist->at(p_i).acol();
      mass = plist->at(p_i).mass();
      // append( id, status, col, acol, p, mass )
      pythia->event.append(id,1,col,acol,pvec,mass); // copy existing parton into event record
    }
  //pythia->event.list(); 
  //cout << "before hadronlevel" << endl;
  if(pythia->forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest
  //cout << "after hadronlevel" << endl;
  return numEvents;
  //pythia->event.list(); 
}

// test the thermal distribution sampling
void thermalPlot(double T, Random * random)
{
  double runs = 100000;
  const int bins = 50;
  double scale = 4;
  int binning[bins];
  int binning_g[bins];
  for(int iy=0; iy<bins; iy++)
    {
      binning[iy]=0;
      binning_g[iy]=0;
    }
  int posy;
  
  int p_imax;
  double p;
  
  for (int i=0; i<runs; i++)
    {
      // second argument of thermal: -1 Bose, 0 Boltzmann, +1 Fermi
      p = random->thermal(T,0).pAbs();
      posy = floor(p*(bins/scale));
      if(posy>0 && posy<bins)
	{
	  binning[posy]+=1;
	}
    }
  double sum = 0.;
  for(int iy=0; iy<bins; iy++)
    {
      cout << (iy)/(bins/scale) << " " << bins*static_cast<double>(binning[iy])/scale/runs << endl;
      sum += bins*static_cast<double>(binning[iy])/scale/runs*scale/bins;
    }
  cout << "sum=" << sum << endl;
}

int main(int argc, char* argv[]) 
{
  double dt = 0.01; // fm 
  double T = 0.3;   // GeV
  double alpha_s = 0.3;
  int fixedEnergy = 0; // switch between phythia generated initial partons (0) or fixed energy initial partons (1)
  int fragmentationSwitch = 1;
  int evolution = 1;
  int numEvents = 0;
  int Nf = 3;
  double maxTime = 2; // fm
  const int runs = 10000;
  int counter = 1;
  int countHardGluons = 0;
  int countSplits = 0;
  int countSplitsOfInitialPartons = 0;
  double avSplits = 0.;
  double p;
  int pisum = 0;
  
  int mt; // maximal time steps
  HardPartons initialPartons;
  Vec4 vecp, vecp2; // this is a pythia four-vector object
  
  // initialize new Import and Random object
  Import * import = new Import();
  Random * random = new Random();
  
  ///read in data file with transition rate
  import->init();
  
    /// initialize random number generator with current time as seed
  random->init_genrand64(time(0));
  
  /// --- do things -------------------------------------------------//
  
  // output spectrum using many samples (jet energy = first value)
    //generateSpectrum(50.9957, T, alpha_s, Nf, random, import,2); 
  
  // calculate distribution of quark going through brick
      //quarkBrick(p, T, alpha_s, Nf, dt, 2, random, import, 50000);
  
  // sample thermal distribution as test
  //thermalPlot(T,random);
  
  // initialize pythia object
  
  
  Pythia * pythia = new Pythia();
  
  initPythia(pythia);
  Parton jp1;
  Parton jp2;
  
  
  vector<Parton> ** plist;
  plist = new vector<Parton> *[runs]; 
  
  // make them both quarks with 11 GeV for testing:
  
  if ( fixedEnergy == 1)
    {
      jp1.id(1);
      jp2.id(1);
      
      jp1.p(10.,0.,0.);
      jp2.p(-10.,0.,0.);
      
      jp1.col(101);
      jp1.acol(102);
      jp2.col(103);
      jp2.acol(104);
      
      jp1.splits(0);
      jp2.splits(0);
    }
  
  // init binning
  const int bins = 20;
  double scale = 10;
  int binning[bins];
  int binning_g[bins];
  for(int iy=0; iy<bins; iy++)
    {
      binning[iy]=0;
      binning_g[iy]=0;
    }
  int posy;
  
  int p_imax;
  int id;
  
  mt = static_cast<int>(maxTime/dt+0.0001); // compute number of steps
  
  for (int j=0; j<runs; j++) // loop over all events
    {
      plist[j] = new vector<Parton>;
      
      if( fixedEnergy == 0)
	{
	  //generate pythia event and fill produced _hard_ partons in plist[j]
	  generateEventWithShower(pythia,plist[j]);
	  //generateEvent(pythia,plist[j]);
	  //pythia->event.list();
	  //	  int plsize=plist[j]->size();
	  //for (int i=0; i<plsize; i++)
	  //  cout << "list(" << i << "): " << plist[j]->at(i).id() << " with p=" << plist[j]->at(i).p();
	}
      else
	{
	  plist[j]->push_back( jp1 );
	  plist[j]->push_back( jp2 );
	}
      
      // evolve in medium if settings allow it
      if (evolution == 1)
	{
	  for(int i=0; i<mt; i++) // loop over all time steps 
	    {
	      evolve(plist[j], T, alpha_s, Nf, dt, random, import, counter);
	      counter+=1;
	    }
	}
      
      
      // fragmentation
      if (fragmentationSwitch == 1)
	{
	  numEvents+=fragmentation( plist[j], pythia, T, random, evolution );
	  double pl;
	  double En;
	  double y; // rapidity
	  double eta; //pseudo-rapidity
	  p_imax = pythia->event.size();
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	    {
	      id = pythia->event[p_i].id();
	      // count pi_0s (111) or pi+ (211)
	      if ( id == 111)
		{		
		  p = pythia->event[p_i].pT(); //p_trans
		  posy = floor(p*(bins/scale));
		  pl = pythia->event[p_i].pz(); //p_long
		  En = pythia->event[p_i].e(); //energy
		  y = 0.5*log((En+pl)/(En-pl)); //rapidity
		  eta = 0.5*log((sqrt(pl*pl+p*p)+pl)/(sqrt(pl*pl+p*p)-pl)); //pseudo-rapidity
		  //cout << "E=" << En << ", pl=" << pl << ", y=" << y << endl;
		  if(posy>=0 && posy<bins && abs(eta)<=1.)
		    {
		      binning[posy] += 1;
		      pisum += 1;
		    }
		}
	    }
	}
      else
	{
	  p_imax = plist[j]->size();
	  for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons
	    {
	      double pl;
	      double En;
	      double y; // rapidity
	      id = plist[j]->at(p_i).id();
	      p = sqrt(pow(plist[j]->at(p_i).p().px(),2.)+pow(plist[j]->at(p_i).p().py(),2.));
	      posy = floor(p*(bins/scale));
	      pl = pythia->event[p_i].pz(); //p_long
	      En = pythia->event[p_i].e(); //energy
	      y = 0.5*log((En+pl)/(En-pl)); //rapidity
	      
	      countSplits+=plist[j]->at(p_i).splits();
	      avSplits+=plist[j]->at(p_i).splits()/plist[j]->size();
	      if (p_i==0 ||p_i==1) countSplitsOfInitialPartons+=plist[j]->at(p_i).splits();
	      
	      if (p>500.) countHardGluons++;
	      if ( abs(y)<=1. )
		{
		  if( id > 0 && id < 4 ) if(posy>0 && posy<bins) binning[posy]+=1;
		  if( id == 21 ) if(posy>0 && posy<bins) binning_g[posy]+=1;
		}
	      //cout << "p. " << p_i << " in list " << j 
	      //<< " has ID=" << plist[j]->at(p_i).id() << " p=" << plist[j]->at(p_i).p();
	    }
	}
      if(j%200==0) cout << "#" << j  << endl;

    } // end loop over all events
  
  if (fragmentationSwitch==0)
    {
      for(int iy=0; iy<bins; iy++)
	{
	  cout << (iy)/(bins/scale) << " " << bins*static_cast<double>(binning[iy])/scale/2./runs //devide by runs
	       << " " << bins*static_cast<double>(binning_g[iy])/scale/2./runs << endl; //devide by runs
	  //divide by 2 since we have 2 initial particles.
	}
      
      cout << "total number of splittings per run=" << countSplits/runs 
	   << ", splittings per parton=" << avSplits/runs << endl;
      cout << "total number of splittings of the initial partons per run=" << countSplitsOfInitialPartons/runs << endl;
      cout << "hard gluons with p>500 GeV = " << countHardGluons << endl;
    }
  else
    {
      for(int iy=0; iy<bins; iy++)
	{
	  cout.precision(8);  
	  cout << ((iy)/(bins/scale)+(iy+1)/(bins/scale))/2. << " " << bins*static_cast<double>(binning[iy])/scale/numEvents/2. << endl; 
	  // divide by 2 to get 1 unit of rapidity (\Delta y = 2). (scale/bins)=\Delta p_t = bin width. 
	}
    
    }  
  double sum = 0.;
  for(int iy=0; iy<bins; iy++)
    sum += bins*static_cast<double>(binning[iy])/scale/runs*scale/bins;
  cout << "sum=" << sum << endl;
  cout << "pisum=" << pisum << endl;
  
  
  cout << "time steps=" << mt << endl;
  
  delete pythia;
  delete random;
  delete import;
  for (int j=0;j<runs;j++) delete plist[j];
  delete plist;
  
}

