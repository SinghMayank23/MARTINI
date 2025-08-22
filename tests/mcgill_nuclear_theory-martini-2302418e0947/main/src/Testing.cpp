// Testing.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains several testing routines

#include "Testing.h"

void Testing::generateSpectrum(double p, double T, double alpha_s, int Nf,
                               Random *random, Import *import, Rates * rates, int process)
{
  cout.precision(10);
  const int bins = 1000;
  double scale = 2000;
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
      f = rates->findValuesRejection(p, T, alpha_s, Nf, random, import, process);
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

void Testing::sampleElasticRate(double p, double T, double alpha_s, int Nf,
                               Random *random, Import *import, Elastic * elastic, int process)
{
  cout << " sampling elastic rate" << endl;
  const int bins = 800;
  double scale = 100;
  int binning[bins];
  double y;
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
  cout << "starting" << endl;
  for(int a=0; a<runs; a++)
    {
      y = elastic->findValuesRejection(p, T, alpha_s, Nf, random, import, process);
      posy = floor((y+12.)*(bins/scale));
      if(posy<bins) binning[posy]+=1;
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
  cout << "time needed for " 
       << runs << " samplings: " << (end.time*1000+end.millitm-start.time*1000-start.millitm) << " ms." <<  endl;
  cout << "time per sampling: " 
       << static_cast<double>(end.time*1000+end.millitm-start.time*1000-start.millitm)/runs << " ms." <<  endl;
}

void Testing::sampleElasticRateOmegaQ(double p, double omega, double T, double alpha_s, int Nf,
				      Random *random, Import *import, Elastic * elastic, int process)
{
  cout.precision(6);
  const int bins = 200;
  double scale = 400;
  int binning[bins];
  double y;
  struct timeb start, end;
  for(int iy=0; iy<bins; iy++)
    {
      binning[iy]=0;
    }
  int posy;
  int runs=5000;
  int totrej=0;
  int totacc=runs;
  
  ftime( &start );
  cout << "starting" << endl;
  for(int a=0; a<runs; a++)
    {
      if (a%5000==0) cout << "# " << a << endl;
      y = elastic->findValuesRejectionOmegaQ(p, omega, T, alpha_s, Nf, random, import, process);
      posy = floor((y)*(bins/scale));
      if(posy<bins) binning[posy]+=1;
    }
  ftime( &end );
  
  
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale) << " " << bins*static_cast<double>(binning[iy])/runs << endl;
    }
  double sum=0;
  for(int iy=0; iy<bins; iy++)
    {
      sum+=binning[iy];
    }
  cout << endl;
  cout << "sum=" << sum/runs << endl;
  cout << "time needed for " 
       << runs << " samplings: " << (end.time*1000+end.millitm-start.time*1000-start.millitm) << " ms." <<  endl;
  cout << "time per sampling: " 
       << static_cast<double>(end.time*1000+end.millitm-start.time*1000-start.millitm)/runs << " ms." <<  endl;
}

void Testing::quarkBrick(double initp, double T, double alpha_s, int Nf, double dtfm, double maxTime, Random * random, Import * import, Rates * rates, int runs)
{
  ReturnValue f;
  struct timeb start, end;
  double p = initp;
  double dt = dtfm/hbarc;
  double mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  cout << "mt=" << mt << endl;
  Norms n = rates->integrate(p/T,T,alpha_s, Nf,0,0,0);
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
      n = rates->integrate(p/T,T,alpha_s, Nf,0,0,0);
      for(int i=0; i<mt; i++)
	{
	  if(random->genrand64_real1() < dt*(n.Gamma))
	    {
	      if(p/T>4.01) 
		{
		  f = rates->findValuesRejection(p/T, T, alpha_s, Nf, random, import, 1);
		  p = p - f.y*T;
		  n = rates->integrate(p/T,T,alpha_s, Nf,0,0,0);
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

// test the thermal distribution sampling
void Testing::thermalPlot(double T, Random * random)
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

void Testing::initialDist(Random *random, Glauber *glauber)
{
  const int bins = 30;
  double scale = 10.;
  int binning[bins][bins];
  ReturnValue f;
  struct timeb start, end;
  for(int ix=0; ix<bins; ix++)
    for(int iy=0; iy<bins; iy++)
      {
	binning[ix][iy]=0.;
      }
  int posx, posy;
  int runs=200000;
  int totrej=0;
  int totacc=runs;
  
  ftime( &start );
  for(int a=0; a<runs; a++)
    {
      if(a%1000==0) cout << a << endl;
      f = glauber->SamplePABRejection(random);
      posx = floor((f.x)*(bins/scale));
      posy = floor((f.y)*(bins/scale));
      if(posx<bins && posy<bins && posx>=0 && posy>=0)
	{
	  binning[posx][posy]+=1;
	}
    }
  ftime( &end );
  
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale) << " " << bins*static_cast<double>(binning[0][iy])/runs << endl;
    }
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale) << " " << bins*static_cast<double>(binning[iy][0])/runs << endl;
    }
  exit(1);
}

void Testing::sampleTA(Random *random, Glauber *glauber)
{
  const int bins = 30;
  double scale = 15.;
  int binning[bins][bins];
  ReturnValue f;
  struct timeb start, end;
  for(int ix=0; ix<bins; ix++)
    for(int iy=0; iy<bins; iy++)
      {
	binning[ix][iy]=0.;
      }
  int posx, posy;
  int runs=200000;
  int totrej=0;
  int totacc=runs;
  
  ftime( &start );
  for(int a=0; a<runs; a++)
    {
      if(a%1000==0) cout << a << endl;
      f = glauber->SampleTARejection(random);
      posx = floor((f.x)*(bins/scale));
      posy = floor((f.y)*(bins/scale));
      //cout << "x=" << f.x << " y=" << f.y << endl;
      if(posx<bins && posy<bins && posx>=0 && posy>=0)
	{
	  binning[posx][posy]+=1;
	}
    }
  ftime( &end );
  
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale) << " " << bins*static_cast<double>(binning[0][iy])/runs << endl;
    }
  for(int iy=0; iy<bins; iy++)
    {
      cout << ((iy))/(bins/scale) << " " << bins*static_cast<double>(binning[iy][0])/runs << endl;
    }
  exit(1);
}

