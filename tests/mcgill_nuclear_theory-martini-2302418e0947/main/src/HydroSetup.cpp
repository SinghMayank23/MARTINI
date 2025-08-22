// HydroSetup.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions that return interpolated data at a given space-time point

#include "HydroSetup.h"


// all hydro data is stored in tau steps (not t) - the t and z in the MARTINI evolution is converted to tau when accessing the hydro data
void HydroSetup::readHydroData(double A, double tau0, double taumax, double dtau, 
			       double xmax, double zmax, double dx, double dz, 
			       int whichHydro, int subset, bool viscous, double b, vector<HydroCell> *lattice, string evolution_name)
{
  //Clear the lattice, in case it has already been instantiated:
  lattice->clear();

  if (whichHydro==1)
    {
      cout << "Using 2+1D Kolb Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[4];
      string path     = "";
      string subfolder;
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if ( A == 197 ) //Au
	{
	  if ( b == 0 ) subfolder = "/hydro/kolb/b0/";
	  else if ( b == 2.4 ) subfolder = "/hydro/kolb/b2.4/";
	  else if ( b == 7.5 ) subfolder = "/hydro/kolb/b7.5/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=0, 2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      else if ( A == 63 ) //Cu
	{
	  if ( b == 0 ) subfolder = "/hydro/kolb/b0CuCu/";
	  else if ( b == 1.7 ) subfolder = "/hydro/kolb/b1.7CuCu/";
	  else if ( b == 2.4 ) subfolder = "/hydro/kolb/b2.4CuCu/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=0, 1.7, 2.4 fm. " << endl;
	      exit(1);
	    }
	}
      else
	{
	  cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given nucleus-nucleus collision. A=" << A << endl;
	  cout << " available nuclei are are A=197 (Au) and A=63 (Cu)." << endl;
	  exit(1);
	}

      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += subfolder;
	}
      else path = ".."+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double T, vx, vy, tau, epsilon, QGPfrac;
      HydroCell newCell;  
      int position=0;
      // position = ix+ixmax*(iy+iymax*itau)
      
      // make sure you that MARTINI won't try files that aren't there:
      if (b==7.5 && itaumax>302)
	itaumax=302;
      if (b==2.4 && itaumax>398)
	itaumax=398;
      if (b==0. && itaumax>415)
	itaumax=415;

      // open files with hydro data to read in:
      for(int i=0; i<itaumax; i++)
	{
	  tau=tau0+i*dtau;
	  temp.str("");
	  temp << i+1;
	  istr = temp.str();
	  if (i%5==0) cout << "." << flush;
	  if (i+1<10) 
	    {
	      file[0] = path + "T00" + istr + ".dat";
	      file[1] = path + "VX00" + istr + ".dat";
	      file[2] = path + "VY00" + istr + ".dat";
	      file[3] = path + "EPS00" + istr + ".dat";
	    }
	  else if (i+1<100)
	    {
	      file[0] = path + "T0" + istr + ".dat";
	      file[1] = path + "VX0" + istr + ".dat";
	      file[2] = path + "VY0" + istr + ".dat";
	      file[3] = path + "EPS0" + istr + ".dat";
	    }
	  else 
	    {
	      file[0] = path + "T" + istr + ".dat";
	      file[1] = path + "VX" + istr + ".dat";
	      file[2] = path + "VY" + istr + ".dat";
	      file[3] = path + "EPS" + istr + ".dat";
	    }
	  // temperature
	  fin.open(file[0].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	      exit(1);
	    }
	  int ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> T;
	      newCell.T(T);
	      if (ik<=ixmax*ixmax) lattice->push_back(newCell);
	    }
	  fin.close();
	  // flow in x direction
	  fin.open(file[1].c_str(),ios::in);
	 
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> vx;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).vx(vx);
		  position++;
		}
	    }
	  fin.close();
	  // flow in y direction
	  fin.open(file[2].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[2] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> vy;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).vy(vy);
		  position++;
		}
	    }
	  fin.close();
	  // energy density / QGP fraction
	  fin.open(file[3].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[3] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> epsilon;
	      QGPfrac=(epsilon-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
	      if (QGPfrac>1.) QGPfrac = 1;
	      else if (QGPfrac<0.) QGPfrac=0.;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).QGPfrac(QGPfrac);
		  position++;
		}
	    }
	  fin.close();
	}
      cout << " ok." << endl;
    }
  else if (whichHydro == 2) // read in Eskola hydro
    {
      cout << "Using 2+1D Eskola Hydro: reading data ..." << endl;
      // r changes first (inner loop), tau (outer loop)
      int static const itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int static const ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[3];
      double r, x ,y;
      double TofRandTau[itaumax][300];
      double VrofRandTau[itaumax][300];
      string path     = "";
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/RHIC200/";
	}
      else path = "../hydro/RHIC200/";
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double T, vx, vy, vr;
      HydroCell newCell;  
      int position=0;
      // position = ix+ixmax*(iy+iymax*itau)
      
      // open files with hydro data to read in:
      file[0] = path + "TMATin";
      file[1] = path + "VrMATin";
      // temperature
      fin.open(file[0].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      int i=0;
      int ir, itau;
      cout.precision(6);
      while ( !fin.eof() )
	{
	  fin >> T;
	  ir=(i)%ixmax;
	  itau=floor((i)/ixmax);
	  //if (i%1000==0) cout << "i=" << i << " itau=" << itau << " ir=" << ir << " T=" << T << endl;
	  if (itau<itaumax && ir<ixmax) TofRandTau[itau][ir] = T;
	  i++;
	}
      fin.close();
      
      for(itau=0; itau<itaumax; itau++)
	for(int ix=0; ix<ixmax; ix++)
	  for(int iy=0; iy<ixmax; iy++)
	    {
	      x = -xmax + ix*dx;
	      y = -xmax + iy*dx;
	      r=sqrt(x*x+y*y);
	      ir=floor(r/dx);
	      if ( ir<ixmax ) T=TofRandTau[itau][ir];
	      else T=0.;
	      newCell.T(T);
	      lattice->push_back(newCell);
	      //cout << "ix=" << ix << " iy=" << iy << " ir=" << ir << " r=" << r << " T=" << T << endl; 
	    }
      // flow in r direction
      fin.open(file[1].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	  exit(1);
	}
      cout.precision(6);
      i=0;
      while ( !fin.eof() )
	{
	  fin >> vr;
	  ir=(i)%ixmax;
	  itau=floor((i)/ixmax);
	  //if (i%1000==0) cout << "i=" << i << " itau=" << itau << " ir=" << ir << " Vr=" << vr << endl;
	  if (itau<itaumax && ir<ixmax) VrofRandTau[itau][ir] = vr;
	  i++;
	}
      fin.close();
      
      for(itau=0; itau<itaumax; itau++)
	for(int ix=0; ix<ixmax; ix++)
	  for(int iy=0; iy<ixmax; iy++)
	    {
	      x = -xmax + ix*dx;
	      y = -xmax + iy*dx;
	      r=sqrt(x*x+y*y);
	      ir=floor(r/dx);
	      if ( ir<ixmax && r>0 )
		{
		  vx=VrofRandTau[itau][ir]*(x/r);
		  vy=VrofRandTau[itau][ir]*(y/r);
		}
	      else
		{
		  vx=0.;
		  vy=0.;
		}
	      position = ix+ixmax*(iy+ixmax*itau);
	      if(position < lattice->size()) lattice->at(position).vx(vx);
	      if(position < lattice->size()) lattice->at(position).vy(vy);
	      
	      //if (itau==299) cout << "itau=" << itau << " ix=" << ix << " iy=" << iy << " ir=" 
	      //	      << ir << " r=" << r << " vr=" << VrofRandTau[itau][ir] 
	      //	      << " myvr=" << sqrt(vx*vx+vy*vy) << " vx=" 
	      //	      << lattice->at(position).vx() << " vy=" <<  lattice->at(position).vy() << endl; 
	    }
      cout << " ok." << endl;
    }
  else if (whichHydro==3) // optimized 3+1D Nonaka and Bass Hydro
    {
      ifstream fin;
      string sent;
      ostringstream temp;
      string istr;
      double T, vx, vy, veta, QGPfrac;
      int position=0;
      HydroCell newCell;  
      string subfolder;
 
      cout << "Using 3+1D Nonaka Bass Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      //cout << "itaumax=" << itaumax << endl;
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      //cout << "ixmax=" << ixmax << endl;
      int izmax = static_cast<int>(2*zmax/dz+0.001);
      //cout << "izmax=" << izmax << endl;
      
      string file[2];
      // check if path to MARTINI folder is set and if so, use it
      string path = "";
      if ( b == 2.4 ) subfolder = "b2.4/";
      else if ( b == 4.5 ) subfolder = "b4.5/";
      else if ( b == 6.3 ) subfolder = "b6.3/";
      else if ( b == 7.5 ) subfolder = "b7.5/";
      else 
	{
	  cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	  cout << " available impact parameters are b= 2.4, 4.5, 6.3, and 7.5 fm." << endl;
	  exit(1);
	}
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/nonaka/"+subfolder;
	}
      // if environment does not know about MARTINI, use relative folder position
      // else path = "../hydro/nonaka/"+subfolder;
      
      //path = "./output/nonaka_quarter/"+subfolder;
      cout << "Path to hydro data: " << path << endl;

      // position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
      
      // open files with hydro data to read in:
    
      for(int i=0; i<itaumax; i++)
	{
	  // define file names
	  temp.str("");
	  temp << i;
	  istr = temp.str();
	  file[0] = path + "T" + istr + ".dat";
	  file[1] = path + "V" + istr + ".dat";
	  
	  // read in temperature and QGP fraction
	  fin.open(file[0].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	      exit(1);
	    }
	  if (i%5==0) cout << "." << flush; // show progress ...
	  int ik = 0;
	  while ( !fin.eof() ) // like while() which caused problems 
	    {
	      ik++;
	      fin >> T;
	      fin >> QGPfrac;
	      newCell.T(T);
	      newCell.QGPfrac(QGPfrac);
	      //cout << T << " " << QGPfrac << endl;
	      if (ik<=ixmax/2*ixmax/2*izmax) lattice->push_back(newCell);
	      //if (ik%1000==0) cout << i << " " << ik << endl;
	    }
	  fin.close();
	  // read in flow in x, y, eta direction
	  fin.open(file[1].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	      exit(1);
	    }
	  ik = 0;
	  position=ixmax/2*ixmax/2*izmax*i;
	  while ( ik<=ixmax/2*ixmax/2*izmax )
	    {
	      ik++;
	      fin >> vx;
	      fin >> vy;
	      fin >> veta;
	      //	      cout << scientific << vx << " " << vy << endl;
	      //if (ik%1000==0 && ik>ixmax/2*ixmax/2*izmax) cout << "eof" << i << " " << ik << endl;
	      if (ik<=ixmax/2*ixmax/2*izmax)
		{
		  if (position<lattice->size())
		    {
		      lattice->at(position).vx(vx);
		      lattice->at(position).vy(vy);
		      lattice->at(position).vz(veta);
		    }
		  else 
		    {
		      cout << "MARTINI:HydroSetup ERROR while reading file." << endl;
		      cout << position << endl;
		    }
		  position++;
		}
	    }
	  fin.close();
	}
      cout << " ok." << endl;
      //testing:
      /*
      int itau, ix, iy, iz;
      ix=ixmax/2;
      cout.precision(12);
      cout <<"xmax=" << xmax << " ixmax=" << ixmax << endl;
      for(itau=3;itau<4;itau++)
	for(iz=izmax/2;iz<izmax/2+1;iz++)
	  for(iy=0;iy<ixmax;iy++)
	    {
	      position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
	      cout << "tau=" << itau*dtau+tau0 << " x=" << -xmax+ix*dx 
		   << " y=" << -xmax+iy*dx << " z=" << -zmax+iz*dz
		   << " T=" << lattice->at(position).T() << endl;
	    }
      exit(1);
      */
    }
  else if (whichHydro==4 || whichHydro==6) // 3+1D MUSIC hydro (Schenke, Jeon, Gale) (for 4, just 1 quarter of the transverse plane, 6 all)
    {
      ifstream fin;
      string sent;
      ostringstream temp;
      string istr, bstr;
      double T, vx, vy, vz, QGPfrac;
      int position=0;
      HydroCell newCell;  
      string subfolder;
      cout << "Using 3+1D Jeon Schenke Hydro subset " << subset << " at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      cout << "itaumax=" << itaumax << endl;
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      //cout << "ixmax=" << ixmax << endl;
      int izmax = static_cast<int>(2*zmax/dz+0.001);
      //cout << "izmax=" << izmax << endl;
      
      string file[2];
      temp.str("");
      temp << subset;
      istr = temp.str();
      temp.str("");
      temp << b;
      bstr = temp.str();
     
      // check if path to MARTINI folder is set and if so, use it
      string path = "";
      subfolder = "b";
      subfolder += bstr;
      subfolder += "/";
      
      if ( b == 2.4 ) subfolder = "b2.4/";
      else if ( b == 1.76 ) subfolder = "b1.76/";
      else if ( b == 2.28 ) subfolder = "b2.28/";
      else if ( b == 4.16 ) subfolder = "b4.16/";
      else if ( b == 4.5 ) subfolder = "b4.5/";
      else if ( b == 5.4 ) subfolder = "b5.4/";
      else if ( b == 6.3 ) subfolder = "b6.3/";
      else if ( b == 7.5 ) subfolder = "b7.5/";
      else 
	{
	  cout << "[HydroSetup::readHydroData]: WARNING: Are you sure that hydro data for given impact parameter " << b << 
	    " exists?" << endl;
	  cout << " available impact parameters are b= 1.76 (Cu), 2.28(Cu), 2.4, 4.16(Cu), 4.5, 5.4(Cu), 6.3, and 7.5 fm." << endl;
	}
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/schenke/"+subfolder;
	}
      // if environment does not know about MARTINI, use relative folder position
      else path = "../hydro/schenke/"+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      // position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
      
      
      // open files with hydro data to read in:
      if (subset==1)
	file[0] = path + "evolution.dat";
      else
	file[0] = path + "evolution" + istr + ".dat";

      // read in temperature, QGP fraction , flow velocity
      //fin.open(file[0].c_str(),ios::in);
      //By hand, I am simply changing the input file to evolution.dat:
      //The name of the evolution file:
      string evolution_file_name = evolution_name;
      cout << "Evolution file name = " << evolution_file_name << endl;
      fin.open(evolution_file_name.c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      int ik = 0;
      while ( !fin.eof() )
	{
	  ik++;
	  fin >> T;
	  fin >> QGPfrac;
	  fin >> vx;
	  fin >> vy;
	  fin >> vz;
	  newCell.T(T);
	  newCell.QGPfrac(QGPfrac);
	  newCell.vx(vx);
	  newCell.vy(vy);
	  newCell.vz(vz);
	  //	  	  cout << T << " " << QGPfrac << " " << vx << " " << vy << " " << vz << endl;
	  lattice->push_back(newCell);
	  if (ik%50000==0) cout << "o" << flush;
	}
      cout << ik << endl;
      fin.close();
      //testing:
      
//        int itau, ix, iy, iz;
//        iy=0;
//        iz=0;
//        cout.precision(12);
//        cout <<"xmax=" << xmax << " ixmax=" << ixmax << endl;
//        for(itau=0;itau<1;itau++)
// 	 for(ix=0;ix<ixmax/2;ix++)
// 	   {
// 	     position = ix + ixmax/2*(iy+ixmax/2*(iz+izmax*itau));
// 	     cout << "tau=" << itau*dtau+tau0 << " x=" << ix*dx 
// 		  << " y=" << iy*dx << " z=" << -zmax+iz*dz
// 		  << " T=" << lattice->at(position).T() << endl;
// 	   }
//        exit(1);
       
    }
  else if (whichHydro==5)
    {
      cout << "Using 2+1D Moreland-Song-Heinz Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[4];
      string path     = "";
      string subfolder;
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (viscous)
	{
	  if ( b == 2.4 ) subfolder = "/hydro/moreland/fKLNb2.4etas0.08/";
	  else if ( b == 7.5 ) subfolder = "/hydro/moreland/fKLNb7.5etas0.08/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      else
	{
	  if ( b == 2.4 ) subfolder = "/hydro/moreland/fKLNb2.4ideal/";
	  else if ( b == 7.5 ) subfolder = "/hydro/moreland/fKLNb7.5ideal/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += subfolder;
	}
      else path = ".."+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double x, y, T, vx, vy, tau, epsilon, QGPfrac;
      HydroCell newCell;  
      int position=0;
      // position = iy+iymax*(ix+ixmax*itau) y is the INNER LOOP here
      
      // open file with hydro data to read in:
	  
      file[0] = path + "HydroEvol.dat";
      
      fin.open(file[0].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      
      int ik = 0;
      while ( !fin.eof() && ik<=ixmax*ixmax*itaumax )
	{
	  ik++;
	  fin >> x;
	  fin >> y;
	  fin >> tau;
	  fin >> epsilon;
	  fin >> T;
	  fin >> vx;
	  fin >> vy;
	  
	  T*=hbarc;
	  epsilon*=hbarc;
	  
	  QGPfrac=(epsilon-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
	  //cout << T << " " << QGPfrac << endl;
	  if (QGPfrac>1.) QGPfrac = 1;
	  else if (QGPfrac<0.) QGPfrac=0.;
	  
	  //cout << tau << " " << itaumax << " " << ik << " " << ixmax*ixmax*itaumax << endl;

	  newCell.QGPfrac(QGPfrac);
	  newCell.T(T);
	  newCell.T(T);
	  newCell.vx(vx);
	  newCell.vy(vy);
	  if (ik<=ixmax*ixmax*itaumax && (x<xmax && y<xmax)) lattice->push_back(newCell);
	}
      fin.close();
      // flow in x direction
      cout << " ok." << endl;
    }
}

void HydroSetup::readHydroData(double A, double tau0, double taumax, double dtau, 
			       double xmax, double zmax, double dx, double dz, 
			       int whichHydro, int file_number, int subset, bool viscous, double b, vector<HydroCell> *lattice, string evolution_name)
{
  if (whichHydro==1)
    {
      cout << "Using 2+1D Kolb Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[4];
      string path     = "";
      string subfolder;
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if ( A == 197 ) //Au
	{
	  if ( b == 0 ) subfolder = "/hydro/kolb/b0/";
	  else if ( b == 2.4 ) subfolder = "/hydro/kolb/b2.4/";
	  else if ( b == 7.5 ) subfolder = "/hydro/kolb/b7.5/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=0, 2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      else if ( A == 63 ) //Cu
	{
	  if ( b == 0 ) subfolder = "/hydro/kolb/b0CuCu/";
	  else if ( b == 1.7 ) subfolder = "/hydro/kolb/b1.7CuCu/";
	  else if ( b == 2.4 ) subfolder = "/hydro/kolb/b2.4CuCu/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=0, 1.7, 2.4 fm. " << endl;
	      exit(1);
	    }
	}
      else
	{
	  cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given nucleus-nucleus collision. A=" << A << endl;
	  cout << " available nuclei are are A=197 (Au) and A=63 (Cu)." << endl;
	  exit(1);
	}

      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += subfolder;
	}
      else path = ".."+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double T, vx, vy, tau, epsilon, QGPfrac;
      HydroCell newCell;  
      int position=0;
      // position = ix+ixmax*(iy+iymax*itau)
      
      // make sure you that MARTINI won't try files that aren't there:
      if (b==7.5 && itaumax>302)
	itaumax=302;
      if (b==2.4 && itaumax>398)
	itaumax=398;
      if (b==0. && itaumax>415)
	itaumax=415;

      // open files with hydro data to read in:
      for(int i=0; i<itaumax; i++)
	{
	  tau=tau0+i*dtau;
	  temp.str("");
	  temp << i+1;
	  istr = temp.str();
	  if (i%5==0) cout << "." << flush;
	  if (i+1<10) 
	    {
	      file[0] = path + "T00" + istr + ".dat";
	      file[1] = path + "VX00" + istr + ".dat";
	      file[2] = path + "VY00" + istr + ".dat";
	      file[3] = path + "EPS00" + istr + ".dat";
	    }
	  else if (i+1<100)
	    {
	      file[0] = path + "T0" + istr + ".dat";
	      file[1] = path + "VX0" + istr + ".dat";
	      file[2] = path + "VY0" + istr + ".dat";
	      file[3] = path + "EPS0" + istr + ".dat";
	    }
	  else 
	    {
	      file[0] = path + "T" + istr + ".dat";
	      file[1] = path + "VX" + istr + ".dat";
	      file[2] = path + "VY" + istr + ".dat";
	      file[3] = path + "EPS" + istr + ".dat";
	    }
	  // temperature
	  fin.open(file[0].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	      exit(1);
	    }
	  int ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> T;
	      newCell.T(T);
	      if (ik<=ixmax*ixmax) lattice->push_back(newCell);
	    }
	  fin.close();
	  // flow in x direction
	  fin.open(file[1].c_str(),ios::in);
	   
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> vx;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).vx(vx);
		  position++;
		}
	    }
	  fin.close();
	  // flow in y direction
	  fin.open(file[2].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[2] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> vy;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).vy(vy);
		  position++;
		}
	    }
	  fin.close();
	  // energy density / QGP fraction
	  fin.open(file[3].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[3] << endl;
	      exit(1);
	    }
	  position=ixmax*ixmax*i;
	  ik = 0;
	  while ( !fin.eof() )
	    {
	      ik++;
	      fin >> epsilon;
	      QGPfrac=(epsilon-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
	      if (QGPfrac>1.) QGPfrac = 1;
	      else if (QGPfrac<0.) QGPfrac=0.;
	      if (ik<=ixmax*ixmax) 
		{
		  lattice->at(position).QGPfrac(QGPfrac);
		  position++;
		}
	    }
	  fin.close();
	}
      cout << " ok." << endl;
    }
  else if (whichHydro == 2) // read in Eskola hydro
    {
      cout << "Using 2+1D Eskola Hydro: reading data ..." << endl;
      // r changes first (inner loop), tau (outer loop)
      int static const itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int static const ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[3];
      double r, x ,y;
      double TofRandTau[itaumax][300];
      double VrofRandTau[itaumax][300];
      string path     = "";
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/RHIC200/";
	}
      else path = "../hydro/RHIC200/";
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double T, vx, vy, vr;
      HydroCell newCell;  
      int position=0;
      // position = ix+ixmax*(iy+iymax*itau)
      
      // open files with hydro data to read in:
      file[0] = path + "TMATin";
      file[1] = path + "VrMATin";
      // temperature
      fin.open(file[0].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      int i=0;
      int ir, itau;
      cout.precision(6);
      while ( !fin.eof() )
	{
	  fin >> T;
	  ir=(i)%ixmax;
	  itau=floor((i)/ixmax);
	  //if (i%1000==0) cout << "i=" << i << " itau=" << itau << " ir=" << ir << " T=" << T << endl;
	  if (itau<itaumax && ir<ixmax) TofRandTau[itau][ir] = T;
	  i++;
	}
      fin.close();
      
      for(itau=0; itau<itaumax; itau++)
	for(int ix=0; ix<ixmax; ix++)
	  for(int iy=0; iy<ixmax; iy++)
	    {
	      x = -xmax + ix*dx;
	      y = -xmax + iy*dx;
	      r=sqrt(x*x+y*y);
	      ir=floor(r/dx);
	      if ( ir<ixmax ) T=TofRandTau[itau][ir];
	      else T=0.;
	      newCell.T(T);
	      lattice->push_back(newCell);
	      //cout << "ix=" << ix << " iy=" << iy << " ir=" << ir << " r=" << r << " T=" << T << endl; 
	    }
      // flow in r direction
      fin.open(file[1].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	  exit(1);
	}
      cout.precision(6);
      i=0;
      while ( !fin.eof() )
	{
	  fin >> vr;
	  ir=(i)%ixmax;
	  itau=floor((i)/ixmax);
	  //if (i%1000==0) cout << "i=" << i << " itau=" << itau << " ir=" << ir << " Vr=" << vr << endl;
	  if (itau<itaumax && ir<ixmax) VrofRandTau[itau][ir] = vr;
	  i++;
	}
      fin.close();
      
      for(itau=0; itau<itaumax; itau++)
	for(int ix=0; ix<ixmax; ix++)
	  for(int iy=0; iy<ixmax; iy++)
	    {
	      x = -xmax + ix*dx;
	      y = -xmax + iy*dx;
	      r=sqrt(x*x+y*y);
	      ir=floor(r/dx);
	      if ( ir<ixmax && r>0 )
		{
		  vx=VrofRandTau[itau][ir]*(x/r);
		  vy=VrofRandTau[itau][ir]*(y/r);
		}
	      else
		{
		  vx=0.;
		  vy=0.;
		}
	      position = ix+ixmax*(iy+ixmax*itau);
	      if(position < lattice->size()) lattice->at(position).vx(vx);
	      if(position < lattice->size()) lattice->at(position).vy(vy);
	            
	      //if (itau==299) cout << "itau=" << itau << " ix=" << ix << " iy=" << iy << " ir=" 
	      //      << ir << " r=" << r << " vr=" << VrofRandTau[itau][ir] 
	      //      << " myvr=" << sqrt(vx*vx+vy*vy) << " vx=" 
	      //      << lattice->at(position).vx() << " vy=" <<  lattice->at(position).vy() << endl; 
	    }
      cout << " ok." << endl;
    }
  else if (whichHydro==3) // optimized 3+1D Nonaka and Bass Hydro
    {
      ifstream fin;
      string sent;
      ostringstream temp;
      string istr;
      double T, vx, vy, veta, QGPfrac;
      int position=0;
      HydroCell newCell;  
      string subfolder;
 
      cout << "Using 3+1D Nonaka Bass Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      //cout << "itaumax=" << itaumax << endl;
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      //cout << "ixmax=" << ixmax << endl;
      int izmax = static_cast<int>(2*zmax/dz+0.001);
      //cout << "izmax=" << izmax << endl;
      
      string file[2];
      // check if path to MARTINI folder is set and if so, use it
      string path = "";
      if ( b == 2.4 ) subfolder = "b2.4/";
      else if ( b == 4.5 ) subfolder = "b4.5/";
      else if ( b == 6.3 ) subfolder = "b6.3/";
      else if ( b == 7.5 ) subfolder = "b7.5/";
      else 
	{
	  cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	  cout << " available impact parameters are b= 2.4, 4.5, 6.3, and 7.5 fm." << endl;
	  exit(1);
	}
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/nonaka/"+subfolder;
	}
      // if environment does not know about MARTINI, use relative folder position
      // else path = "../hydro/nonaka/"+subfolder;
      
      //path = "./output/nonaka_quarter/"+subfolder;
      cout << "Path to hydro data: " << path << endl;

      // position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
      
      // open files with hydro data to read in:
    
      for(int i=0; i<itaumax; i++)
	{
	  // define file names
	  temp.str("");
	  temp << i;
	  istr = temp.str();
	  file[0] = path + "T" + istr + ".dat";
	  file[1] = path + "V" + istr + ".dat";
	    
	  // read in temperature and QGP fraction
	  fin.open(file[0].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	      exit(1);
	    }
	  if (i%5==0) cout << "." << flush; // show progress ...
	  int ik = 0;
	  while ( !fin.eof() ) // like while() which caused problems 
	    {
	      ik++;
	      fin >> T;
	      fin >> QGPfrac;
	      newCell.T(T);
	      newCell.QGPfrac(QGPfrac);
	      //cout << T << " " << QGPfrac << endl;
	      if (ik<=ixmax/2*ixmax/2*izmax) lattice->push_back(newCell);
	      //if (ik%1000==0) cout << i << " " << ik << endl;
	    }
	  fin.close();
	  // read in flow in x, y, eta direction
	  fin.open(file[1].c_str(),ios::in);
	  if(!fin)
	    {
	      cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[1] << endl;
	      exit(1);
	    }
	  ik = 0;
	  position=ixmax/2*ixmax/2*izmax*i;
	  while ( ik<=ixmax/2*ixmax/2*izmax )
	    {
	      ik++;
	      fin >> vx;
	      fin >> vy;
	      fin >> veta;
	      //      cout << scientific << vx << " " << vy << endl;
	      //if (ik%1000==0 && ik>ixmax/2*ixmax/2*izmax) cout << "eof" << i << " " << ik << endl;
	      if (ik<=ixmax/2*ixmax/2*izmax)
		{
		  if (position<lattice->size())
		    {
		      lattice->at(position).vx(vx);
		      lattice->at(position).vy(vy);
		      lattice->at(position).vz(veta);
		    }
		  else 
		    {
		      cout << "MARTINI:HydroSetup ERROR while reading file." << endl;
		      cout << position << endl;
		    }
		  position++;
		}
	    }
	  fin.close();
	}
      cout << " ok." << endl;
      //testing:
      /*
      int itau, ix, iy, iz;
      ix=ixmax/2;
      cout.precision(12);
      cout <<"xmax=" << xmax << " ixmax=" << ixmax << endl;
      for(itau=3;itau<4;itau++)
      for(iz=izmax/2;iz<izmax/2+1;iz++)
        for(iy=0;iy<ixmax;iy++)
	    {
	          position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
		        cout << "tau=" << itau*dtau+tau0 << " x=" << -xmax+ix*dx 
			   << " y=" << -xmax+iy*dx << " z=" << -zmax+iz*dz
			      << " T=" << lattice->at(position).T() << endl;
			          }
      exit(1);
      */
    }
  else if (whichHydro==4 || whichHydro==6) // 3+1D MUSIC hydro (Schenke, Jeon, Gale) (for 4, just 1 quarter of the transverse plane, 6 all)
    {
      ifstream fin;
      string sent;
      ostringstream temp;
      string istr, bstr;
      double T, vx, vy, vz, QGPfrac;
      int position=0;
      HydroCell newCell;  
      string subfolder;
      cout << "Using 3+1D Jeon Schenke Hydro subset " << subset << " at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      cout << "itaumax=" << itaumax << endl;
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      //cout << "ixmax=" << ixmax << endl;
      int izmax = static_cast<int>(2*zmax/dz+0.001);
      //cout << "izmax=" << izmax << endl;
      
      string file[2];
      temp.str("");
      temp << subset;
      istr = temp.str();
      temp.str("");
      temp << b;
      bstr = temp.str();
     
      // check if path to MARTINI folder is set and if so, use it
      string path = "";
      subfolder = "b";
      subfolder += bstr;
      subfolder += "/";
      
      if ( b == 2.4 ) subfolder = "b2.4/";
      else if ( b == 1.76 ) subfolder = "b1.76/";
      else if ( b == 2.28 ) subfolder = "b2.28/";
      else if ( b == 4.16 ) subfolder = "b4.16/";
      else if ( b == 4.5 ) subfolder = "b4.5/";
      else if ( b == 5.4 ) subfolder = "b5.4/";
      else if ( b == 6.3 ) subfolder = "b6.3/";
      else if ( b == 7.5 ) subfolder = "b7.5/";
      else 
	{
	  cout << "[HydroSetup::readHydroData]: WARNING: Are you sure that hydro data for given impact parameter " << b << 
	    " exists?" << endl;
	  cout << " available impact parameters are b= 1.76 (Cu), 2.28(Cu), 2.4, 4.16(Cu), 4.5, 5.4(Cu), 6.3, and 7.5 fm." << endl;
	}
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += "/hydro/schenke/"+subfolder;
	}
      // if environment does not know about MARTINI, use relative folder position
      else path = "../hydro/schenke/"+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      // position = ix + ixmax*(iy+ixmax*(iz+izmax*itau));
      
      
      // open files with hydro data to read in:
      if (subset==1)
	file[0] = path + "evolution.dat";
      else
	file[0] = path + "evolution" + istr + ".dat";

      // read in temperature, QGP fraction , flow velocity
      // For now, I'm changing this by hand. -CFY
      //fin.open(file[0].c_str(),ios::in);
      stringstream s;
      s << file_number;
      string file_name =evolution_name +"_"+s.str()+".dat";
      fin.open(file_name.c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      int ik = 0;
      while ( !fin.eof() )
	{
	  ik++;
	  fin >> T;
	  fin >> QGPfrac;
	  fin >> vx;
	  fin >> vy;
	  fin >> vz;
	  newCell.T(T);
	  newCell.QGPfrac(QGPfrac);
	  newCell.vx(vx);
	  newCell.vy(vy);
	  newCell.vz(vz);
	  //    cout << T << " " << QGPfrac << " " << vx << " " << vy << " " << vz << endl;
	  lattice->push_back(newCell);
	  if (ik%50000==0) cout << "o" << flush;
	  //cout << ik << endl;
	}
      fin.close();
      //testing:
      
      //        int itau, ix, iy, iz;
      //        iy=0;
      //        iz=0;
      //        cout.precision(12);
      //        cout <<"xmax=" << xmax << " ixmax=" << ixmax << endl;
      //        for(itau=0;itau<1;itau++)
      //  for(ix=0;ix<ixmax/2;ix++)
      //    {
      //      position = ix + ixmax/2*(iy+ixmax/2*(iz+izmax*itau));
      //      cout << "tau=" << itau*dtau+tau0 << " x=" << ix*dx 
      //   << " y=" << iy*dx << " z=" << -zmax+iz*dz
      //   << " T=" << lattice->at(position).T() << endl;
      //    }
      //        exit(1);
       
    }
  else if (whichHydro==5)
    {
      cout << "Using 2+1D Moreland-Song-Heinz Hydro at impact parameter b=" << b << " fm. reading data ..." << endl;
      int itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
      int ixmax = static_cast<int>(2*xmax/dx+0.001);
      string file[4];
      string path     = "";
      string subfolder;
      const char* HYDROPATH = "HYDROPATH";
      char* envPath = getenv(HYDROPATH);
      if (viscous)
	{
	  if ( b == 2.4 ) subfolder = "/hydro/moreland/fKLNb2.4etas0.08/";
	  else if ( b == 7.5 ) subfolder = "/hydro/moreland/fKLNb7.5etas0.08/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      else
	{
	  if ( b == 2.4 ) subfolder = "/hydro/moreland/fKLNb2.4ideal/";
	  else if ( b == 7.5 ) subfolder = "/hydro/moreland/fKLNb7.5ideal/";
	  else 
	    {
	      cout << "[HydroSetup::readHydroData]:ERROR: no hydro data for given impact parameter " << b << endl;
	      cout << " available impact parameters are b=2.4, and 7.5 fm." << endl;
	      exit(1);
	    }
	}
      if (envPath != 0 && *envPath != '\0') 
	{
	  int i = 0;
	  while (*(envPath+i) != '\0') path += *(envPath+(i++));
	  path += subfolder;
	}
      else path = ".."+subfolder;
      cout << "Path to hydro data: " << path << endl;
      
      ifstream fin;
      ostringstream temp;
      string istr;
      double x, y, T, vx, vy, tau, epsilon, QGPfrac;
      HydroCell newCell;  
      int position=0;
      // position = iy+iymax*(ix+ixmax*itau) y is the INNER LOOP here
      
      // open file with hydro data to read in:
        
      file[0] = path + "HydroEvol.dat";
      
      fin.open(file[0].c_str(),ios::in);
      if(!fin)
	{
	  cerr << "[HydroSetup::readHydroData]: ERROR: Unable to open file " << file[0] << endl;
	  exit(1);
	}
      
      int ik = 0;
      while ( !fin.eof() && ik<=ixmax*ixmax*itaumax )
	{
	  ik++;
	  fin >> x;
	  fin >> y;
	  fin >> tau;
	  fin >> epsilon;
	  fin >> T;
	  fin >> vx;
	  fin >> vy;
	    
	  T*=hbarc;
	  epsilon*=hbarc;
	    
	  QGPfrac=(epsilon-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
	  //cout << T << " " << QGPfrac << endl;
	  if (QGPfrac>1.) QGPfrac = 1;
	  else if (QGPfrac<0.) QGPfrac=0.;
	    
	  //cout << tau << " " << itaumax << " " << ik << " " << ixmax*ixmax*itaumax << endl;

	  newCell.QGPfrac(QGPfrac);
	  newCell.T(T);
	  newCell.T(T);
	  newCell.vx(vx);
	  newCell.vy(vy);
	  if (ik<=ixmax*ixmax*itaumax && (x<xmax && y<xmax)) lattice->push_back(newCell);
	}
      fin.close();
      // flow in x direction
      cout << " ok." << endl;
    }
}

// // find hydro values at current position doing interpolation
// // veta is read in and vz is returned and stored in the HydroInfo structure
// HydroInfo HydroSetup::getHydroValues2(double x, double y, double z, double t, double hydroXmax, double hydroZmax, double hydroTau0, 
// 				      double hydroDx, double hydroDz, double hydroDtau, int hydroWhichHydro, vector<HydroCell>* lattice)
// {
//   HydroInfo info;
//   double tau;
//   int position, positionTauUp;        // cell number in the hydro lattice ...
//   int positionZUp, positionZUpTauUp;  // the 16 corners of a 4D rectangle :) (in x,y,z and tau)
//   int positionXUp, positionXUpTauUp;  
//   int positionYUp, positionYUpTauUp;  
//   int positionXUpYUp;
//   int positionXUpYUpTauUp;            
//   int positionXUpYUpZUp;              
//   int positionXUpZUp;                 
//   int positionYUpZUp;                 
//   int positionXUpZUpTauUp;            
//   int positionYUpZUpTauUp;            
//   int positionXUpYUpZUpTauUp;         

//   double taufrac, zfrac;              // used for linear interpolation between tau steps 
//   double xfrac, yfrac;                // used for linear interpolation between tau steps 

//   double Tz, Tztau;                   // temperature
//   double Tx, Ty, Txy;                 
//   double Txz, Tyz, Txyz;              
//   double Tzt;                         
//   double Txt, Tyt, Txyt;              
//   double Txzt, Tyzt, Txyzt;           

//   double Qz, Qztau;                   // QGP fraction
//   double Qx, Qy, Qxy;                 
//   double Qxz, Qyz, Qxyz;              
//   double Qzt;                         
//   double Qxt, Qyt, Qxyt;              
//   double Qxzt, Qyzt, Qxyzt;           

//   double VXz, VXztau;                 // three components of the flow velocity
//   double VXx, VXy, VXxy;                 
//   double VXxz, VXyz, VXxyz;              
//   double VXzt;                         
//   double VXxt, VXyt, VXxyt;              
//   double VXxzt, VXyzt, VXxyzt;           

//   double VYz, VYztau;                
//   double VYx, VYy, VYxy;                 
//   double VYxz, VYyz, VYxyz;              
//   double VYzt;                         
//   double VYxt, VYyt, VYxyt;              
//   double VYxzt, VYyzt, VYxyzt;           
 
//   double VETAz, VETAztau;                
//   double VETAx, VETAy, VETAxy;                 
//   double VETAxz, VETAyz, VETAxyz;              
//   double VETAzt;                         
//   double VETAxt, VETAyt, VETAxyt;              
//   double VETAxzt, VETAyzt, VETAxyzt;           
 
//   double vetaY;                       // used to transform 3D hydro data
//   double veta, eta;                   // used in transform of 3D hydro data

//   int ix, iy, iz, itau;               // parton's position in cell coordinates
//   int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction

//   ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
//   if ( hydroWhichHydro > 2 && hydroWhichHydro < 5 ) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);

//   ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
//   iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
//   if (hydroWhichHydro > 2 && hydroWhichHydro < 5) 
//     iz = floor((hydroZmax+z)/hydroDz+0.0001);               // z-coordinate of the cell we are in now
//   // note that x and y run from -hydroXmax to +hydroXmax
//   // and ix and iy from 0 to 2*hydroXmax
//   // hence the (hydroXmax+x or y) for both
  
//   if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
//   itau = floor((tau-hydroTau0)/hydroDtau+0.0001);     	  
//   taufrac = (tau-hydroTau0)/hydroDtau-itau;
  
//   //cout << "x=" << x << " y=" << y << " z=" << z << " ix=" << ix << " iy=" << iy << " iz=" << iz << endl;

//   if ( hydroWhichHydro < 3 || hydroWhichHydro > 4 )                            // for 2D hydro
//     {
//       xfrac = x/hydroDx-(ix)+ixmax/2.;
//       yfrac = y/hydroDx-(iy)+ixmax/2.;
//       position=ix+ixmax*(iy+ixmax*itau);                // position of the parton on the bg lattice (cell number)
//       positionTauUp=ix+ixmax*(iy+ixmax*(itau+1));   
//       positionXUp=(ix+1)+ixmax*(iy+ixmax*(itau));	     
//       positionYUp=(ix)+ixmax*((iy+1)+ixmax*(itau));
//       positionXUpTauUp=(ix+1)+ixmax*(iy+ixmax*(itau+1));	     
//       positionYUpTauUp=(ix)+ixmax*((iy+1)+ixmax*(itau+1));
//       positionXUpYUp=(ix+1)+ixmax*((iy+1)+ixmax*((itau)));	     
//       positionXUpYUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*((itau+1)));
//     }
//   else if ( hydroWhichHydro > 2 && hydroWhichHydro < 5 )                      // for 3D hydro
//     {
//       xfrac = x/hydroDx-(ix)+ixmax/2.;
//       yfrac = y/hydroDx-(iy)+ixmax/2.;
//       zfrac = z/hydroDz-(iz)+izmax/2.;
//       position=ix+ixmax*(iy+ixmax*(iz+izmax*itau));     // position of the parton on the bg lattice (cell number)
//       positionXUp=(ix+1)+ixmax*(iy+ixmax*(iz+izmax*itau));	     
//       positionYUp=(ix)+ixmax*((iy+1)+ixmax*(iz+izmax*itau));
//       positionZUp=ix+ixmax*(iy+ixmax*((iz+1)+izmax*(itau)));
//       positionTauUp=ix+ixmax*(iy+ixmax*(iz+izmax*(itau+1)));
//       positionXUpZUp=(ix+1)+ixmax*(iy+ixmax*((iz+1)+izmax*itau));	     
//       positionYUpZUp=(ix)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*itau));
//       positionXUpTauUp=(ix+1)+ixmax*(iy+ixmax*(iz+izmax*(itau+1)));	     
//       positionYUpTauUp=(ix)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau+1)));
//       positionZUpTauUp=ix+ixmax*(iy+ixmax*((iz+1)+izmax*(itau+1)));
//       positionXUpYUp=(ix+1)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau)));	     
//       positionXUpYUpZUp=(ix+1)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau)));
//       positionXUpYUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau+1)));
//       positionXUpZUpTauUp=(ix+1)+ixmax*((iy)+ixmax*((iz+1)+izmax*(itau+1)));
//       positionYUpZUpTauUp=(ix)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau+1)));
//       positionXUpYUpZUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau+1)));
//     }
  
//   if ( ix < 0 || ix >= ixmax ) 
//     {
//       cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - x out of range x=" << x << ", ix=" << ix << ", ixmax=" << ixmax << endl;
//     }
//   if ( iy < 0 || iy >= ixmax )
//     {
//       cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - y out of range, y=" << y << ", iy="  << iy << ", iymax=" << ixmax << endl;
//     }
//   if ( (hydroWhichHydro > 2 && hydroWhichHydro < 5) && ( iz < 0 || iz >= izmax ) ) 
//     {
//       cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - z out of range, iz=" << iz << ", izmax=" << izmax << endl;
//     }
  
//   if (hydroWhichHydro  > 2) // for 3D hydro interpolate in 4D!
//     {
//       Tx=(1.-xfrac)*lattice->at(position).T()+xfrac*lattice->at(positionXUp).T();
//       Txy=(1.-xfrac)*lattice->at(positionYUp).T()+xfrac*lattice->at(positionXUpYUp).T();
//       Txz=(1.-xfrac)*lattice->at(positionZUp).T()+xfrac*lattice->at(positionXUpZUp).T();
//       Txyz=(1.-xfrac)*lattice->at(positionYUpZUp).T()+xfrac*lattice->at(positionXUpYUpZUp).T();
      
//       Ty=(1.-yfrac)*Tx+yfrac*Txy;
//       Tyz=(1.-yfrac)*Txz+yfrac*Txyz;
      
//       Tz=(1.-zfrac)*Ty+zfrac*Tyz;
      
//       Txt=(1.-xfrac)*lattice->at(positionTauUp).T()+xfrac*lattice->at(positionXUpTauUp).T();
//       Txyt=(1.-xfrac)*lattice->at(positionYUpTauUp).T()+xfrac*lattice->at(positionXUpYUpTauUp).T();
//       Txzt=(1.-xfrac)*lattice->at(positionZUpTauUp).T()+xfrac*lattice->at(positionXUpZUpTauUp).T();
//       Txyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).T()+xfrac*lattice->at(positionXUpYUpZUpTauUp).T();
      
//       Tyt=(1.-yfrac)*Txt+yfrac*Txyt;
//       Tyzt=(1.-yfrac)*Txzt+yfrac*Txyzt;
      
//       Tzt=(1.-zfrac)*Tyt+zfrac*Tyzt;
//       Tztau=(1.-zfrac)*Tyt+zfrac*Tyzt;
//       info.T=(1.-taufrac)*Tz+taufrac*Tztau; // get temperature at current pos.
//     }
//   else // for 2D hydro interpolate in tau only 
//     {
//       Tx=(1.-xfrac)*lattice->at(position).T()+xfrac*lattice->at(positionXUp).T();
//       Txy=(1.-xfrac)*lattice->at(positionYUp).T()+xfrac*lattice->at(positionXUpYUp).T();
//       Ty=(1.-yfrac)*Tx+yfrac*Txy;
      
//       Txt=(1.-xfrac)*lattice->at(positionTauUp).T()+xfrac*lattice->at(positionXUpTauUp).T();
//       Txyt=(1.-xfrac)*lattice->at(positionYUpTauUp).T()+xfrac*lattice->at(positionXUpYUpTauUp).T();
//       Tyt=(1.-yfrac)*Txt+yfrac*Txyt;
      
//       info.T=(1.-taufrac)*Ty+taufrac*Tyt; // get temperature at current pos.
    
//       Qx=(1.-xfrac)*lattice->at(position).QGPfrac()+xfrac*lattice->at(positionXUp).QGPfrac();
//       Qxy=(1.-xfrac)*lattice->at(positionYUp).QGPfrac()+xfrac*lattice->at(positionXUpYUp).QGPfrac();
//       Qy=(1.-yfrac)*Qx+yfrac*Qxy;
    
//       Qxt=(1.-xfrac)*lattice->at(positionTauUp).QGPfrac()+xfrac*lattice->at(positionXUpTauUp).QGPfrac();
//       Qxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpTauUp).QGPfrac();
//       Qyt=(1.-yfrac)*Qxt+yfrac*Qxyt;

//       info.QGPfrac=(1.-taufrac)*Qy+taufrac*Qyt; // get QGP fraction at current pos.
//     }
  
//   if (hydroWhichHydro > 2 && hydroWhichHydro < 5) // for 3D hydro interpolate in 4D!
//     {
//       Qx=(1.-xfrac)*lattice->at(position).QGPfrac()+xfrac*lattice->at(positionXUp).QGPfrac();
//       Qxy=(1.-xfrac)*lattice->at(positionYUp).QGPfrac()+xfrac*lattice->at(positionXUpYUp).QGPfrac();
//       Qxz=(1.-xfrac)*lattice->at(positionZUp).QGPfrac()+xfrac*lattice->at(positionXUpZUp).QGPfrac();
//       Qxyz=(1.-xfrac)*lattice->at(positionYUpZUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpZUp).QGPfrac();
      
//       Qy=(1.-yfrac)*Qx+yfrac*Qxy;
//       Qyz=(1.-yfrac)*Qxz+yfrac*Qxyz;
      
//       Qz=(1.-zfrac)*Qy+zfrac*Qyz;
      
//       Qxt=(1.-xfrac)*lattice->at(positionTauUp).QGPfrac()+xfrac*lattice->at(positionXUpTauUp).QGPfrac();
//       Qxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpTauUp).QGPfrac();
//       Qxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpZUpTauUp).QGPfrac();
//       Qxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpZUpTauUp).QGPfrac();
      
//       Qyt=(1.-yfrac)*Qxt+yfrac*Qxyt;
//       Qyzt=(1.-yfrac)*Qxzt+yfrac*Qxyzt;
      
//       Qzt=(1.-zfrac)*Qyt+zfrac*Qyzt;
//       Qztau=(1.-zfrac)*Qyt+zfrac*Qyzt;
//       info.QGPfrac=(1.-taufrac)*Qz+taufrac*Qztau; // get temperature at current pos.
//     }

//   // note that beta_x=vhydro_x/cosh(eta)=vhydro_x*tau/t and same for beta_y.
  
//   // get flow velocity 
//   if ( hydroWhichHydro < 3 || hydroWhichHydro > 4 ) // convert for 2D hydro
//     {
//       info.vx=(1.-taufrac)*lattice->at(position).vx()+taufrac*lattice->at(positionTauUp).vx(); 
//       info.vy=(1.-taufrac)*lattice->at(position).vy()+taufrac*lattice->at(positionTauUp).vy(); 
      
//       info.vx *= (tau/t);   // convert vx to betax by multiplying with 1/cosh(eta)=tau/t
//       info.vy *= (tau/t);   // same for vy
//       info.vz = z/t;        // longitudinal flow velocity for 2D Hydro (assuming Bjorken)
//     }
//   else if ( hydroWhichHydro > 2  && hydroWhichHydro < 5) // convert for 3D hydro incl. 4D interpolation
//     {
//       VXx=(1.-xfrac)*lattice->at(position).vx()+xfrac*lattice->at(positionXUp).vx();
//       VXxy=(1.-xfrac)*lattice->at(positionYUp).vx()+xfrac*lattice->at(positionXUpYUp).vx();
//       VXxz=(1.-xfrac)*lattice->at(positionZUp).vx()+xfrac*lattice->at(positionXUpZUp).vx();
//       VXxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vx()+xfrac*lattice->at(positionXUpYUpZUp).vx();
      
//       VXy=(1.-yfrac)*VXx+yfrac*VXxy;
//       VXyz=(1.-yfrac)*VXxz+yfrac*VXxyz;
      
//       VXz=(1.-zfrac)*VXy+zfrac*VXyz;
      
//       VXxt=(1.-xfrac)*lattice->at(positionTauUp).vx()+xfrac*lattice->at(positionXUpTauUp).vx();
//       VXxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vx()+xfrac*lattice->at(positionXUpYUpTauUp).vx();
//       VXxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vx()+xfrac*lattice->at(positionXUpZUpTauUp).vx();
//       VXxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vx()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vx();
      
//       VXyt=(1.-yfrac)*VXxt+yfrac*VXxyt;
//       VXyzt=(1.-yfrac)*VXxzt+yfrac*VXxyzt;
      
//       VXzt=(1.-zfrac)*VXyt+zfrac*VXyzt;
//       VXztau=(1.-zfrac)*VXyt+zfrac*VXyzt;
//       info.vx=(1.-taufrac)*VXz+taufrac*VXztau; // get vx at current pos.
      
//       VYx=(1.-xfrac)*lattice->at(position).vy()+xfrac*lattice->at(positionXUp).vy();
//       VYxy=(1.-xfrac)*lattice->at(positionYUp).vy()+xfrac*lattice->at(positionXUpYUp).vy();
//       VYxz=(1.-xfrac)*lattice->at(positionZUp).vy()+xfrac*lattice->at(positionXUpZUp).vy();
//       VYxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vy()+xfrac*lattice->at(positionXUpYUpZUp).vy();
      
//       VYy=(1.-yfrac)*VYx+yfrac*VYxy;
//       VYyz=(1.-yfrac)*VYxz+yfrac*VYxyz;
      
//       VYz=(1.-zfrac)*VYy+zfrac*VYyz;
//       VYxt=(1.-xfrac)*lattice->at(positionTauUp).vy()+xfrac*lattice->at(positionXUpTauUp).vy();
//       VYxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vy()+xfrac*lattice->at(positionXUpYUpTauUp).vy();
//       VYxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vy()+xfrac*lattice->at(positionXUpZUpTauUp).vy();
//       VYxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vy()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vy();
      
//       VYyt=(1.-yfrac)*VYxt+yfrac*VYxyt;
//       VYyzt=(1.-yfrac)*VYxzt+yfrac*VYxyzt;
      
//       VYzt=(1.-zfrac)*VYyt+zfrac*VYyzt;
//       VYztau=(1.-zfrac)*VYyt+zfrac*VYyzt;
//       info.vy=(1.-taufrac)*VYz+taufrac*VYztau; // get vy at current pos.
      
//       VETAx=(1.-xfrac)*lattice->at(position).vz()+xfrac*lattice->at(positionXUp).vz();
//       VETAxy=(1.-xfrac)*lattice->at(positionYUp).vz()+xfrac*lattice->at(positionXUpYUp).vz();
//       VETAxz=(1.-xfrac)*lattice->at(positionZUp).vz()+xfrac*lattice->at(positionXUpZUp).vz();
//       VETAxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vz()+xfrac*lattice->at(positionXUpYUpZUp).vz();
      
//       VETAy=(1.-yfrac)*VETAx+yfrac*VETAxy;
//       VETAyz=(1.-yfrac)*VETAxz+yfrac*VETAxyz;
      
//       VETAz=(1.-zfrac)*VETAy+zfrac*VETAyz;
//       VETAxt=(1.-xfrac)*lattice->at(positionTauUp).vz()+xfrac*lattice->at(positionXUpTauUp).vz();
//       VETAxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vz()+xfrac*lattice->at(positionXUpYUpTauUp).vz();
//       VETAxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vz()+xfrac*lattice->at(positionXUpZUpTauUp).vz();
//       VETAxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vz()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vz();
      
//       VETAyt=(1.-yfrac)*VETAxt+yfrac*VETAxyt;
//       VETAyzt=(1.-yfrac)*VETAxzt+yfrac*VETAxyzt;
      
//       VETAzt=(1.-zfrac)*VETAyt+zfrac*VETAyzt;
//       VETAztau=(1.-zfrac)*VETAyt+zfrac*VETAyzt;
//       veta=(1.-taufrac)*VETAz+taufrac*VETAztau; // get veta at current pos.
      
//       vetaY = 0.5 * log ( (1.+veta)/(1.-veta) );
//       eta = 0.5 * log ( (t+z)/(t-z) );
//       info.vx *= cosh(vetaY)/cosh(vetaY+eta); 
//       info.vy *= cosh(vetaY)/cosh(vetaY+eta); 
//       info.vz = (t*(veta)+z)/(t+(veta)*z);
//       info.veta = veta;
//     }
//   return info;
// }

// find hydro values at current position doing interpolation
// veta is read in and vz is returned and stored in the HydroInfo structure
HydroInfo HydroSetup::getHydroValues(double x, double y, double z, double t, double hydroXmax, double hydroZmax, double hydroTauMax,
				     double hydroTau0, double hydroDx, double hydroDz, double hydroDtau, int hydroWhichHydro, 
				     int fixedDistribution, vector<HydroCell>* lattice, bool trackHistory)
{
  //  cout << "getHydroValues" << endl;
  //cout << "hydroXmax=" << hydroXmax << endl; 
  //cout << "hydroZmax=" << hydroZmax << endl; 
  //cout << "hydroDx" << hydroDx << endl; 
  //cout << "hydroDz" << hydroDz << endl;
  //cout << x << " " << y << " " << z << " " << t << endl;
  //cout << "hydroTau0 = " << hydroTau0 << endl;
  //cout << "hydroTauMax = " << hydroTauMax << endl;
  //cout << "hydroDtau = " << hydroDtau << endl;

   HydroInfo info;
  double tau;
  double tempx;
  int position, positionTauUp;        // cell number in the hydro lattice ...
  int positionZUp, positionZUpTauUp;  // the 16 corners of a 4D rectangle :) (in x,y,z and tau)
  int positionXUp, positionXUpTauUp;  
  int positionYUp, positionYUpTauUp;  
  int positionXUpYUp;
  int positionXUpYUpTauUp;            
  int positionXUpYUpZUp;              
  int positionXUpZUp;                 
  int positionYUpZUp;                 
  int positionXUpZUpTauUp;            
  int positionYUpZUpTauUp;            
  int positionXUpYUpZUpTauUp;         
  double xFull, yFull;                // x and y including the sign (till the end)

  double taufrac, zfrac;              // used for linear interpolation between tau steps 
  double xfrac, yfrac;                // used for linear interpolation between tau steps 

  double Tz, Tztau;                   // temperature
  double Tx, Ty, Txy;                 
  double Txz, Tyz, Txyz;              
  double Tzt;                         
  double Txt, Tyt, Txyt;              
  double Txzt, Tyzt, Txyzt;           

  double Qz, Qztau;                   // QGP fraction
  double Qx, Qy, Qxy;                 
  double Qxz, Qyz, Qxyz;              
  double Qzt;                         
  double Qxt, Qyt, Qxyt;              
  double Qxzt, Qyzt, Qxyzt;           

  double VXz, VXztau;                 // three components of the flow velocity
  double VXx, VXy, VXxy;                 
  double VXxz, VXyz, VXxyz;              
  double VXzt;                         
  double VXxt, VXyt, VXxyt;              
  double VXxzt, VXyzt, VXxyzt;           

  double VYz, VYztau;                
  double VYx, VYy, VYxy;                 
  double VYxz, VYyz, VYxyz;              
  double VYzt;                         
  double VYxt, VYyt, VYxyt;              
  double VYxzt, VYyzt, VYxyzt;           
 
  double VETAz, VETAztau;                
  double VETAx, VETAy, VETAxy;                 
  double VETAxz, VETAyz, VETAxyz;              
  double VETAzt;                         
  double VETAxt, VETAyt, VETAxyt;              
  double VETAxzt, VETAyzt, VETAxyzt;           
 
  double vetaY;                       // used to transform 3D hydro data
  double veta, eta;                   // used in transform of 3D hydro data

  int ix, iy, iz, itau;               // parton's position in cell coordinates
  int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction
  int itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);

  ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
  if ( hydroWhichHydro > 2 && hydroWhichHydro != 5 ) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
  
  if ( hydroWhichHydro == 5 ) // switch x and y for the Morland hydro data (cause y is the inner loop there)
    {
      tempx=x;
      x=y;
      y=tempx;
    }

  ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
  iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
 
  if (hydroWhichHydro > 2 && hydroWhichHydro != 5 && hydroWhichHydro != 6)
    {
      ixmax = static_cast<int>(hydroXmax/hydroDx+0.0001);   // only the positive x,y are saved
      iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
      xFull = x;
      yFull = y;
      x = abs(x);                                           // using the symmetry in the transverse plane 
      y = abs(y);                                           // we only need to store data in one quarter of the plane
      ix = floor((x)/hydroDx+0.0001);             // x-coordinate of the cell we are in now
      iy = floor((y)/hydroDx+0.0001);             // y-coordinate of the cell we are in now
    }
  if (hydroWhichHydro == 6)
    iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
  
  //cout << "izmax=" << izmax << endl;
  //cout << "ixmax=" << ixmax << endl;
  //cout << "size=" << ixmax*ixmax*izmax*itaumax << endl;

  // note that x and y run from -hydroXmax to +hydroXmax
  // and ix and iy from 0 to 2*hydroXmax
  // hence the (hydroXmax+x or y) for both
  
  if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
  if (fixedDistribution) tau = hydroTau0;               // if the distribution is fixed (no evolution) use tau=tau0 always
  itau = floor((tau-hydroTau0)/hydroDtau+0.0001);     	  
  taufrac = (tau-hydroTau0)/hydroDtau-itau;
  
  //cout << "x=" << x << " y=" << y << " z=" << z << " ix=" << ix << " iy=" << iy << " iz=" << iz << endl;
  //cout << "t=" << t << " tau=" << tau << " itau=" << itau << " itaumax=" << itaumax << endl;

  if ( hydroWhichHydro < 3 || hydroWhichHydro == 5 )    // for 2D hydro
    {
      xfrac = x/hydroDx-(ix)+ixmax/2.;
      yfrac = y/hydroDx-(iy)+ixmax/2.;
      position=ix+ixmax*(iy+ixmax*itau);                // position of the parton on the bg lattice (cell number)
      positionTauUp=ix+ixmax*(iy+ixmax*(itau+1));   
      positionXUp=(ix+1)+ixmax*(iy+ixmax*(itau));	     
      positionYUp=(ix)+ixmax*((iy+1)+ixmax*(itau));
      positionXUpTauUp=(ix+1)+ixmax*(iy+ixmax*(itau+1));	     
      positionYUpTauUp=(ix)+ixmax*((iy+1)+ixmax*(itau+1));
      positionXUpYUp=(ix+1)+ixmax*((iy+1)+ixmax*((itau)));	     
      positionXUpYUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*((itau+1)));
    }
  else if ( hydroWhichHydro > 2 && hydroWhichHydro != 5 ) // for 3D hydro
    {
      xfrac = x/hydroDx-(ix);                             // x runs from 0 to xmax
      yfrac = y/hydroDx-(iy);
      if ( hydroWhichHydro == 6 )                         // x runs over all space
	{
	  xfrac = x/hydroDx-(ix)+ixmax/2.;
	  yfrac = y/hydroDx-(iy)+ixmax/2.;
	}
      zfrac = z/hydroDz-(iz)+izmax/2.;
      position=ix+ixmax*(iy+ixmax*(iz+izmax*itau));     // position of the parton on the bg lattice (cell number)
      if (ix+1<ixmax) positionXUp=(ix+1)+ixmax*(iy+ixmax*(iz+izmax*itau));
      else positionXUp = position;
      if (iy+1<ixmax) positionYUp=(ix)+ixmax*((iy+1)+ixmax*(iz+izmax*itau));
      else positionYUp = position;
      if (iz+1<izmax) positionZUp=ix+ixmax*(iy+ixmax*((iz+1)+izmax*(itau)));
      else positionZUp = position;
      if (itau+1<itaumax) positionTauUp=ix+ixmax*(iy+ixmax*(iz+izmax*(itau+1)));
      else positionTauUp = position;
      if (ix+1<ixmax && iz+1<izmax) positionXUpZUp=(ix+1)+ixmax*(iy+ixmax*((iz+1)+izmax*itau));
      else if (ix+1>=ixmax && iz+1<izmax) positionXUpZUp = positionZUp;
      else if (ix+1<ixmax && iz+1>=izmax) positionXUpZUp = positionXUp;
      else if (ix+1>=ixmax && iz+1>=izmax) positionXUpZUp = position;
      if (iy+1<ixmax && iz+1<izmax) positionYUpZUp=(ix)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*itau));
      else if (iy+1>=ixmax && iz+1<izmax) positionYUpZUp = positionZUp;
      else if (iy+1<ixmax && iz+1>=izmax) positionYUpZUp = positionYUp;
      else if (iy+1>=ixmax && iz+1>=izmax) positionYUpZUp = position;
      if (ix+1<ixmax && itau+1<itaumax) positionXUpTauUp=(ix+1)+ixmax*(iy+ixmax*(iz+izmax*(itau+1)));
      else if (ix+1>=ixmax && itau+1<itaumax) positionXUpTauUp=positionTauUp;
      else if (ix+1<ixmax && itau+1>=itaumax) positionXUpTauUp=positionXUp;
      else if (ix+1>=ixmax && itau+1>=itaumax) positionXUpTauUp=position;
      if (iy+1<ixmax && itau+1<itaumax) positionYUpTauUp=(ix)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau+1)));
      else if (iy+1>=ixmax && itau+1<itaumax) positionYUpTauUp=positionTauUp;
      else if (iy+1<ixmax && itau+1>=itaumax) positionYUpTauUp=positionYUp;
      else if (iy+1>=ixmax && itau+1>=itaumax) positionYUpTauUp=position;
      if (iz+1<izmax && itau+1<itaumax) positionZUpTauUp=ix+ixmax*(iy+ixmax*((iz+1)+izmax*(itau+1)));
      else if (iz+1>=izmax && itau+1<itaumax) positionZUpTauUp=positionTauUp;
      else if (iz+1<izmax && itau+1>=itaumax) positionZUpTauUp=positionZUp;
      else if (iz+1>=izmax && itau+1>=itaumax) positionZUpTauUp=position;
      if (ix+1<ixmax && iy+1<ixmax) positionXUpYUp=(ix+1)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau)));
      else if (ix+1>=ixmax && iy+1<ixmax) positionXUpYUp=positionYUp;
      else if (ix+1<ixmax && iy+1>=ixmax) positionXUpYUp=positionXUp;
      else if (ix+1>=ixmax && iy+1>=ixmax) positionXUpYUp=position;
      if(ix+1<ixmax && iy+1<ixmax && iz+1<izmax) positionXUpYUpZUp=(ix+1)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau)));
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1<izmax) positionXUpYUpZUp=positionZUp;
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1>=izmax) positionXUpYUpZUp=position;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1<izmax) positionXUpYUpZUp=positionYUpZUp;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1>=izmax) positionXUpYUpZUp=positionYUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1<izmax) positionXUpYUpZUp=positionXUpZUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1>=izmax) positionXUpYUpZUp=positionXUp;
      else if (ix+1<ixmax && iy+1<ixmax && iz+1>=izmax) positionXUpYUpZUp=positionXUpYUp;
      if (ix+1<ixmax && iy+1<ixmax && itau+1<itaumax) positionXUpYUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*(iz+izmax*(itau+1)));
      else if (ix+1>=ixmax && iy+1>=ixmax && itau+1<itaumax) positionXUpYUpTauUp=positionTauUp;
      else if (ix+1>=ixmax && iy+1>=ixmax && itau+1>=itaumax) positionXUpYUpTauUp=position;
      else if (ix+1>=ixmax && iy+1<ixmax && itau+1<itaumax) positionXUpYUpTauUp=positionYUpTauUp;
      else if (ix+1>=ixmax && iy+1<ixmax && itau+1>=itaumax) positionXUpYUpTauUp=positionYUp;
      else if (ix+1<ixmax && iy+1>=ixmax && itau+1<itaumax) positionXUpYUpTauUp=positionXUpTauUp;
      else if (ix+1<ixmax && iy+1>=ixmax && itau+1>=itaumax) positionXUpYUpTauUp=positionXUp;
      else if (ix+1<ixmax && iy+1<ixmax && itau+1>=itaumax) positionXUpYUpTauUp=positionXUpYUp;
      if (ix+1<ixmax && iz+1<izmax && itau+1<itaumax) positionXUpZUpTauUp=(ix+1)+ixmax*((iy)+ixmax*((iz+1)+izmax*(itau+1)));
      else if (ix+1>=ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpZUpTauUp=positionTauUp;
      else if (ix+1>=ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpZUpTauUp=position;
      else if (ix+1>=ixmax && iz+1<izmax && itau+1<itaumax) positionXUpZUpTauUp=positionZUpTauUp;
      else if (ix+1>=ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpZUpTauUp=positionZUp;
      else if (ix+1<ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpZUpTauUp=positionXUpTauUp;
      else if (ix+1<ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpZUpTauUp=positionXUp;
      else if (ix+1<ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpZUpTauUp=positionXUpZUp;
      if (iy+1<ixmax && iz+1<izmax && itau+1<itaumax) positionYUpZUpTauUp=(ix)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau+1)));
      else if (iy+1>=ixmax && iz+1>=izmax && itau+1<itaumax) positionYUpZUpTauUp=positionTauUp;
      else if (iy+1>=ixmax && iz+1>=izmax && itau+1>=itaumax) positionYUpZUpTauUp=position;
      else if (iy+1>=ixmax && iz+1<izmax && itau+1<itaumax) positionYUpZUpTauUp=positionZUpTauUp;
      else if (iy+1>=ixmax && iz+1<izmax && itau+1>=itaumax) positionYUpZUpTauUp=positionZUp;
      else if (iy+1<ixmax && iz+1>=izmax && itau+1<itaumax) positionYUpZUpTauUp=positionYUpTauUp;
      else if (iy+1<ixmax && iz+1>=izmax && itau+1>=itaumax) positionYUpZUpTauUp=positionYUp;
      else if (iy+1<ixmax && iz+1<izmax && itau+1>=itaumax) positionYUpZUpTauUp=positionYUpZUp;
      if (ix+1<ixmax && iy+1<ixmax && iz+1<izmax && itau+1<itaumax) 
	positionXUpYUpZUpTauUp=(ix+1)+ixmax*((iy+1)+ixmax*((iz+1)+izmax*(itau+1)));
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=position;
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionTauUp;      
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionZUp;
      else if (ix+1>=ixmax && iy+1>=ixmax && iz+1<izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionZUpTauUp;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionYUp;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionYUpTauUp;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionYUpZUp;
      else if (ix+1>=ixmax && iy+1<ixmax && iz+1<izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionYUpZUpTauUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionXUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionXUpTauUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionXUpZUp;
      else if (ix+1<ixmax && iy+1>=ixmax && iz+1<izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionXUpZUpTauUp;
      else if (ix+1<ixmax && iy+1<ixmax && iz+1>=izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionXUpYUp;
      else if (ix+1<ixmax && iy+1<ixmax && iz+1>=izmax && itau+1<itaumax) positionXUpYUpZUpTauUp=positionXUpYUpTauUp;
      else if (ix+1<ixmax && iy+1<ixmax && iz+1<izmax && itau+1>=itaumax) positionXUpYUpZUpTauUp=positionXUpYUpZUp;
      //cout << "OK up to here" << endl;

//       cout << "position=" << position << endl;
//       cout << "position=" << position << " " << positionTauUp << " " << positionZUp << " " << positionZUpTauUp
//          << " " << positionXUp << " " << positionXUpTauUp << " " << positionYUp << " " << positionYUpTauUp 
//          << " " << positionXUpYUp << " " << positionXUpYUpTauUp << " " << positionXUpYUpZUp 
//          << " " << positionXUpZUp << " " << positionYUpZUp << " " << positionXUpZUpTauUp
//          << " " << positionYUpZUpTauUp << " " << positionXUpYUpZUpTauUp << endl;         
    }
  if (!trackHistory)
    {
      if ( ix < 0 || ix >= ixmax ) 
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - x out of range x=" << x << ", ix=" << ix << ", ixmax=" << ixmax << endl;
	  cout << "x=" << x << " y=" << y << " z=" << z << " ix=" << ix << " iy=" << iy << " iz=" << iz << endl;
	  cout << "t=" << t << " tau=" << tau << " itau=" << itau << " itaumax=" << itaumax << endl;
	}
      if ( iy < 0 || iy >= ixmax )
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - y out of range, y=" << y << ", iy="  << iy << ", iymax=" << ixmax << endl;
	  cout << "x=" << x << " y=" << y << " z=" << z << " ix=" << ix << " iy=" << iy << " iz=" << iz << endl;
	  cout << "t=" << t << " tau=" << tau << " itau=" << itau << " itaumax=" << itaumax << endl;
	}
      if ( itau < 0 || itau >= itaumax )
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - tau out of range, itau=" << itau << ", itaumax=" << itaumax << endl;
	}
      if ( (hydroWhichHydro > 2 && hydroWhichHydro != 5) && ( iz < 0 || iz >= izmax ) ) 
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - z out of range, iz=" << iz << ", izmax=" << izmax << endl;
	}
    }
  else
    {
      if ( ix < 0 || ix >= ixmax ) 
	{
	  info.vx = 0.; 
	  info.vy = 0.; 
	  info.vz = 0.;
	  info.veta = 0.;
	  info.T = 0.;
	  info.QGPfrac = 0.;
	  return info;
	}
      if ( iy < 0 || iy >= ixmax )
	{
	  info.vx = 0.; 
	  info.vy = 0.; 
	  info.vz = 0.;
	  info.veta = 0.;
	  info.T = 0.;
	  info.QGPfrac = 0.;
	  return info;
	}
      if ( itau < 0 || itau >= itaumax )
	{
	  info.vx = 0.; 
	  info.vy = 0.; 
	  info.vz = 0.;
	  info.veta = 0.;
	  info.T = 0.;
	  info.QGPfrac = 0.;
	  return info;
	}
      if ( (hydroWhichHydro > 2 && hydroWhichHydro != 5) && ( iz < 0 || iz >= izmax ) ) 
	{
	  info.vx = 0.; 
	  info.vy = 0.; 
	  info.vz = 0.;
	  info.veta = 0.;
	  info.T = 0.;
	  info.QGPfrac = 0.;
	  return info;
	}
    }
  if (hydroWhichHydro > 2 && hydroWhichHydro != 5) // for 3D hydro interpolate in 4D!
    {
      Tx=(1.-xfrac)*lattice->at(position).T()+xfrac*lattice->at(positionXUp).T();
      Txy=(1.-xfrac)*lattice->at(positionYUp).T()+xfrac*lattice->at(positionXUpYUp).T();
      Txz=(1.-xfrac)*lattice->at(positionZUp).T()+xfrac*lattice->at(positionXUpZUp).T();
      Txyz=(1.-xfrac)*lattice->at(positionYUpZUp).T()+xfrac*lattice->at(positionXUpYUpZUp).T();
      Ty=(1.-yfrac)*Tx+yfrac*Txy;
      Tyz=(1.-yfrac)*Txz+yfrac*Txyz;

      Tz=(1.-zfrac)*Ty+zfrac*Tyz;

      Txt=(1.-xfrac)*lattice->at(positionTauUp).T()+xfrac*lattice->at(positionXUpTauUp).T();
      Txyt=(1.-xfrac)*lattice->at(positionYUpTauUp).T()+xfrac*lattice->at(positionXUpYUpTauUp).T();
      Txzt=(1.-xfrac)*lattice->at(positionZUpTauUp).T()+xfrac*lattice->at(positionXUpZUpTauUp).T();
      Txyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).T()+xfrac*lattice->at(positionXUpYUpZUpTauUp).T();

      Tyt=(1.-yfrac)*Txt+yfrac*Txyt;
      Tyzt=(1.-yfrac)*Txzt+yfrac*Txyzt;
      
      Tzt=(1.-zfrac)*Tyt+zfrac*Tyzt;
      Tztau=(1.-zfrac)*Tyt+zfrac*Tyzt;
      info.T=(1.-taufrac)*Tz+taufrac*Tztau; // get temperature at current pos.
    }
  else // for 2D hydro interpolate in tau only 
    {
      Tx=(1.-xfrac)*lattice->at(position).T()+xfrac*lattice->at(positionXUp).T();
      Txy=(1.-xfrac)*lattice->at(positionYUp).T()+xfrac*lattice->at(positionXUpYUp).T();
      Ty=(1.-yfrac)*Tx+yfrac*Txy;
      
      Txt=(1.-xfrac)*lattice->at(positionTauUp).T()+xfrac*lattice->at(positionXUpTauUp).T();
      Txyt=(1.-xfrac)*lattice->at(positionYUpTauUp).T()+xfrac*lattice->at(positionXUpYUpTauUp).T();
      Tyt=(1.-yfrac)*Txt+yfrac*Txyt;
      
      info.T=(1.-taufrac)*Ty+taufrac*Tyt; // get temperature at current pos.

      Qx=(1.-xfrac)*lattice->at(position).QGPfrac()+xfrac*lattice->at(positionXUp).QGPfrac();
      Qxy=(1.-xfrac)*lattice->at(positionYUp).QGPfrac()+xfrac*lattice->at(positionXUpYUp).QGPfrac();
      Qy=(1.-yfrac)*Qx+yfrac*Qxy;
    
      Qxt=(1.-xfrac)*lattice->at(positionTauUp).QGPfrac()+xfrac*lattice->at(positionXUpTauUp).QGPfrac();
      Qxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpTauUp).QGPfrac();
      Qyt=(1.-yfrac)*Qxt+yfrac*Qxyt;

      info.QGPfrac=(1.-taufrac)*Qy+taufrac*Qyt; // get QGP fraction at current pos.
    }
  
  if (hydroWhichHydro > 2 && hydroWhichHydro != 5) // for 3D hydro interpolate in 4D!
    {
      Qx=(1.-xfrac)*lattice->at(position).QGPfrac()+xfrac*lattice->at(positionXUp).QGPfrac();
      Qxy=(1.-xfrac)*lattice->at(positionYUp).QGPfrac()+xfrac*lattice->at(positionXUpYUp).QGPfrac();
      Qxz=(1.-xfrac)*lattice->at(positionZUp).QGPfrac()+xfrac*lattice->at(positionXUpZUp).QGPfrac();
      Qxyz=(1.-xfrac)*lattice->at(positionYUpZUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpZUp).QGPfrac();
      
      Qy=(1.-yfrac)*Qx+yfrac*Qxy;
      Qyz=(1.-yfrac)*Qxz+yfrac*Qxyz;
      
      Qz=(1.-zfrac)*Qy+zfrac*Qyz;
      
      Qxt=(1.-xfrac)*lattice->at(positionTauUp).QGPfrac()+xfrac*lattice->at(positionXUpTauUp).QGPfrac();
      Qxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpTauUp).QGPfrac();
      Qxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpZUpTauUp).QGPfrac();
      Qxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).QGPfrac()+xfrac*lattice->at(positionXUpYUpZUpTauUp).QGPfrac();
      
      Qyt=(1.-yfrac)*Qxt+yfrac*Qxyt;
      Qyzt=(1.-yfrac)*Qxzt+yfrac*Qxyzt;
      
      Qzt=(1.-zfrac)*Qyt+zfrac*Qyzt;
      Qztau=(1.-zfrac)*Qyt+zfrac*Qyzt;
      info.QGPfrac=(1.-taufrac)*Qz+taufrac*Qztau; // get QGP fraction at current pos.
    }

  // note that beta_x=vhydro_x/cosh(eta)=vhydro_x*tau/t and same for beta_y.
  
  // get flow velocity 
  if ( hydroWhichHydro < 3 || hydroWhichHydro == 5 ) // convert for 2D hydro
    {
      VXx=(1.-xfrac)*lattice->at(position).vx()+xfrac*lattice->at(positionXUp).vx();
      VXxy=(1.-xfrac)*lattice->at(positionYUp).vx()+xfrac*lattice->at(positionXUpYUp).vx();
      VXy=(1.-yfrac)*VXx+yfrac*VXxy;

      VXxt=(1.-xfrac)*lattice->at(positionTauUp).vx()+xfrac*lattice->at(positionXUpTauUp).vx();
      VXxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vx()+xfrac*lattice->at(positionXUpYUpTauUp).vx();
      VXyt=(1.-yfrac)*VXxt+yfrac*VXxyt;

      info.vx=(1.-taufrac)*VXy+taufrac*VXyt; // get vx at current pos.

      VYx=(1.-xfrac)*lattice->at(position).vy()+xfrac*lattice->at(positionXUp).vy();
      VYxy=(1.-xfrac)*lattice->at(positionYUp).vy()+xfrac*lattice->at(positionXUpYUp).vy();
      VYy=(1.-yfrac)*VYx+yfrac*VYxy;

      VYxt=(1.-xfrac)*lattice->at(positionTauUp).vy()+xfrac*lattice->at(positionXUpTauUp).vy();
      VYxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vy()+xfrac*lattice->at(positionXUpYUpTauUp).vy();
      VYyt=(1.-yfrac)*VYxt+yfrac*VYxyt;

      info.vy=(1.-taufrac)*VYy+taufrac*VYyt; // get vy at current pos.
      
      info.vx *= (tau/t);   // convert vx to betax by multiplying with 1/cosh(eta)=tau/t
      info.vy *= (tau/t);   // same for vy
      info.vz = z/t;        // longitudinal flow velocity for 2D Hydro (assuming Bjorken)
    }
  else if ( hydroWhichHydro > 2 && hydroWhichHydro != 5 ) // convert for 3D hydro incl. 4D interpolation
    {
      VXx=(1.-xfrac)*lattice->at(position).vx()+xfrac*lattice->at(positionXUp).vx();
      VXxy=(1.-xfrac)*lattice->at(positionYUp).vx()+xfrac*lattice->at(positionXUpYUp).vx();
      VXxz=(1.-xfrac)*lattice->at(positionZUp).vx()+xfrac*lattice->at(positionXUpZUp).vx();
      VXxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vx()+xfrac*lattice->at(positionXUpYUpZUp).vx();
      
      VXy=(1.-yfrac)*VXx+yfrac*VXxy;
      VXyz=(1.-yfrac)*VXxz+yfrac*VXxyz;
      
      VXz=(1.-zfrac)*VXy+zfrac*VXyz;
      
      VXxt=(1.-xfrac)*lattice->at(positionTauUp).vx()+xfrac*lattice->at(positionXUpTauUp).vx();
      VXxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vx()+xfrac*lattice->at(positionXUpYUpTauUp).vx();
      VXxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vx()+xfrac*lattice->at(positionXUpZUpTauUp).vx();
      VXxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vx()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vx();
      
      VXyt=(1.-yfrac)*VXxt+yfrac*VXxyt;
      VXyzt=(1.-yfrac)*VXxzt+yfrac*VXxyzt;
      
      VXzt=(1.-zfrac)*VXyt+zfrac*VXyzt;
      VXztau=(1.-zfrac)*VXyt+zfrac*VXyzt;
      info.vx=(1.-taufrac)*VXz+taufrac*VXztau; // get vx at current pos.
      
      VYx=(1.-xfrac)*lattice->at(position).vy()+xfrac*lattice->at(positionXUp).vy();
      VYxy=(1.-xfrac)*lattice->at(positionYUp).vy()+xfrac*lattice->at(positionXUpYUp).vy();
      VYxz=(1.-xfrac)*lattice->at(positionZUp).vy()+xfrac*lattice->at(positionXUpZUp).vy();
      VYxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vy()+xfrac*lattice->at(positionXUpYUpZUp).vy();
      
      VYy=(1.-yfrac)*VYx+yfrac*VYxy;
      VYyz=(1.-yfrac)*VYxz+yfrac*VYxyz;
      
      VYz=(1.-zfrac)*VYy+zfrac*VYyz;
      VYxt=(1.-xfrac)*lattice->at(positionTauUp).vy()+xfrac*lattice->at(positionXUpTauUp).vy();
      VYxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vy()+xfrac*lattice->at(positionXUpYUpTauUp).vy();
      VYxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vy()+xfrac*lattice->at(positionXUpZUpTauUp).vy();
      VYxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vy()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vy();
      
      VYyt=(1.-yfrac)*VYxt+yfrac*VYxyt;
      VYyzt=(1.-yfrac)*VYxzt+yfrac*VYxyzt;
      
      VYzt=(1.-zfrac)*VYyt+zfrac*VYyzt;
      VYztau=(1.-zfrac)*VYyt+zfrac*VYyzt;
      info.vy=(1.-taufrac)*VYz+taufrac*VYztau; // get vy at current pos.
      
      VETAx=(1.-xfrac)*lattice->at(position).vz()+xfrac*lattice->at(positionXUp).vz();
      VETAxy=(1.-xfrac)*lattice->at(positionYUp).vz()+xfrac*lattice->at(positionXUpYUp).vz();
      VETAxz=(1.-xfrac)*lattice->at(positionZUp).vz()+xfrac*lattice->at(positionXUpZUp).vz();
      VETAxyz=(1.-xfrac)*lattice->at(positionYUpZUp).vz()+xfrac*lattice->at(positionXUpYUpZUp).vz();
      
      VETAy=(1.-yfrac)*VETAx+yfrac*VETAxy;
      VETAyz=(1.-yfrac)*VETAxz+yfrac*VETAxyz;
      
      VETAz=(1.-zfrac)*VETAy+zfrac*VETAyz;
      VETAxt=(1.-xfrac)*lattice->at(positionTauUp).vz()+xfrac*lattice->at(positionXUpTauUp).vz();
      VETAxyt=(1.-xfrac)*lattice->at(positionYUpTauUp).vz()+xfrac*lattice->at(positionXUpYUpTauUp).vz();
      VETAxzt=(1.-xfrac)*lattice->at(positionZUpTauUp).vz()+xfrac*lattice->at(positionXUpZUpTauUp).vz();
      VETAxyzt=(1.-xfrac)*lattice->at(positionYUpZUpTauUp).vz()+xfrac*lattice->at(positionXUpYUpZUpTauUp).vz();
      
      VETAyt=(1.-yfrac)*VETAxt+yfrac*VETAxyt;
      VETAyzt=(1.-yfrac)*VETAxzt+yfrac*VETAxyzt;
      
      VETAzt=(1.-zfrac)*VETAyt+zfrac*VETAyzt;
      VETAztau=(1.-zfrac)*VETAyt+zfrac*VETAyzt;
      veta=(1.-taufrac)*VETAz+taufrac*VETAztau; // get veta at current pos.
      
      if (hydroWhichHydro==3)
	{
	  vetaY = 0.5 * log ( (1.+veta)/(1.-veta) );
	  eta = 0.5 * log ( (t+z)/(t-z) );
	  info.vx *= cosh(vetaY)/cosh(vetaY+eta); 
	  info.vy *= cosh(vetaY)/cosh(vetaY+eta); 
	  info.vz = (t*(veta)+z)/(t+(veta)*z);
	  info.veta = veta;
	}
      else
	{
	  // for Jeon/Schenke hydro, veta is actually already vz!
	  info.vz = veta;
	  info.veta = veta;
	}
     
      if (xFull<0) info.vx*=(-1.);
      if (yFull<0) info.vy*=(-1.);
    }
//   cout << "T=" << info.T << endl; 
//   cout << "QGPfrac=" << info.QGPfrac << endl;
//   cout << "vx=" << info.vx << endl; 
//   cout << "vy=" << info.vy << endl; 
//   cout << "vz=" << info.vz << endl; 
//  cout << " leaving getHydroValues" << endl;
  return info;
}

//For interpolation of evolution files in tau-eta coordinates. Only the reading of MUSIC's 
//evolution_xyeta.dat file is implemented here. For simplicity, hydroZmax refers to MUSIC's 
//eta_size, and similarly for hydroDz; however, x, y, z, and t are as usual to stay compatible 
//with MARTINI:
HydroInfo HydroSetup::getHydroValuesTauEta(double x, double y, double z, double t, double hydroXmax, double hydroZmax, double hydroTauMax,
				     double hydroTau0, double hydroDx, double hydroDz, double hydroDtau, int hydroWhichHydro, 
				     int fixedDistribution, vector<HydroCell>* lattice, bool trackHistory)
{
  HydroInfo info;
  double tau, eta;

  int position[2][2][2][2];
  double T;
  double QGPfrac;
  double vx;
  double vy;
  double vz;

  int ixmax, ietamax, itaumax;          // Maximum cell number in x-, y-, eta-,  and tau-directions.
  int ix, iy, ieta, itau;               // Parton's position in cell coordinates

  double xfrac, yfrac, etafrac, taufrac;

  //cout << "lattice size=" << lattice->size() << endl;
  ixmax = static_cast<int>(1.+2.*hydroXmax/hydroDx+0.0001);
  ietamax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
  itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
  
  ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
  iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
 
  if(t*t>z*z){
    tau = sqrt(t*t-z*z);
    eta = 0.5*log((t+z)/(t-z));
  }
  else{
    tau = 0.;
    eta = 0.;
  }

  ieta = floor((hydroZmax+eta)/hydroDz+0.0001);
  itau = floor((tau-hydroTau0)/hydroDtau+0.0001);       

  xfrac = (x-( (double)ix*hydroDx-hydroXmax ))/hydroDx;
  yfrac = (y-( (double)iy*hydroDx-hydroXmax ))/hydroDx;
  etafrac = eta/hydroDz-(double)ieta+0.5*(double)ietamax;
  taufrac = (tau-hydroTau0)/hydroDtau-(double)itau;

  if ( ix < 0 || ix >= ixmax ) 
    {
      cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - x out of range x=" << x << ", ix=" << ix << ", ixmax=" << ixmax << endl;
      cout << "x=" << x << " y=" << y << " eta=" << eta << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
      cout << "t=" << t << " tau=" << tau << " itau=" << itau << " itaumax=" << itaumax << endl;
	}
      if ( iy < 0 || iy >= ixmax )
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - y out of range, y=" << y << ", iy="  << iy << ", iymax=" << ixmax << endl;
	  cout << "x=" << x << " y=" << y << " eta=" << eta << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
	  cout << "t=" << t << " tau=" << tau << " itau=" << itau << " itaumax=" << itaumax << endl;
	}
      if ( itau < 0 || itau >= itaumax )
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - tau out of range, itau=" << itau << ", itaumax=" << itaumax << endl;
	}
      if ( ieta < 0 || ieta >= ietamax  ) 
	{
	  cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - eta out of range, ieta=" << ieta << ", ietamax=" << ietamax << endl;
	}


  //The array of positions on the 4-dimensional rectangle:
  for(int ipx=0; ipx<2; ipx++){
    int px;
    if(ipx==0 || ix==ixmax-1 ) px = ix;
    else px = ix+1;
    for(int ipy=0; ipy<2; ipy++){
      int py;
      if(ipy==0 || iy==ixmax-1 ) py = iy;
      else py = iy+1;
      for(int ipeta=0; ipeta<2; ipeta++){
	int peta;
	if(ipeta==0 || ieta==ietamax-1 ) peta = ieta;
	else peta = ieta+1;
	for(int iptau=0; iptau<2; iptau++){
	  int ptau;
	  if(iptau==0 || itau==itaumax-1 ) ptau = itau;
	  else ptau = itau+1;
	  position[ipx][ipy][ipeta][iptau] = px+ixmax*(py+ixmax*(peta+ietamax*ptau));
	}
      }
    }
  }

  //And now, the interpolation:
  T = QGPfrac = vx = vy = vz = 0.;
  for(int ipx=0; ipx<2; ipx++){
    double xfactor;
    if(ipx==0) xfactor = 1.-xfrac;
    else xfactor = xfrac;
    for(int ipy=0; ipy<2; ipy++){
      double yfactor;
      if(ipy==0) yfactor = 1.-yfrac;
      else yfactor = yfrac;
      for(int ipeta=0; ipeta<2; ipeta++){
	double etafactor;
	if(ipeta==0) etafactor = 1.-etafrac;
	else etafactor = etafrac;
	for(int iptau=0; iptau<2; iptau++){
	  double taufactor;
	  if(iptau==0) taufactor = 1.-taufrac;
	  else taufactor = taufrac;

	  T += xfactor*yfactor*etafactor*taufactor*lattice->at(position[ipx][ipy][ipeta][iptau]).T();
	  QGPfrac += xfactor*yfactor*etafactor*taufactor*lattice->at(position[ipx][ipy][ipeta][iptau]).QGPfrac();
	  vx += xfactor*yfactor*etafactor*taufactor*lattice->at(position[ipx][ipy][ipeta][iptau]).vx();
	  vy += xfactor*yfactor*etafactor*taufactor*lattice->at(position[ipx][ipy][ipeta][iptau]).vy();
	  vz += xfactor*yfactor*etafactor*taufactor*lattice->at(position[ipx][ipy][ipeta][iptau]).vz();
	}
      }
    }
  }

  info.T = T;
  info.QGPfrac = QGPfrac;
  info.vx = vx;
  info.vy = vy;
  info.vz = vz;
  info.veta = vz;
  
  return info;
}
