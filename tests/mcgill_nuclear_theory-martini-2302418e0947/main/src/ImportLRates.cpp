// ImportLRates.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in length dependent transition rates for radiative and elastic processes from files

#include "ImportLRates.h"

ImportLRates::ImportLRates()//constructor
{
}
ImportLRates::~ImportLRates(){}//destructor

// void ImportLRates::write_header(FILE *wfile, Gamma_info *dat, int profile, double T, int nlines)
// {
//   int tt=0;
//   fwrite ( (char *) &(dat->df) , sizeof ( double ) , 1 , wfile );
//   fwrite ( (char *) &(dat->da) , sizeof ( double ) , 1 , wfile );
//   fwrite ( (char *) &(dat->cf) , sizeof ( double ) , 1 , wfile );
//   fwrite ( (char *) &(dat->ca) , sizeof ( double ) , 1 , wfile );
//   fwrite ( (char *) &(dat->Nc) , sizeof ( int ) , 1 , wfile );
//   fwrite ( (char *) &(dat->Nf) , sizeof ( int ) , 1 , wfile );
//   fwrite ( (char *) &(dat->Bethe_Heitler) , sizeof ( int ) , 1 , wfile );
  
//   fwrite ( (char *) &(nlines), sizeof(int),1,wfile);
//   fwrite ( (char *) &(dat->n_k), sizeof(int),1,wfile);
  
//   fwrite ( (char *) &T, sizeof(double), 1,wfile);
  
//   fwrite ( (char *) &(profile), sizeof ( int ) , 1 , wfile );
//   fwrite( &tt, sizeof(int),1,wfile);  // filler
//   // Ideally, the profile information should be more complete than
//   // just a single integer...  I don't know yet what the most
//   // useful way to index profiles should be.
// }

// void ImportLRates::write_one_line(FILE *out, int nk, lineofdata line)
// {
//   /* what the data line will be talking about */
//   fwrite( &(line.p), sizeof(double), 1, out);
//   fwrite( &(line.L), sizeof(double), 1, out);

//   // what values of k are contained
//   fwrite( line.k, sizeof(double), nk, out);

//   /* the actual data */
//   fwrite( line.gam[0], sizeof(double), nk, out);
//   fwrite( line.gam[1], sizeof(double), nk, out);
//   fwrite( line.gam[2], sizeof(double), nk, out);
//   fwrite( line.gam[3], sizeof(double), nk, out);
//   fflush(out);
// }

int ImportLRates::read_raw(Gamma_info *dat, rawdata *raw)
{
  int i, j, tt;
  char fname[100];
  double *ptr1;

  FILE *in;
  string path     = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
      path += "/main/data/";
    }
  else path = "./data/";
  string pathandfile = path+"L2dat";
  const char* filename = pathandfile.c_str();
  cout << endl;
  cout << "Reading rates of length dependent inelastic collisions from file " << endl;
  cout << pathandfile.c_str() << " ... " << endl;
  size_t bytes_read;

  in = fopen ( filename , "rb" ); 
  if(!in)
    {
      cerr << "[ImportLRates::read_raw]: ERROR: Unable to open file " << filename << endl;
      exit(1);
    }
  
  // 
  fread ( (char *) (&dat->df) , sizeof ( double ) , 1 , in );
  fread ( (char *) (&dat->da) , sizeof ( double ) , 1 , in );
  fread ( (char *) (&dat->cf) , sizeof ( double ) , 1 , in );
  fread ( (char *) (&dat->ca) , sizeof ( double ) , 1 , in );
  fread ( (char *) (&dat->Nc) , sizeof ( int ) , 1 , in );
  fread ( (char *) (&dat->Nf) , sizeof ( int ) , 1 , in );
  fread ( (char *) (&dat->Bethe_Heitler) , sizeof ( int ) , 1 , in );

  // For these two things, "dat" is really just a placeholder
  fread ( (char *) &(raw->nlines), sizeof(int),1,in);
  fread ( (char *) &(raw->nk), sizeof(int),1,in);

  fread ( (char *) &(raw->T), sizeof(double),1,in);
  fread ( (char *) &(raw->profile), sizeof(int),1,in);
  fread ( (char *) &tt, sizeof(int),1,in);

  // Minimal sanity check
  if(raw->nlines< 0 || raw->nk< 0 || raw->nlines> 1e5 || raw->nk>1e3) {
    fprintf(stderr, "Data file is corrupted.\n");
    return -1;
  }

  cout << "before memory alloc..." << endl;

  /* memory allocation mess */
  const double size = raw->nlines * sizeof(lineofdata);
 
  raw->dat = new lineofdata[raw->nlines];
  //    raw->dat[0].k = malloc(5*raw->nlines*raw->nk * sizeof(double));
  raw->dat[0].k = new double[5*raw->nlines*raw->nk];

  cout << "after memory alloc..." << endl;


  raw->dat[0].gam[0]= raw->dat[0].k + 1* raw->nlines *raw->nk;
  raw->dat[0].gam[1]= raw->dat[0].k + 2* raw->nlines *raw->nk;
  raw->dat[0].gam[2]= raw->dat[0].k + 3* raw->nlines *raw->nk;
  raw->dat[0].gam[3]= raw->dat[0].k + 4* raw->nlines *raw->nk;
  for(i= 1; i< raw->nlines; i++) 
    {
      raw->dat[i].k= raw->dat[i-1].k + raw->nk;
      raw->dat[i].gam[0]= raw->dat[i-1].gam[0] + raw->nk;
      raw->dat[i].gam[1]= raw->dat[i-1].gam[1] + raw->nk;
      raw->dat[i].gam[2]= raw->dat[i-1].gam[2] + raw->nk;
      raw->dat[i].gam[3]= raw->dat[i-1].gam[3] + raw->nk;
    }
  
  // read the data
  for(i=0; i< raw->nlines; i++) 
    {
      /* what the data line will be talking about */
      fread( &(raw->dat[i].p), sizeof(double), 1, in);
      fread( &(raw->dat[i].L), sizeof(double), 1, in);
      fread( raw->dat[i].k, sizeof(double), raw->nk, in);
      
      /* actual data */
      fread( raw->dat[i].gam[0], sizeof(double), raw->nk, in);
      fread( raw->dat[i].gam[1], sizeof(double), raw->nk, in);
      fread( raw->dat[i].gam[2], sizeof(double), raw->nk, in);
      if(fread( raw->dat[i].gam[3], sizeof(double), raw->nk, in)<= 0)
	break;
    }
  
  //  printf("%d\n", i);
  
  raw->nlines= i;  /* the data file may not be full */
  return (raw->nlines == 0);  /* return value !=0 means error */
}

void ImportLRates::free_raw(rawdata *raw)
{
  free(raw->dat[0].k);
  free(raw->dat);
}

double ImportLRates::interpolate_k(rawdata *raw, int p, double k, int process)
{
  int k1,k2;
  double x;
  
  k1=k2=0;
  while(raw->dat[p].k[k2]< k && k2< raw->nk-1) { k1=k2; k2++; }
  
  if(k1==k2) // will only happen if k is smaller than leftmost point
    return( 0 );
  
  x= (k - raw->dat[p].k[k1]) / (raw->dat[p].k[k2] - raw->dat[p].k[k1]);
  return (
	  (1-x)*raw->dat[p].gam[process][k1]
	  + x*raw->dat[p].gam[process][k2]
	  );
}

// note that p,k are in GeV and L is in fermis.
double ImportLRates::use_raw ( double p , double k , double L,
			       rawdata *raw, int process)
/* Uses the raw lookup table and simple interpolation to get the value
   of dGamma/dk/g^4  at some value of p,L,k. */
/* L is in fermis */
/*   p1,p2,p3,p4 are: (p1,L1),(p1,L2),(p2,L1),(p2,L2) */
{
  double x, result;     /* fraction of way from corner of box. */
  int p1,p2,p3,p4;
  double r1,r2,r3;
  
  if( p< raw->dat[0].p || p> 1e5)   return 0;
  /* Out of range. */
  if ( ( process % 3 ) && ( k > p/2 ) )
    k = p - k;  /* Take advantage of symmetry in these cases */
  
  // the code assumes the lines of data are ordered according to
  // [increasing p][increasing L].

  // So, the first thing to do is set p1,p3 to the beginning of
  // the array of p's.

  p1=p3=0;
  while(raw->dat[p3].p< p && p3< raw->nlines-1) { p1=p3; p3++;  }
  /* p3 is always past p, or is last */
  while(p1>0 && raw->dat[p1].p==raw->dat[p1-1].p) p1--;
  while(p3>0 && raw->dat[p3].p==raw->dat[p3-1].p) p3--;

  // then, sort things out according to L.
  // I'll have p1=(smallest p, smallest L),
  //           p2=(smallest p, largest L),
  //           p3=(largest p, smallest L),
  //           p4=(largest p, largest L)

  p2=p1;
  while(raw->dat[p2].L< L && p2< raw->nlines-1 &&
        raw->dat[p2+1].p== raw->dat[p2].p) { p1=p2; p2++; }
  p4=p3;
  while(raw->dat[p4].L< L && p4< raw->nlines-1 &&
        raw->dat[p4+1].p== raw->dat[p4].p) { p3=p4; p4++; }

  // first, interpolate w/r to L
  r1= interpolate_k(raw,p1,k,process);
  if(p2!= p1) {  // actually have something to do
    r3= interpolate_k(raw,p2,k,process);
    x= (L - raw->dat[p1].L) / (raw->dat[p2].L - raw->dat[p1].L);
    r1= (1-x)*r1 + x*r3;
  }
  else // leftmost point. interpolate linearly to zero.
    r1*= L/raw->dat[p1].L;
  
  // next, interpolate in p, if needed.
  if(p3!= p1) {
    r2= interpolate_k(raw,p3,k,process);
    if(p4!= p3) {  // w/r to L first, for largest p.
      r3=interpolate_k(raw,p4,k,process);
      x= (L - raw->dat[p3].L) / (raw->dat[p4].L - raw->dat[p3].L);
      r2= (1-x)*r2 + x*r3;
    }
    else // leftmost point. interpolate linearly to zero.
      r2*= L/raw->dat[p3].L;
    x= (p - raw->dat[p1].p) / (raw->dat[p3].p - raw->dat[p1].p);
    r1= (1-x)*r1 + x*r2;
  }
  
  result= r1;
  
  // The "temperature" used here is raw->T which is written in the
  // output file; it doesn't have to be equal to the actual temperature
  // of the plasma but must agree with that used in "produce_raw()".
  
  if ( ABS(k) > 0.001 ) /* Avoid division by 0, should never get asked for */
    {
      k/= raw->T;
      p/= raw->T;
      
      switch ( process )
        {
        case 0:
          result /= k;
          if ( k < 20 )
            result /= 1 - exp(-k);
          if ( k > p - 20 )
            result /= 1 + exp(k-p);
          break;
        case 1:
          result /= p;
          if ( k < 20 )
            result /= 1 + exp(-k);
          if ( k > p - 20 )
            result /= 1 + exp(k-p);
          break;
        case 2:
          result /= k * (p-k) / p;
          if ( k < 20 )
            result /= 1 - exp(-k);
          if ( k > p - 20 )
	    result /= 1 - exp(k-p);
          break;
        case 3:
          result /= k;
          if ( k < 0 ) result = 0;
          if ( k > p-20 )
            result /= 1 + exp(k-p);
          break;
        }
      
    }
  
  return ( result );
}

void ImportLRates::init()
{
  cout << "Initializing length dependent radiative rates..."  << endl;
  read_raw( &dat, &Gam );
}

double ImportLRates::getRate(double p, double k, double L)
{
  return use_raw ( p , k , L,  &Gam , 0 );
}

double ImportLRates::getRate_em(double p, double k, double L)
{
  return use_raw ( p , k , L, &Gam , 3 );
}

double ImportLRates::getRate_gqq(double p, double k, double L)
{
  if(k<p/2) return use_raw ( p , k , L, &Gam , 1 );
  else return 0.;
}

double ImportLRates::getRate_ggg(double p, double k, double L)
{
  if(k<p/2) return use_raw ( p , k , L, &Gam , 2 );
  else return 0.;
}
