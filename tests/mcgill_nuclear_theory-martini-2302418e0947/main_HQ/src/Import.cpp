// Import.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in transition rates for radiative and elastic processes from files

#include "Import.h"

Import::Import()//constructor
{
  //arrays for elastic rates:
  dGamma_qq = new double[Nalphas*Nomega];
  dGamma_qg = new double[Nalphas*Nomega];
  dGamma_qq_q = new double[Nalphas*Nomega*Nq];
  dGamma_qg_q = new double[Nalphas*Nomega*Nq];
}
Import::~Import(){}//destructor

void Import::read_table ( Gamma_info *dat , dGammas *Gam, int rateSelector )
     /* Reads in the binary stored file of dGamma values. */
{
  FILE *rfile;
  string path     = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
      path += "/main_HQ/data/";
    }
  else path = "./data/";
  string pathandfile;
  if (rateSelector == 1)
    pathandfile = path+"radgamma";
  else if (rateSelector == 2 )
    pathandfile = path+"radgammaCut03";
  else if (rateSelector == 3 )
    pathandfile = path+"radgammaCut036";
  else if (rateSelector == 4 )
    pathandfile = path+"radgammaCut035";
  else if (rateSelector == 5 )
    pathandfile = path+"radgammaCutPiOver4";
  else if (rateSelector == 6 )
    pathandfile = path+"radgammaCut04PiOver4";
  else if (rateSelector == 7 )
    pathandfile = path+"radgammaCut038PiOver4";
  else
    {
      cout << "[ERROR in Import.cpp::read_table]: selected rate " << rateSelector << " not available." << endl;
      exit(1);
    } 
  const char* filename = pathandfile.c_str();
  cout << endl;
  cout << "Reading rates of inelastic collisions from file " << endl;
  cout << pathandfile.c_str() << " ... " << endl;
  size_t bytes_read;

  rfile = fopen ( filename , "rb" ); 
  bytes_read=fread ( (char *) (&dat->ddf) , sizeof ( double ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->dda) , sizeof ( double ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->dcf) , sizeof ( double ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->dca) , sizeof ( double ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->Nc) , sizeof ( int ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->Nf) , sizeof ( int ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->BetheHeitler) , sizeof ( int ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->BDMPS) , sizeof ( int ) , 1 , rfile );
  bytes_read=fread ( (char *) (&dat->include_gluons) , sizeof ( int ) , 1 , rfile );
  bytes_read=fread ( (char *) Gam->dGamma , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->tau , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->dGamma_gqq , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->tau_gqq , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->dGamma_ggg , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->tau_ggg , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->dGamma_em , sizeof ( double ) , NP * NK , rfile );
  bytes_read=fread ( (char *) Gam->tau_em , sizeof ( double ) , NP * NK , rfile );
  fclose ( rfile );
  cout << " ok." << endl;
  dat->Nf=Nf;
  dat->dp=0.05;
  dat->p_max=20;
  dat->p_min=0;//exp(LogEmin); // set to zero because my array starts at zero!...
  //cout << "p_min=" << dat->p_min << endl;
  

  dat->n_p = static_cast<int>( (1.001 + dat->p_max / dat->dp ) );	//	n_p = int(0.4 + 121 / 0.5) = 241
  dat->p_max = dat->dp * dat->n_p;	//	p_max = 0.5 * 241 = 120.5
  dat->n_pmin = static_cast<int>( (1.001 + dat->p_min / dat->dp ) );	//	np_min = int(0.4 + 3.3 / 0.5) = 7
  dat->n_p -= dat->n_pmin - 1;	//	n_p = 241 - (7 - 1) = 235
  dat->p_min = dat->dp * dat->n_pmin;	//	p_min = 0.5 * 7 = 3.5
 
  dat->n_kmin = 1 + 2 * ( static_cast<int>( (2.0 / dat->dp ) ) );

  dat->k_min = -dat->dp * dat->n_kmin;

  dat->n_k = static_cast<int>( ( 8 + dat->p_max ) / ( 2 * dat->dp ) );
  dat->k_max = 2 * dat->dp * ( dat->n_k - 1 ) + dat->k_min;

}

double Import::use_table ( double p , double k , double dGamma[NP][NK] , int which_kind )
  /* Uses the lookup table and simple interpolation to get the value
     of dGamma/dk dx at some value of p,k. */
  /* This works by inverting the relations between (p,k) and (n_p,n_k)
     used in building the table, to find out what continuous values
     of n_p, n_k should be considered; then linearly interpolates. */
{
  double a , b , result;     /* fraction of way from corner of box. */
  int    n_p , n_k; /* location of corner of box. */
  //if (p<4.01) cout << "small p=" << p << endl;
  if ( (p<4.01) | (p>46000) | (k<-12) | (k>p+12) ) return ( 0 );
  /* Out of range. */
  if ( ( which_kind % 3 ) && ( k > p/2 ) )
    k = p - k;  /* Take advantage of symmetry in these cases */
  a = 24.7743737154026 * log ( p * .2493765586034912718l );
  n_p = (int) a;
  a -= n_p;
  if ( k < 2 )
    {
      if ( k < -1 )
	{
	  if ( k < -2 )
	    b = 60 + 5*k;
	  else
	    b = 70+10*k;
	}
      else
	{
	  if ( k < 1 )
	    b = 80 + 20*k;
	  else
	    b = 90 + 10*k;
	}
    }
  else if ( k < p-2 )
    { /* This is that tricky middle ground. */
      b = 190 - 10*log ( 1.000670700260932956l / 
			 ( 0.0003353501304664781l + (k-2) / (p-4) ) - 1 );
    }
  else
    {
      if ( k < p+1 )
	{
	  if ( k < p-1 )
	    b = 290 + 10*(k-p);
	  else
	    b = 300 + 20*(k-p);
	}
      else
	{
	  if ( k < p+2 )
	    b = 310 + 10*(k-p);
	  else
	    b = 320 + 5*(k-p);
	}
    }
  n_k = (int) b;
  b -= n_k;
  result = (1-a) * ( (1-b) * dGamma[n_p][n_k] + b * dGamma[n_p][n_k+1] )
    +        a * ( (1-b) * dGamma[n_p+1][n_k] + b * dGamma[n_p+1][n_k+1] );
  if ( ABS(k) > 0.001 ) /* Avoid division by 0, should never get asked for */
    {
      switch ( which_kind )
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

void Import::list_dGamma ( Gamma_info *dat , dGammas *gam , 
			   double T , double beta , double cos_phi )
{
  double e_scale , jacobian , gamma , junk;
  double tot_dGamma , tot_dGammag , p , k;
  int    i_p , i_k , here;

  double compton;
  
  gamma = 1.0 / sqrt ( 1 - beta * beta );
  jacobian = 1 - beta * cos_phi;
  e_scale = gamma * jacobian / T;

  dat->dx_max = 0; 
  cout << "p" << " " << "k" << " " << "dGamma" << " " << "dGamma_gqq" << " " << "dGamma_ggg" << "dGamma_em" << endl;

//   for ( i_p = 0 ; i_p < NP ; i_p++ )
//     {
//       if ( i_p == 0 )
// 	p = 4.01;
//       else
//	p = p * 1.04119; /* spaced so 6---1000 is 0--127 */
//       for ( i_k = 0 ; i_k < NK ; i_k++ )
// 	{
// 	  if ( i_k < 50 )        /* spaced by 0.2  from -12 to -2 */
// 	    k = -12 + i_k * 0.2;
// 	  else if ( i_k < 60 )   /* spaced by 0.1  from -2  to -1 */
// 	    k = -2 + (i_k-50) * 0.1;
// 	  else if ( i_k < 100 )  /* spaced by 0.05 from -1  to +1 */
// 	    k = -1 + (i_k-60) * 0.05;
// 	  else if ( i_k < 110 )  /* spaced by 0.1  from +1  to +2 */
// 	    k = 1 + (i_k-100) * 0.1;
// 	  else if ( i_k < 270 )  /* spaced complicated, +2 to p-2 */
// 	    {
// 	      k = 0.1 * (i_k-190);
// 	      k = 2 + (p-4) * ( -0.0003353501304664781l
// 				+ 1.000670700260932956l / (1+exp(-k)) );
// 	    }
// 	  else if ( i_k < 280 )  /* spaced by 0.1  from p-2 to p-1 */
// 	    k = p - 2 + 0.1 * (i_k-270);
// 	  else if ( i_k < 320 )  /* spaced by 0.05 from p-1 to p+1 */
// 	    k = p + 0.05 * (i_k - 300);
// 	  else if ( i_k < 330 )  /* spaced by 0.1  from p+1 to p+2 */
// 	    k = p + 0.1 * (i_k - 310);
// 	  else                   /* spaced by 0.2  from p+2 to p+12 */
// 	    k = p + 0.2 * (i_k - 320);
	  
// 	  here = i_p + i_k * NP;
// 	  dGamma[here] = use_table ( p , k , gam->dGamma , 0 );
// 	  dGamma_em[here] = use_table ( p , k , gam->dGamma_em , 3 ); // 3 in the last argument! was a bug in old program


// 	  if ( 4 * i_k - 2 * dat->n_kmin <= 3 + i_p + dat->n_pmin )
// 	    {
// 	      dGamma_gqq[here] = use_table ( p , k , gam->dGamma_gqq , 1 );
// 	      dGamma_ggg[here] = use_table ( p , k , gam->dGamma_ggg , 2 );
// 	    }
// 	  else
// 	      {
// 		dGamma_gqq[here] = 0;
// 		dGamma_ggg[here] = 0;
// 	      }
//     }
      
  fstream fout1("dGamma4cut038PiOver4.dat",ios::out); 
      fout1.precision(12);  
  
      p = 4.01;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout1 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout1 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout1.close();

      fstream fout2("dGamma7cut038PiOver4.dat",ios::out); 
      fout2.precision(12);  
  
      p = 7.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout2 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout2 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout2.close();



      fstream fout2a("dGamma9cut038PiOver4.dat",ios::out); 
      fout2a.precision(12);  
  
      p = 9.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout2a << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout2a << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout2a.close();

      fstream fout3("dGamma11cut038PiOver4.dat",ios::out); 
      fout3.precision(12);  
  
      p = 11.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout3 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout3 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout3.close();

      fstream fout3b("dGammap20T0.235.dat",ios::out); 
      fout3b.precision(12);  
  
      p = 20./0.235;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout3b << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout3b << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout3b.close();




      fstream fout3a("dGamma50cut038PiOver4.dat",ios::out); 
      fout3a.precision(12);  
  
      p = 50.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout3a << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout3a << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout3a.close();


      fstream fout4("dGamma100cut038PiOver4.dat",ios::out); 
      fout4.precision(12);  
  
      p = 100.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout4 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout4 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout4.close();

      fstream fout5("dGamma1000cut038PiOver4.dat",ios::out); 
      fout5.precision(12);  
  
      p = 1000.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout5 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout5 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout5.close();

      fstream fout6("dGamma5000cut038PiOver4.dat",ios::out); 
      fout6.precision(12);  
  
      p = 5000.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout6 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout6 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout6.close();

      fstream fout7("dGamma10000cut038PiOver4.dat",ios::out); 
      fout7.precision(12);  
  
      p = 10000;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout7 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout7 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	}
      fout7.close();

      fstream fout8("dGamma40000cut038PiOver4.dat",ios::out); 
      fout8.precision(12);  
  
      p = 40000.;

      for ( i_k = 0 ; i_k < 40900 ; i_k++ )
	{
	  if ( i_k < 100 )    
	    k = -12 + i_k * 0.1;
	  else if ( i_k < 500 )
	    k = -2 + ( i_k-100 ) * 0.01;
	  else if ( i_k < 2480 )
	    k = 2 + (i_k-500) * 0.1;
	  else if ( i_k < 40900 ) 
	    k = 200 + (i_k-2480) * 1.;
	  
	  if (k<p/2.)
	    fout8 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << use_table ( p , k , gam->dGamma_gqq , 1 )
		  << " " << use_table ( p , k , gam->dGamma_ggg , 2 ) << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;
	  else
	    fout8 << scientific << p << " " << k << " " <<  use_table ( p , k , gam->dGamma , 0 ) << " "  << 0.
		  << " " << 0. << " " <<  use_table ( p , k , gam->dGamma_em , 3 ) << endl;

	}
      fout8.close();

      exit(1);
}

void Import::prep_dGamma ( Gamma_info *dat , dGammas *gam , 
			   double T , double beta , double cos_phi )
  /* Assumes that the Gamma_info and dGammas have been gotten already. */
  /* Evaluates the rates of emission of hard collinear particles,
     when in a medium of temperature T, moving with velocity beta
     WRT the lab frame, for a hard parton moving at angle phi WRT
     the boost direction of the medium.
     
     For applications, the boost direction is the radial direction,
     and phi is the azimuthal angle (cylindrical coordinates) of
     the hard parton.
     
     The basic formula is that the rest-frame energy is given by
     
     p_rest = p_lab (1 - beta * cos_phi) / sqrt( 1 - beta * beta )
     call 1/sqrt(1-beta*beta) = gamma.
     
     so we want to know
     
     dGamma(k_lab * gamma * (1-beta * cos_phi)) * (1-beta * cos_phi)
     
     where the last bit is a Jacobian for dkdt between frames. 
     ** so that is the lab-frame dGamma: dGamma/dkdt, with k and t both in the lab frame.

     When the flag PHOTON_22 is turned on, it adds the contribution
     of Compton and pair processes under the approximation that
     the photon momentum is the same as the particle's momentum
     (which is valid at leading-log and receives O(T) corrections)
     and that the particle emitting the photon is therefore
     lost.  We do this by modifying dGamma(p,k=p) or dGamma(p,k=p+-Delta)
     depending on whether p is an even or odd multiple of Delta.
     We use the Next-to-Leading-Log calculation here.
  */
{
  double e_scale , jacobian , gamma , junk;
  double tot_dGamma , tot_dGammag , p , k;
  int    i_p , i_k , here;

  /*** change by qin ***/
  double compton;

  gamma = 1.0 / sqrt ( 1 - beta * beta );
  jacobian = 1 - beta * cos_phi;
  e_scale = gamma * jacobian / T;

  dat->dx_max = 0;
  for ( i_p = 0 ; i_p < dat->n_p ; i_p++ )
    {
      tot_dGamma = 0;
      tot_dGammag = 0;
      for ( i_k = 0 ; i_k < dat->n_k ; i_k++ )
	{
	  p = dat->p_min + dat->dp * i_p;

	  k = dat->k_min + 2 * dat->dp * i_k;
	  
	  p *= e_scale;
	  k *= e_scale;
	  here = i_p + i_k * dat->n_p;
	  dGamma[here] = jacobian * 
	    use_table ( p , k , gam->dGamma , 0 );
	  dGamma_em[here] = jacobian *
	    use_table ( p , k , gam->dGamma_em , 0 ); // why not 3 in the last argument?

	  tot_dGamma += dGamma[here];
	  if ( 4 * i_k - 2 * dat->n_kmin <= 3 + i_p + dat->n_pmin )
	    {
	      dGamma_gqq[here] = Nf * jacobian *
		use_table ( p , k , gam->dGamma_gqq , 1 );
	      dGamma_ggg[here] = jacobian *
		use_table ( p , k , gam->dGamma_ggg , 2 );
	      tot_dGammag += dGamma_gqq[here] + dGamma_ggg[here];
	    }
	  else
	      {
		dGamma_gqq[here] = 0;
		dGamma_ggg[here] = 0;
	      }
	}

      if ( tot_dGamma > dat->dx_max )
	dat->dx_max = tot_dGamma;
      if ( dat->include_gluons && tot_dGammag > dat->dx_max )
	dat->dx_max = tot_dGammag;
    }
  dat->dx_max *= 2 * dat->dp;
  dat->dx_max = 0.8 / dat->dx_max;
  /* Value considered sufficient to avoid sign alternation problems. */
}

void Import::readElasticRate()
{
  cout.precision(10);
  string file[2];
  string path     = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
      path += "/main_HQ/data/elastic/";
    }
  else path = "./data/elastic/";

  ifstream fin;
  int iAlphas, iOmega;
  double as, omega;
  // position = iOmega+Nomega*(iAlphas)
      
  // open files with data to read in:
  file[0] = path + "logEnDtrqq";
  file[1] = path + "logEnDtrqg";
  
  cout << "Reading rates of elastic collisions from files" << endl;
  cout << file[0] << endl;
  cout << file[1] << " ..." << endl;

  fin.open(file[0].c_str(),ios::in);
  if(!fin)
    {
      cerr << "[Import::readElasticRate]: ERROR: Unable to open file " << file[0] << endl;
      exit(1);
    }
  int ik = 0;
  while ( !fin.eof() )
    {
      fin >> as;
      fin >> omega;
      fin >> dGamma_qq[ik];
      ik++;
    }
  fin.close();

  fin.open(file[1].c_str(),ios::in);
  if(!fin)
    {
      cerr << "[Import::readElasticRate]: ERROR: Unable to open file " << file[1] << endl;
      exit(1);
    }
  ik = 0;
  while ( !fin.eof() )
    {
      fin >> as;
      fin >> omega;
      fin >> dGamma_qg[ik];
      ik++;
    }
  fin.close();
  
  cout << " ok." << endl;
}

void Import::readElasticRateOmegaQ()
{
  cout.precision(10);
  string file[2];
  string path     = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
      path += "/main_HQ/data/elastic/";
    }
  else path = "./data/elastic/";
      
  ifstream fin;
  int iAlphas, iOmega;
  double as, omega, q;
      
  // open files with data to read in:
  file[0] = path + "logEnDqtrqq";
  file[1] = path + "logEnDqtrqg";

  cout << "Reading rates in omega and q for elastic collisions from files" << endl;
  cout << file[0] << endl;
  cout << file[1] << " ..." << endl;
  
  fin.open(file[0].c_str(),ios::in);
  if(!fin)
    {
      cerr << "[Import::readElasticRateOmegaQ]: ERROR: Unable to open file " << file[0] << endl;
      exit(1);
    }
  int ik = 0;
  while ( !fin.eof() )
    {
      fin >> as;
      fin >> omega;
      fin >> q;
      fin >> dGamma_qq_q[ik];
      ik++;
    }
  fin.close();

  fin.open(file[1].c_str(),ios::in);
  if(!fin)
    {
      cerr << "[Import::readElasticRateOmegaQ]: ERROR: Unable to open file " << file[1] << endl;
      exit(1);
    }
  ik = 0;
  while ( !fin.eof() )
    {
      fin >> as;
      fin >> omega;
      fin >> q;
      fin >> dGamma_qg_q[ik];
      ik++;
    }
  fin.close();
  
  cout << " ok." << endl;
}

double Import::use_elastic_table ( double omega , double alpha_s , int which_kind )
  /* Uses the lookup table and simple interpolation to get the value
     of dGamma/domega at some value of omega and alpha_s. */
{
  double result;
  double alphaFrac, omegaFrac;
  int iOmega;
  int iAlphas;
  int position, positionAlphaUp, positionOmegaUp, positionAlphaUpOmegaUp;
  double rate, rateAlphaUp, rateOmegaUp, rateAlphaUpOmegaUp;
  double rateOmegaAv, rateAlphaUpOmegaAv;

  if (omega>0) iOmega = Nomega/2+floor((log(omega)+5)/omegaStep);
  else if (omega<0) iOmega = Nomega/2-ceil((log(-omega)+5)/omegaStep)-1;
  iAlphas = floor((alpha_s-0.15)/alphaStep);

  position = iOmega+Nomega*(iAlphas);
  positionAlphaUp = iOmega+Nomega*(iAlphas+1);
  positionOmegaUp = iOmega+1+Nomega*(iAlphas);
  positionAlphaUpOmegaUp = iOmega+1+Nomega*(iAlphas+1);

  alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
  if (omega>0) 
    {
      if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (omega - (exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  else if (omega<0) 
    {
      if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (-omega - (exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }

  if ( which_kind < 3 )
    {
      if (position > 0 && iAlphas<Nalphas && iOmega<Nomega) rate = dGamma_qq[position];
      else rate = 0.;
      if (iAlphas+1<Nalphas) 
	rateAlphaUp = dGamma_qq[positionAlphaUp];
      else rateAlphaUp = rate;
      if (iOmega+1<Nomega)
	rateOmegaUp = dGamma_qq[positionOmegaUp];
      else rateOmegaUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega)
	rateAlphaUpOmegaUp = dGamma_qq[positionAlphaUpOmegaUp];
      else rateAlphaUpOmegaUp = rate;
   }
  else 
    {
      if (position > 0 && iAlphas<Nalphas && iOmega<Nomega) rate = dGamma_qg[position];
      else rate = 0.;
      if (iAlphas+1<Nalphas) 
	rateAlphaUp = dGamma_qg[positionAlphaUp];
      else rateAlphaUp = rate;
      if (iOmega+1<Nomega)
	rateOmegaUp = dGamma_qg[positionOmegaUp];
      else rateOmegaUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega)
	rateAlphaUpOmegaUp = dGamma_qg[positionAlphaUpOmegaUp];
      else rateAlphaUpOmegaUp = rate;
    }
  
  if (omega > 0)
    {
      rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
      rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
    }
  if (omega < 0)
    {
      rateOmegaAv = (omegaFrac)*rate + (1.-omegaFrac)*rateOmegaUp;
      rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp + (1.-omegaFrac)*rateAlphaUpOmegaUp;
    }      
  result = (1.-alphaFrac)*rateOmegaAv + alphaFrac*rateAlphaUpOmegaAv;

  return result; 
  // leave out the * 9./4. for processes 3 and 4 to use the same envelope later
}

double Import::use_elastic_table_omega_q ( double omega , double q, double alpha_s , int which_kind )
  /* Uses the lookup table and simple interpolation to get the value
     of dGamma/domegadq at some value of omega, q, and alpha_s. */
{
  cout.setf(ios::scientific);
  cout.precision(12);
  double result;
  double alphaFrac, omegaFrac, qFrac;
  int iOmega;
  int iAlphas;
  int iQ;
  int position, positionAlphaUp, positionOmegaUp, positionAlphaUpOmegaUp;
  int positionQUp, position2QUp, positionAlphaUpQUp, positionOmegaUpQUp, positionAlphaUpOmegaUpQUp;
  double rate, rateAlphaUp, rateOmegaUp, rateQUp, rate2QUp;
  double rateAlphaUpOmegaUp, rateAlphaUpQUp, rateAlphaUp2QUp, rateOmegaUpQUp, rateOmegaUp2QUp;
  double rateAlphaUpOmegaUpQUp, rateAlphaUpOmegaUp2QUp;
  double rateOmegaAv, rateAlphaUpOmegaAv, rateQUpOmegaAv, rate2QUpOmegaAv, rateAlphaUpQUpOmegaAv, rateAlphaUp2QUpOmegaAv;
  double rateQAv, rateAlphaUpQAv;
  double slope, slopeAlphaUp;

  if (omega>0) iOmega = Nomega/2.+floor((log(omega)+5.)/omegaStep);
  else if (omega<0) iOmega = Nomega/2.-ceil((log(-omega)+5.)/omegaStep)-1.;
  iQ = floor((log(q)+5.)/qStep+0.0001);
  iAlphas = floor((alpha_s-0.15)/alphaStep+0.0001);

  position = iQ + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUp = iQ + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas));
  positionQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas));
  position2QUp = iQ+2 + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUpOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas+1));
  positionAlphaUpQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas));
  positionAlphaUpOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas+1));

  alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
  if (omega>0.) 
    {
      if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (omega - (exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  else if (omega<0.) 
    {
      if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (-omega - (exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  if ( omega > 20. ) // interpolate the logs linearly for large omegas (since there the spectrum is dropping exp in q) 
    {
      qFrac = (log(q)-(floor((log(q)+5.)/qStep)*qStep-5.))/qStep;
    }
  else // linear interpolation in q for small omegas
    {
      if (exp(ceil((log(q)+5.)/qStep)*qStep-5.)!=exp(floor((log(q)+5.)/qStep)*qStep-5.))
	qFrac = (q - (exp(floor((log(q)+5.)/qStep)*qStep-5.)))
	  /((exp(ceil((log(q)+5.)/qStep)*qStep-5.))-exp(floor((log(q)+5.)/qStep)*qStep-5.));
      else qFrac = 0.;
    }
   if ( which_kind == 1 || which_kind == 2 )
    {
      if (position >= 0 && iAlphas<Nalphas && iOmega<Nomega && iQ<Nq ) rate = dGamma_qq_q[position];
      else rate = 0.;
      if (iAlphas+1<Nalphas) rateAlphaUp = dGamma_qq_q[positionAlphaUp];
      else rateAlphaUp = rate;
      if (iOmega+1<Nomega) rateOmegaUp = dGamma_qq_q[positionOmegaUp];
      else rateOmegaUp = rate;
      if (iQ+1<Nq) rateQUp = dGamma_qq_q[positionQUp];
      else rateQUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega) rateAlphaUpOmegaUp = dGamma_qq_q[positionAlphaUpOmegaUp];
      else rateAlphaUpOmegaUp = rate;
      if (iAlphas<Nalphas && iQ<Nq) rateAlphaUpQUp = dGamma_qq_q[positionAlphaUpQUp];
      else rateAlphaUpQUp = rate;
      if (iOmega+1<Nomega && iQ+1<Nq) rateOmegaUpQUp = dGamma_qq_q[positionOmegaUpQUp];
      else rateOmegaUpQUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega && iQ<Nq) rateAlphaUpOmegaUpQUp = dGamma_qq_q[positionAlphaUpOmegaUpQUp];
      else rateAlphaUpOmegaUpQUp = rate;
      if (omega>20) // used for extrapolation when the data points are too far apart
	{ 
	  if (iQ+2<Nq ) rate2QUp = dGamma_qq_q[position2QUp];
	  else rate2QUp = rateQUp;
	  if (iAlphas<Nalphas && iQ+2<Nq ) rateAlphaUp2QUp = dGamma_qq_q[positionAlphaUpQUp+1];
	  else rateAlphaUp2QUp = rateAlphaUpQUp;
	  if (iOmega<Nomega && iQ+2<Nq ) rateOmegaUp2QUp = dGamma_qq_q[positionOmegaUpQUp+1];
	  else rateOmegaUp2QUp = rateOmegaUpQUp;
	  if (iAlphas<Nalphas && iOmega<Nomega && iQ+2<Nq ) rateAlphaUpOmegaUp2QUp = dGamma_qq_q[positionAlphaUpOmegaUpQUp+1];
	  else rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  else if ( which_kind == 3 || which_kind == 4)
    {
      if (position > 0 && iAlphas<Nalphas && iOmega<Nomega && iQ<Nq ) rate = dGamma_qg_q[position];
      else rate = 0.;
      if (iAlphas+1<Nalphas) rateAlphaUp = dGamma_qg_q[positionAlphaUp];
      else rateAlphaUp = rate;
      if (iOmega+1<Nomega) rateOmegaUp = dGamma_qg_q[positionOmegaUp];
      else rateOmegaUp = rate;
      if (iQ+1<Nq) rateQUp = dGamma_qg_q[positionQUp];
      else rateQUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega) rateAlphaUpOmegaUp = dGamma_qg_q[positionAlphaUpOmegaUp];
      else rateAlphaUpOmegaUp = rate;
      if (iAlphas<Nalphas && iQ<Nq) rateAlphaUpQUp = dGamma_qg_q[positionAlphaUpQUp];
      else rateAlphaUpQUp = rate;
      if (iOmega+1<Nomega && iQ+1<Nq) rateOmegaUpQUp = dGamma_qg_q[positionOmegaUpQUp];
      else rateOmegaUpQUp = rate;
      if (iAlphas<Nalphas && iOmega<Nomega && iQ<Nq) rateAlphaUpOmegaUpQUp = dGamma_qg_q[positionAlphaUpOmegaUpQUp];
      else rateAlphaUpOmegaUpQUp = rate;
      if (omega>20.) // used for extrapolation when the data points are too far apart
	{ 
	  if (iQ+2<Nq ) rate2QUp = dGamma_qg_q[position2QUp];
	  else rate2QUp = rateQUp;
	  if (iAlphas<Nalphas && iQ+2<Nq ) rateAlphaUp2QUp = dGamma_qg_q[positionAlphaUpQUp+1];
	  else rateAlphaUp2QUp = rateAlphaUpQUp;
	  if (iOmega<Nomega && iQ+2<Nq ) rateOmegaUp2QUp = dGamma_qg_q[positionOmegaUpQUp+1];
	  else rateOmegaUp2QUp = rateOmegaUpQUp;
	  if (iAlphas<Nalphas && iOmega<Nomega && iQ+2<Nq ) rateAlphaUpOmegaUp2QUp = dGamma_qg_q[positionAlphaUpOmegaUpQUp+1];
	  else rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  
  if (omega > 0. && omega <= 20.)
    {
      rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
      rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
      rateQUpOmegaAv = (1.-omegaFrac)*rateQUp + omegaFrac*rateOmegaUpQUp;
      rateAlphaUpQUpOmegaAv = (1.-omegaFrac)*rateAlphaUpQUp + omegaFrac*rateAlphaUpOmegaUpQUp;
    }
  else if (omega>20.)
    {
      if ( rate != 0. && rateOmegaUp != 0. )
	rateOmegaAv = exp((1.-omegaFrac)*log(rate) + omegaFrac*log(rateOmegaUp));
      else if ( rate == 0. )
	rateOmegaAv = rateOmegaUp;
      else if ( rateOmegaUp == 0. )
	rateOmegaAv = rate;
      else 
	rateOmegaAv = 0.;

      
      if ( rateAlphaUpOmegaUp != 0. && rateAlphaUp != 0. )
	rateAlphaUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp) + omegaFrac*log(rateAlphaUpOmegaUp));
      else if ( rateAlphaUp == 0. )
	rateAlphaUpOmegaAv = rateAlphaUpOmegaUp;
      else if ( rateAlphaUpOmegaUp == 0. )
	rateAlphaUpOmegaAv = rateAlphaUp;
      else 
	rateAlphaUpOmegaAv = 0.;


      if ( rateOmegaUpQUp != 0. && rateQUp != 0. )
	rateQUpOmegaAv = exp((1.-omegaFrac)*log(rateQUp) + omegaFrac*log(rateOmegaUpQUp));
      else if ( rateOmegaUpQUp == 0. )
	rateQUpOmegaAv = rateQUp;
      else if ( rateQUp == 0. )
	rateQUpOmegaAv = rateOmegaUpQUp;
      else 
	rateQUpOmegaAv = 0.;


      if ( rateAlphaUpOmegaUpQUp != 0. && rateAlphaUpQUp != 0. )
	rateAlphaUpQUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUpQUp) + omegaFrac*log(rateAlphaUpOmegaUpQUp));
      else if ( rateAlphaUpQUp == 0. )
	rateAlphaUpQUpOmegaAv = rateAlphaUpOmegaUpQUp;
      else if ( rateAlphaUpOmegaUpQUp == 0. )
	rateAlphaUpQUpOmegaAv = rateAlphaUpQUp;
      else 
	rateAlphaUpQUpOmegaAv = 0.;

      rate2QUpOmegaAv = exp((1.-omegaFrac)*log(rate2QUp) + omegaFrac*log(rateOmegaUp2QUp));
      rateAlphaUp2QUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp2QUp) + omegaFrac*log(rateAlphaUpOmegaUp2QUp));
    }
  else if (omega < 0.)
    {
      rateOmegaAv = (omegaFrac)*rate + (1.-omegaFrac)*rateOmegaUp;
      rateQUpOmegaAv = (omegaFrac)*rateQUp + (1.-omegaFrac)*rateOmegaUpQUp;
      rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp + (1.-omegaFrac)*rateAlphaUpOmegaUp;
      rateAlphaUpQUpOmegaAv = (omegaFrac)*rateAlphaUpQUp + (1.-omegaFrac)*rateAlphaUpOmegaUpQUp;
    }      

  if ( omega > 20. ) // interpolate logs for large omega
    {
      if ( rateOmegaAv > 0. ) rateQAv = exp((1.-qFrac)*log(rateOmegaAv) + qFrac*log(rateQUpOmegaAv));
      else // use extrapolation
	{
	  slope = (log(rate2QUpOmegaAv)-log(rateQUpOmegaAv))/qStep;
	  rateQAv = exp(log(rateQUpOmegaAv)-slope*((1.-qFrac)*qStep));
	}
      if ( rateAlphaUpOmegaAv > 0.) rateAlphaUpQAv = exp((1.-qFrac)*log(rateAlphaUpOmegaAv) + qFrac*log(rateAlphaUpQUpOmegaAv));
      else  // use extrapolation
	{
	  slopeAlphaUp = (log(rateAlphaUp2QUpOmegaAv)-log(rateAlphaUpQUpOmegaAv))/qStep;
	  rateAlphaUpQAv = exp(log(rateAlphaUpQUpOmegaAv)-slope*((1.-qFrac)*qStep));
	}
    }
  else // interpolate linearly for small omega
    {
      rateQAv = (1.-qFrac)*rateOmegaAv + qFrac*rateQUpOmegaAv;
      rateAlphaUpQAv = (1.-qFrac)*rateAlphaUpOmegaAv + qFrac*rateAlphaUpQUpOmegaAv;
    }

  result = (1.-alphaFrac)*rateQAv + alphaFrac*rateAlphaUpQAv;
  
  return result;
  // the absolute normalization doesn't matter when sampling the shape - it only matters in "totalRate" etc.
}

void Import::show_dGamma()
{
  //T, beta, cos_phi
  list_dGamma( &dat , &Gam , 1. , 0. , 1. );
}

void Import::init(int rateSelector)
{
  read_table( &dat, &Gam, rateSelector );
  
  readElasticRate();
  readElasticRateOmegaQ();
}

double Import::getRate(double p, double k)
{
  return use_table ( p , k , Gam.dGamma , 0 );
}

double Import::getRate_em(double p, double k)
{
  return use_table ( p , k , Gam.dGamma_em , 3 ); // last argument is 3! was bug in old program ...
}

double Import::getRate_gqq(double p, double k)
{
  if(k<p/2) return use_table ( p , k , Gam.dGamma_gqq , 1 );
  else return 0.;
}

double Import::getRate_ggg(double p, double k)
{
  if(k<p/2) return use_table ( p , k , Gam.dGamma_ggg , 2 );
  else return 0.;
}

double Import::getElasticRate(double p, double omega, double alpha_s, int process)
{
  if ((omega > 0 && omega < p) || (omega < 0 && omega > -p))
    return use_elastic_table ( omega, alpha_s, process);
  else 
    return 0.;
}

double Import::getElasticRateOmegaQ(double p, double omega, double q, double alpha_s, int process)
{
  if (q > sqrt(omega*omega) && ((omega > 0 && omega < p) || (omega < 0 && omega > -p)))
    return use_elastic_table_omega_q ( omega, q, alpha_s, process);
  else 
    return 0.;
}
