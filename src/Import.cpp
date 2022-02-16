// Import.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in transition rates for radiative and elastic processes from files

#include "Import.h"

Import::Import()//constructor
{
    //arrays for elastic rates:
    dGamma_qq   = new double[Nalphas*Nomega];
    dGamma_qg   = new double[Nalphas*Nomega];
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
        while (*(envPath+i) != '\0')
            path += *(envPath+(i++));
        path += "/main/data/";
    }
    else path = "./data/";
    string pathandfile;
    // Pure AMY
    if (rateSelector == 1)
        pathandfile = path+"radgamma_LO";
    else if (rateSelector == 2 )
        pathandfile = path+"radgamma_NLO";
    else if (rateSelector == 3 )
        pathandfile = path+"radgamma_NP";
    else if (rateSelector == 4 )
        pathandfile = path+"radgamma_LO-smallq";
    else
    {
        cout << "[ERROR in Import.cpp::read_table]: selected rate " << rateSelector << " not available." << endl;
        exit(1);
    } 
    const char* filename = pathandfile.c_str();
    cout << endl;
    cout << "[Import.cpp::read_table]: Reading rates of inelastic collisions from file " << endl;
    cout << "\t" << pathandfile.c_str() << " ... " << endl;
    size_t bytes_read;
    
    rfile = fopen ( filename , "rb" ); 
    bytes_read=fread ( (char *) (&dat->ddf)            , sizeof ( double ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->dda)            , sizeof ( double ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->dcf)            , sizeof ( double ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->dca)            , sizeof ( double ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->Nc)             , sizeof ( int ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->Nf)             , sizeof ( int ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->BetheHeitler)   , sizeof ( int ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->BDMPS)          , sizeof ( int ) , 1 , rfile );
    bytes_read=fread ( (char *) (&dat->include_gluons) , sizeof ( int ) , 1 , rfile );
    bytes_read=fread ( (char *) Gam->dGamma            , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->tau               , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->dGamma_gqq        , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->tau_gqq           , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->dGamma_ggg        , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->tau_ggg           , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->dGamma_em         , sizeof ( double ) , NP * NK , rfile );
    bytes_read=fread ( (char *) Gam->tau_em            , sizeof ( double ) , NP * NK , rfile );
    fclose ( rfile );
    cout << "[Import.cpp::read_table]: ok." << endl;
    dat->Nf=Nf;
    dat->dp=0.05;
    dat->p_max=20;
    dat->p_min=0;//exp(LogEmin); // set to zero because my array starts at zero!...
    //cout << "p_min=" << dat->p_min << endl;
  

    dat->n_p    = static_cast<int>( (1.001 + dat->p_max / dat->dp ) );	//	n_p = int(0.4 + 121 / 0.5) = 241
    dat->p_max  = dat->dp * dat->n_p;	//	p_max = 0.5 * 241 = 120.5
    dat->n_pmin = static_cast<int>( (1.001 + dat->p_min / dat->dp ) );	//	np_min = int(0.4 + 3.3 / 0.5) = 7
    dat->n_p   -= dat->n_pmin - 1;	//	n_p = 241 - (7 - 1) = 235
    dat->p_min  = dat->dp * dat->n_pmin;	//	p_min = 0.5 * 7 = 3.5
    
    dat->n_kmin = 1 + 2 * ( static_cast<int>( (2.0 / dat->dp ) ) );
    
    dat->k_min  = -dat->dp * dat->n_kmin;
    
    dat->n_k    = static_cast<int>( ( 8 + dat->p_max ) / ( 2 * dat->dp ) );
    dat->k_max  = 2 * dat->dp * ( dat->n_k - 1 ) + dat->k_min;

}

double Import::use_table ( double p , double k , double dGamma[NP][NK] , int which_kind )
/* 
    Uses the lookup table and simple interpolation to get the value
    of dGamma/dk dx at some value of p,k. 
    This works by inverting the relations between (p,k) and (n_p,n_k)
    used in building the table, to find out what continuous values
    of n_p, n_k should be considered; then linearly interpolates.
*/
{
    double a , b , result;     /* fraction of way from corner of box. */
    int    n_p , n_k; /* location of corner of box. */
    //if (p<4.01) cout << "small p=" << p << endl;
    if ( (p<4.01) || (p>46000) || (k<-12) || (k>p+12) ) return ( 0 );
    /* Out of range. */
    if ( ( which_kind % 3 ) && ( k > p/2 ) )
      k = p - k;  /* Take advantage of symmetry in these cases */
    a = 24.7743737154026 * log ( p * .2493765586034912718l );
    n_p = (int) a;
    a -= n_p;
    if (k < -2)
    {
        b = 60 + 5*k;
    }
    else if (k<-1)
    {
        b = 70 + 10*k;
    }
    else if (k<1)
    {
        b = 80 + 20*k;
    }
    else if (k<2)
    {
        b = 90 + 10*k;
    }
    else if ( k < p-2 )
    {
        b = 190 - 10*log ( 1.000670700260932956l / 
          ( 0.0003353501304664781l + (k-2) / (p-4) ) - 1 );
    }
    else if (k<p-1)
    {
        b = 290 + 10*(k-p);
    }
    else if (k<p+1)
    {
        b = 300 + 20*(k-p);
    }
    else if (k<p+2)
    {
        b = 310 + 10*(k-p);
    }
    else
    {
        b = 320 + 5*(k-p);
    }
    n_k = (int) b;
    b -= n_k;
    result = (1-a) * ( (1-b) * dGamma[n_p][n_k] + b * dGamma[n_p][n_k+1] )
             + a * ( (1-b) * dGamma[n_p+1][n_k] + b * dGamma[n_p+1][n_k+1] );
    //if (which_kind==3) Shuezhe Shi, August 6 2021 -->
    //  result /= ((1-a)*dGamma[n_p][290] + a*dGamma[n_p+1][290]); // Shuzhe -- scale by dGamma(k=p).
    //else
    //  result /= ((1-a)*dGamma[n_p][80] + a*dGamma[n_p+1][80]); // Shuzhe -- scale by dGamma(k=0).
    result /= ((1-a)*dGamma[n_p][80] + a*dGamma[n_p+1][80]); // Shuzhe -- scale by dGamma(k=0). ROUZ: comment this for a test
    
    if ( ABS(k) < 0.001 ) return result; /* Avoid division by 0, should never get asked for */
    
    switch ( which_kind )
    {
        case 0:
        {
            result /= k;
            if ( k < 20 )
                result /= 1 - exp(-k);
            if ( k > p - 20 )
                result /= 1 + exp(k-p);
            break;
        }
        case 1:
        {
            //result /= p; // Shuzhe -- scale by p.
            if ( k < 20 )
                result /= 1 + exp(-k);
            if ( k > p - 20 )
                result /= 1 + exp(k-p);
            break;
        }
        case 2:
        {
            result /= k * (p-k) / p;
            result /= 2./(1.-exp(-0.5*p)); // Shuzhe -- rescaled
            if ( k < 20 )
                result /= 1 - exp(-k);
            if ( k > p - 20 )
                result /= 1 - exp(k-p);
            break;
        }
        case 3:
        {
            //result /= k; // un-comment this for tests
            result /= k/sqrt(0.005*p); //un comment this to get MC behaviour.
            if ( k < 0 )
                result = 0;
            if ( k > p-20 )
                result /= 1 + exp(k-p);
            break;
        }
    }
    return result;
}

void Import::list_dGamma(Gamma_info *dat, dGammas *gam, double T, double beta, double cos_phi)
{
  double p , k;
  int    i_k;

  dat->dx_max = 0; 
  cout << "p" << " " << "k" << " " << "dGamma" << " " << "dGamma_gqq" << " " << "dGamma_ggg" << " " << "dGamma_em" << endl;

  double plist[11] = {4.01, 7.,  9.,  11.,  20.,  50.,  100.,  1000.,  5000.,  10000.,  40000.};
  string nlist[11] = {"4",  "7", "9", "11", "20", "50", "100", "1000", "5000", "10000", "40000"};
  for (int ifile=0; ifile<11; ifile++){
    string filename = "dGamma_p"+nlist[ifile]+".dat";
    p = plist[ifile];
    fstream fout1(filename,ios::out); 
    fout1.precision(12);
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
  }
  exit(1);
}

void Import::prep_dGamma( Gamma_info *dat, dGammas *gam, double T, double beta, double cos_phi)
  /* Assumes that the Gamma_info and dGammas have been gotten already.
   Evaluates the rates of emission of hard collinear particles,
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
    double e_scale, jacobian, gamma;
    double tot_dGamma, tot_dGammag, p, k;
    int    i_p, i_k, here;

    // changed by Qin
    //double compton;

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
            dGamma[here] = jacobian * use_table ( p, k, gam->dGamma, 0 );
            dGamma_em[here] = jacobian * use_table ( p, k, gam->dGamma_em, 0 ); // why not 3 in the last argument?
            
            tot_dGamma += dGamma[here];
            if ( 4 * i_k - 2 * dat->n_kmin <= 3 + i_p + dat->n_pmin )
            {
                dGamma_gqq[here] = Nf * jacobian * use_table ( p , k , gam->dGamma_gqq , 1 );
                dGamma_ggg[here] = jacobian * use_table ( p , k , gam->dGamma_ggg , 2 );
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
    string path = "";
    const char* MARTINIPATH = "MARTINIPATH";
    char* envPath = getenv(MARTINIPATH);
    if (envPath != 0 and *envPath != '\0') 
    {
        int i = 0;
        while (*(envPath+i) != '\0')
            path += *(envPath+(i++));
        path += "/main/data/elastic/";
    }
    else path = "./data/elastic/";

    ifstream fin;
    double as, omega;
    // position = iOmega+Nomega*(iAlphas)
        
    // open files with data to read in:
    file[0] = path + "logEnDtrqq";
    file[1] = path + "logEnDtrqg";
    
    cout << "Reading rates of elastic collisions from files:" << endl;
    cout << "\t 1. " <<file[0] << endl;
    cout << "\t 2. " <<file[1] << " ..." << endl;

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
              
        // Sangyong
        //cout << "ik = " << ik << endl;
        //cout << "as = " << as << endl;
        //cout << "omega = " << omega << endl;
        //cout << "dGamma_qq[" << ik << "] = " << dGamma_qq[ik] << endl;
        
        ik++;
    }
    fin.close();

    fin.open(file[1].c_str(),ios::in);
    if ( !fin )
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
        
        /*
        // Sangyong
        cout << "ik = " << ik << endl;
        cout << "as = " << as << endl;
        cout << "omega = " << omega << endl;
        cout << "dGamma_qg[" << ik << "] = " << dGamma_qg[ik] << endl;
        */
        
        ik++;
    }
    fin.close();
    
    cout << "[Import::readElasticRate]: ok." << endl;
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
        while (*(envPath+i) != '\0')
            path += *(envPath+(i++));
        path += "/main/data/elastic/";
    }
    else 
        path = "./data/elastic/";
      
    ifstream fin;
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
        // Sangyong:
        //cout << "ik = " << ik << endl;
        //cout << "as = " << as << endl;
        //cout << "omega = " << omega << endl;
        //cout << "q = " << q << endl;
        //cout << "dGamma_qq_q[" << ik << "] = " << dGamma_qq_q[ik] << endl;
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
    cout << "[Import::readElasticRateOmegaQ]: ok." << endl;
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
    
    if ( omega > 0 )
        iOmega = Nomega/2+floor((log(omega)+5)/omegaStep);
    else if ( omega < 0 )
        iOmega = Nomega/2-ceil((log(-omega)+5)/omegaStep)-1;
    iAlphas = floor((alpha_s-0.15)/alphaStep);
    
    position = iOmega+Nomega*(iAlphas);
    positionAlphaUp = iOmega+Nomega*(iAlphas+1);
    positionOmegaUp = iOmega+1+Nomega*(iAlphas);
    positionAlphaUpOmegaUp = iOmega+1+Nomega*(iAlphas+1);

    alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
    if (omega > 0) 
    {
        if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
            omegaFrac = (omega - (exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
                        /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
        else
            omegaFrac = 0.;
    }
    else if (omega < 0) 
    {
        if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
            omegaFrac = (-omega - (exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
                       /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
        else 
            omegaFrac = 0.;
    }

    if ( which_kind < 3 )
    {
        if (position > 0 and iAlphas < Nalphas and iOmega < Nomega) 
            rate = dGamma_qq[position];
        else 
            rate = 0.;
        if (iAlphas+1 < Nalphas) 
            rateAlphaUp = dGamma_qq[positionAlphaUp];
        else 
            rateAlphaUp = rate;
        if (iOmega+1 < Nomega)
            rateOmegaUp = dGamma_qq[positionOmegaUp];
        else 
            rateOmegaUp = rate;
        if (iAlphas < Nalphas and iOmega < Nomega)
            rateAlphaUpOmegaUp = dGamma_qq[positionAlphaUpOmegaUp];
        else 
            rateAlphaUpOmegaUp = rate;
    }
    else 
    {
        if (position > 0 and iAlphas<Nalphas and iOmega<Nomega)
            rate = dGamma_qg[position];
        else 
            rate = 0.;
        if (iAlphas+1 < Nalphas) 
            rateAlphaUp = dGamma_qg[positionAlphaUp];
        else 
            rateAlphaUp = rate;
        if (iOmega+1 < Nomega)
            rateOmegaUp = dGamma_qg[positionOmegaUp];
        else 
            rateOmegaUp = rate;
        if (iAlphas<Nalphas and iOmega<Nomega)
            rateAlphaUpOmegaUp = dGamma_qg[positionAlphaUpOmegaUp];
        else 
            rateAlphaUpOmegaUp = rate;
    }

    if (omega > 0)
    {
        rateOmegaAv        = (1.-omegaFrac)*rate        + omegaFrac*rateOmegaUp;
        rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
    }
    if (omega < 0)
    {
        rateOmegaAv        = (omegaFrac)*rate        + (1.-omegaFrac)*rateOmegaUp;
        rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp + (1.-omegaFrac)*rateAlphaUpOmegaUp;
    }      
    result = (1.-alphaFrac)*rateOmegaAv + alphaFrac*rateAlphaUpOmegaAv;
    
    return result; // leave out the * 9./4. for processes 3 and 4 to use the same envelope later
}

double Import::use_elastic_table_omega_q ( double omega , double q, double alpha_s , int which_kind )
  /* Uses the lookup table and simple interpolation to get the value
     of dGamma/domegadq at some value of omega, q, and alpha_s. */
{
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
    double rateOmegaAv, rateAlphaUpOmegaAv, rateQUpOmegaAv, rate2QUpOmegaAv, rateAlphaUpQUpOmegaAv;
    double rateAlphaUp2QUpOmegaAv;
    double rateQAv, rateAlphaUpQAv;
    double slope;
    double slopeAlphaUp;

    if (omega>0)
        iOmega = Nomega/2.+floor((log(omega)+5.)/omegaStep);
    else if (omega<0) 
        iOmega = Nomega/2.-ceil((log(-omega)+5.)/omegaStep)-1.;
    iQ = floor((log(q)+5.)/qStep+0.0001);
    iAlphas = floor((alpha_s-0.15)/alphaStep+0.0001);
    
    position                  = iQ   + Nq*(iOmega+Nomega*(iAlphas));
    positionAlphaUp           = iQ   + Nq*(iOmega+Nomega*(iAlphas+1));
    positionOmegaUp           = iQ   + Nq*(iOmega+1+Nomega*(iAlphas));
    positionQUp               = iQ+1 + Nq*(iOmega+Nomega*(iAlphas));
    position2QUp              = iQ+2 + Nq*(iOmega+Nomega*(iAlphas));
    positionAlphaUpOmegaUp    = iQ   + Nq*(iOmega+1+Nomega*(iAlphas+1));
    positionAlphaUpQUp        = iQ+1 + Nq*(iOmega+Nomega*(iAlphas+1));
    positionOmegaUpQUp        = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas));
    positionAlphaUpOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas+1));
    
    alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
    if (omega > 0.) 
    {
        if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
            omegaFrac = (omega - (exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
                        /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
        else 
            omegaFrac = 0.;
    }
    else if ( omega < 0. ) 
    {
        if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)!=exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
            omegaFrac = (-omega - (exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
                        /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))-exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
        else 
            omegaFrac = 0.;
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
        else 
            qFrac = 0.;
    }
    if ( which_kind == 1 or which_kind == 2 )
    {
        if (position >= 0 and iAlphas<Nalphas and iOmega<Nomega and iQ<Nq )
            rate = dGamma_qq_q[position];
        else 
            rate = 0.;
        if (iAlphas+1<Nalphas) 
            rateAlphaUp = dGamma_qq_q[positionAlphaUp];
        else 
            rateAlphaUp = rate;
        if (iOmega+1<Nomega)
            rateOmegaUp = dGamma_qq_q[positionOmegaUp];
        else
            rateOmegaUp = rate;
        if (iQ+1<Nq)
            rateQUp = dGamma_qq_q[positionQUp];
        else
            rateQUp = rate;
        if (iAlphas<Nalphas and iOmega<Nomega)
            rateAlphaUpOmegaUp = dGamma_qq_q[positionAlphaUpOmegaUp];
        else
            rateAlphaUpOmegaUp = rate;
        if (iAlphas<Nalphas and iQ<Nq)
            rateAlphaUpQUp = dGamma_qq_q[positionAlphaUpQUp];
        else 
            rateAlphaUpQUp = rate;
        if (iOmega+1<Nomega and iQ+1<Nq)
            rateOmegaUpQUp = dGamma_qq_q[positionOmegaUpQUp];
        else 
            rateOmegaUpQUp = rate;
        if (iAlphas<Nalphas and iOmega<Nomega and iQ<Nq)
            rateAlphaUpOmegaUpQUp = dGamma_qq_q[positionAlphaUpOmegaUpQUp];
        else 
            rateAlphaUpOmegaUpQUp = rate;
        if ( omega > 20 ) // used for extrapolation when the data points are too far apart
        { 
            if (iQ+2<Nq )
                rate2QUp = dGamma_qq_q[position2QUp];
            else 
                rate2QUp = rateQUp;
            if (iAlphas<Nalphas and iQ+2<Nq ) 
                rateAlphaUp2QUp = dGamma_qq_q[positionAlphaUpQUp+1];
            else 
                rateAlphaUp2QUp = rateAlphaUpQUp;
            if (iOmega<Nomega and iQ+2<Nq ) 
                rateOmegaUp2QUp = dGamma_qq_q[positionOmegaUpQUp+1];
            else 
                rateOmegaUp2QUp = rateOmegaUpQUp;
            if (iAlphas<Nalphas and iOmega<Nomega and iQ+2<Nq )
                rateAlphaUpOmegaUp2QUp = dGamma_qq_q[positionAlphaUpOmegaUpQUp+1];
            else 
                rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
        }
    }
    else if ( which_kind == 3 or which_kind == 4)
    {
        if (position > 0 and iAlphas < Nalphas and iOmega<Nomega and iQ<Nq )
            rate = dGamma_qg_q[position];
        else
            rate = 0.;
        if (iAlphas+1 < Nalphas)
            rateAlphaUp = dGamma_qg_q[positionAlphaUp];
        else
            rateAlphaUp = rate;
        if (iOmega+1 < Nomega)
            rateOmegaUp = dGamma_qg_q[positionOmegaUp];
        else
            rateOmegaUp = rate;
        if (iQ+1 < Nq)
            rateQUp = dGamma_qg_q[positionQUp];
        else 
            rateQUp = rate;
        if (iAlphas < Nalphas and iOmega < Nomega)
            rateAlphaUpOmegaUp = dGamma_qg_q[positionAlphaUpOmegaUp];
        else
            rateAlphaUpOmegaUp = rate;
        if (iAlphas < Nalphas and iQ < Nq)
            rateAlphaUpQUp = dGamma_qg_q[positionAlphaUpQUp];
        else
            rateAlphaUpQUp = rate;
        if (iOmega+1 < Nomega and iQ+1 < Nq)
            rateOmegaUpQUp = dGamma_qg_q[positionOmegaUpQUp];
        else
            rateOmegaUpQUp = rate;
        if (iAlphas < Nalphas and iOmega < Nomega and iQ < Nq)
            rateAlphaUpOmegaUpQUp = dGamma_qg_q[positionAlphaUpOmegaUpQUp];
        else
            rateAlphaUpOmegaUpQUp = rate;
        if (omega > 20.) // used for extrapolation when the data points are too far apart
        { 
            if (iQ+2 < Nq )
                rate2QUp = dGamma_qg_q[position2QUp];
            else
                rate2QUp = rateQUp;
            if (iAlphas < Nalphas and iQ+2<Nq )
                rateAlphaUp2QUp = dGamma_qg_q[positionAlphaUpQUp+1];
            else 
                rateAlphaUp2QUp = rateAlphaUpQUp;
            if (iOmega < Nomega and iQ+2 < Nq )
                rateOmegaUp2QUp = dGamma_qg_q[positionOmegaUpQUp+1];
            else 
                rateOmegaUp2QUp = rateOmegaUpQUp;
            if (iAlphas < Nalphas and iOmega < Nomega and iQ+2 < Nq )
                rateAlphaUpOmegaUp2QUp = dGamma_qg_q[positionAlphaUpOmegaUpQUp+1];
            else
                rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
        }
    }

    if (omega > 0. and omega <= 20.)
    {
        rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
        rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
        rateQUpOmegaAv = (1.-omegaFrac)*rateQUp + omegaFrac*rateOmegaUpQUp;
        rateAlphaUpQUpOmegaAv = (1.-omegaFrac)*rateAlphaUpQUp + omegaFrac*rateAlphaUpOmegaUpQUp;
    }
    else if (omega > 20.)
    {
        if ( rate != 0. and rateOmegaUp != 0. )//TODO : double comparison, needs to be fixed.
            rateOmegaAv = exp((1.-omegaFrac)*log(rate) + omegaFrac*log(rateOmegaUp));
        else if ( rate == 0. ) //TODO : come back here and fixed the double comparison.
            rateOmegaAv = rateOmegaUp;
        else if ( rateOmegaUp == 0. )
            rateOmegaAv = rate;
        else 
            rateOmegaAv = 0.;
        
        if ( rateAlphaUpOmegaUp != 0. and rateAlphaUp != 0. )//TODO : double comparison
            rateAlphaUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp) + omegaFrac*log(rateAlphaUpOmegaUp));
        else if ( rateAlphaUp == 0. )
            rateAlphaUpOmegaAv = rateAlphaUpOmegaUp;
        else if ( rateAlphaUpOmegaUp == 0. )
            rateAlphaUpOmegaAv = rateAlphaUp;
        else 
            rateAlphaUpOmegaAv = 0.;
        
        
        if ( rateOmegaUpQUp != 0. and rateQUp != 0. )//TODO : double comparison
        rateQUpOmegaAv = exp((1.-omegaFrac)*log(rateQUp) + omegaFrac*log(rateOmegaUpQUp));
        else if ( rateOmegaUpQUp == 0. ) // TODO : double comparison, needs fixing.
            rateQUpOmegaAv = rateQUp;
        else if ( rateQUp == 0. )
            rateQUpOmegaAv = rateOmegaUpQUp;
        else 
            rateQUpOmegaAv = 0.;
    
    
        if ( rateAlphaUpOmegaUpQUp != 0. and rateAlphaUpQUp != 0. )
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
        rateOmegaAv           = (omegaFrac)*rate           + (1.-omegaFrac)*rateOmegaUp;
        rateQUpOmegaAv        = (omegaFrac)*rateQUp        + (1.-omegaFrac)*rateOmegaUpQUp;
        rateAlphaUpOmegaAv    = (omegaFrac)*rateAlphaUp    + (1.-omegaFrac)*rateAlphaUpOmegaUp;
        rateAlphaUpQUpOmegaAv = (omegaFrac)*rateAlphaUpQUp + (1.-omegaFrac)*rateAlphaUpOmegaUpQUp;
    }      

    if ( omega > 20. ) // interpolate logs for large omega
    {
        if ( rateOmegaAv > 0. )
            rateQAv = exp((1.-qFrac)*log(rateOmegaAv) + qFrac*log(rateQUpOmegaAv));
        else // use extrapolation
        {
            slope = (log(rate2QUpOmegaAv)-log(rateQUpOmegaAv))/qStep;
            rateQAv = exp(log(rateQUpOmegaAv)-slope*((1.-qFrac)*qStep));
        }
        if ( rateAlphaUpOmegaAv > 0.)
            rateAlphaUpQAv = exp((1.-qFrac)*log(rateAlphaUpOmegaAv) + qFrac*log(rateAlphaUpQUpOmegaAv));
        else  // use extrapolation
        {
            slopeAlphaUp = (log(rateAlphaUp2QUpOmegaAv)-log(rateAlphaUpQUpOmegaAv))/qStep;
            rateAlphaUpQAv = exp(log(rateAlphaUpQUpOmegaAv)-slope*((1.-qFrac)*qStep));
        }
    }
    else // interpolate linearly for small omega
    {
        rateQAv        = (1.-qFrac)*rateOmegaAv        + qFrac*rateQUpOmegaAv;
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
    return use_table( p, k, Gam.dGamma, 0);
}

double Import::getRate_gqq(double p, double k)
{
    if ( k < p/2)
        return use_table( p, k, Gam.dGamma_gqq, 1);
    else
        return 0.;
}

double Import::getRate_ggg(double p, double k)
{
    if ( k < p/2 )
        return use_table( p, k, Gam.dGamma_ggg, 2);
    else return 0.;
}

double Import::getRate_em(double p, double k)
{
    return use_table( p, k, Gam.dGamma_em, 3); // last argument is 3! was bug in old program ...
}
double Import::getElasticRate(double p, double omega, double alpha_s, int process)
{
    if ((omega > 0 and omega < p) or (omega < 0 and omega > -p))
        return use_elastic_table ( omega, alpha_s, process);
    else 
        return 0.;
}

double Import::getElasticRateOmegaQ(double p, double omega, double q, double alpha_s, int process)
{
    if (q > sqrt(omega*omega) and ((omega > 0 and omega < p) or (omega < 0 and omega > -p)))
        return use_elastic_table_omega_q ( omega, q, alpha_s, process);
    else 
        return 0.;
}
