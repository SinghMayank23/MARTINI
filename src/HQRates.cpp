// MS
// This file contains routines to read HQ energy loss data from files and functions 
// that return interpolated data at a given space-time point

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<math.h>
#include<sstream>

#include "HQRates.h"

using namespace std;

HQRates::HQRates()  // constructor
{
    HQlattice = new HQratecell* [135];
    for (int ii = 0; ii < 135; ii++) HQlattice[ii] = new HQratecell[81];
    Coalesence_Prob_lattice = new double[81]; 
    Coalesence_Meson_Prob_lattice = new double[81]; 
}

HQRates::~HQRates()   // destructor
{
    for (int ii =0; ii < 135; ii++) delete[] HQlattice[ii];
    delete[] HQlattice;
    delete[] Coalesence_Prob_lattice;
    delete[] Coalesence_Meson_Prob_lattice;
}

void HQRates::readHQRates(
    int scat_rad_flag, string rate_path)
{
    for (int ii = 0; ii < 135; ii++)
    for (int jj = 0; jj < 81 ; jj++)
    {
      HQlattice[ii][jj].A0 = 0.;
      HQlattice[ii][jj].B0 = 0.;
      HQlattice[ii][jj].B1 = 0.;
    }

    string filenameA0, filenameB0, filenameB1;
    if (scat_rad_flag == 0 || scat_rad_flag == 1)
    {
      filenameA0 = rate_path + "/A0_scattering.dat";
      filenameB0 = rate_path + "/B0_scattering.dat";
      filenameB1 = rate_path + "/B1_scattering.dat";
      ifstream fiA0, fiB0, fiB1;
      fiA0.open(filenameA0.c_str(), ios::in);
      fiB0.open(filenameB0.c_str(), ios::in);
      fiB1.open(filenameB1.c_str(), ios::in);
      if(!fiA0 || !fiB0 || !fiB1)
      {
          cerr << "[HQRates]: ERROR: Unable to open file " << endl;
          exit(1);
      }

      double A0, B0, B1;
      for (int ii = 0; ii < 135; ii++)
      for (int jj = 0; jj < 81 ; jj++)
      {
        fiA0 >> A0;
        fiB0 >> B0;
        fiB1 >> B1;
        HQlattice[ii][jj].A0 += A0*hbarc            ;//Converted to [GeV]
        HQlattice[ii][jj].B0 += B0*hbarc*hbarc*hbarc;//Converted to [Gev^3]
        HQlattice[ii][jj].B1 += B1*hbarc*hbarc*hbarc;
      }

      fiA0.close();
      fiB0.close();
      fiB1.close();
    }

    if (scat_rad_flag == 0 || scat_rad_flag == 2)
    {
      filenameA0 = rate_path + "/A0_radiation.dat";
      filenameB0 = rate_path + "/B0_radiation.dat";
      filenameB1 = rate_path + "/B1_radiation.dat";
      ifstream fiA0, fiB0, fiB1;
      fiA0.open(filenameA0.c_str(), ios::in);
      fiB0.open(filenameB0.c_str(), ios::in);
      fiB1.open(filenameB1.c_str(), ios::in);
      if(!fiA0 || !fiB0 || !fiB1)
      {
          cerr << "[HQRates]: ERROR: Unable to open file " << endl;
          exit(1);
      }

      double A0, B0, B1;
      for (int ii = 0; ii < 135; ii++)
      for (int jj = 0; jj < 81 ; jj++)
      {
        fiA0 >> A0;
        fiB0 >> B0;
        fiB1 >> B1;
        HQlattice[ii][jj].A0 += A0*hbarc            ;//Converted to [GeV]
        HQlattice[ii][jj].B0 += B0*hbarc*hbarc*hbarc;//Converted to [Gev^3]
        HQlattice[ii][jj].B1 += B1*hbarc*hbarc*hbarc;
      }

      fiA0.close();
      fiB0.close();
      fiB1.close();
    }
}

void HQRates::getTandPdependent_coeff(double T, double pHQ, double& A0, double& B0, double& B1)
{

// First axis of the table HQlattice corresponds to temperature ranging from 0.130 GeV to 0.800 GeV
// in increments of 0.005 GeV
// Second axis corresponds to momentum. The first value HQlattice[xx][0] corresponds to p = 0.05 GeV.
// Next one is 0.25 GeV and then subsequent ones are in increment of 0.25 GeV

    int iT_low, iT_high, ip_low, ip_high;
    double T_low, T_high, p_low, p_high, dp, dT;
    double A0_11, A0_12, A0_21, A0_22;
    double B0_11, B0_12, B0_21, B0_22;
    double B1_11, B1_12, B1_21, B1_22;

    dT = 0.005;
    if (T >= 0.130 && T <= 0.800)
    {
      iT_low  = (T - 0.130)/0.005;
      iT_high = iT_low + 1;
      T_low  = 0.130 + (double)iT_low * 0.005;
      T_high = 0.130 + (double)iT_high* 0.005;
    } else if (T < 0.130)
    {
      iT_low  = 0;
      iT_high = 1;
      T_low  = 0.130;
      T_high = 0.135;
    } else if (T > 0.800)
    {
      iT_low  = 133;
      iT_high = 134;
      T_low  = 0.795;
      T_high = 0.800;
    }

    if (pHQ >= 0.25 && pHQ <= 20.)
    {
      ip_low  = (pHQ - 0.25)/0.25 + 1;
      ip_high = ip_low + 1;
      p_low  = (double)ip_low * 0.25;
      p_high = (double)ip_high* 0.25;
      dp = 0.25;
    } else if (pHQ < 0.25)
    {
      ip_low  = 0;
      ip_high = 1;
      p_low  = 0.05;
      p_high = 0.25;
      dp = 0.2;
    } else if (pHQ > 20.0)
    {
      ip_low  = 79;
      ip_high = 80;
      p_low  = 19.75;
      p_high = 20.00;
      dp = 0.25;
    }

    A0_11 = HQlattice[iT_low ][ip_low ].A0;
    A0_21 = HQlattice[iT_high][ip_low ].A0;
    A0_12 = HQlattice[iT_low ][ip_high].A0;
    A0_22 = HQlattice[iT_high][ip_high].A0;

    B0_11 = HQlattice[iT_low ][ip_low ].B0;
    B0_21 = HQlattice[iT_high][ip_low ].B0;
    B0_12 = HQlattice[iT_low ][ip_high].B0;
    B0_22 = HQlattice[iT_high][ip_high].B0;

    B1_11 = HQlattice[iT_low ][ip_low ].B1;
    B1_21 = HQlattice[iT_high][ip_low ].B1;
    B1_12 = HQlattice[iT_low ][ip_high].B1;
    B1_22 = HQlattice[iT_high][ip_high].B1;

  //And now, the interpolation:
 
  A0  = A0_11*(T_high - T)*(p_high - pHQ) + A0_21*(T - T_low)*(p_high - pHQ);
  A0 += A0_12*(T_high - T)*(pHQ - p_low ) + A0_22*(T - T_low)*(pHQ - p_low );
  A0 /= dT*dp;
 
  B0  = B0_11*(T_high - T)*(p_high - pHQ) + B0_21*(T - T_low)*(p_high - pHQ);
  B0 += B0_12*(T_high - T)*(pHQ - p_low ) + B0_22*(T - T_low)*(pHQ - p_low );
  B0 /= dT*dp;
 
  B1  = B1_11*(T_high - T)*(p_high - pHQ) + B1_21*(T - T_low)*(p_high - pHQ);
  B1 += B1_12*(T_high - T)*(pHQ - p_low ) + B1_22*(T - T_low)*(pHQ - p_low );
  B1 /= dT*dp;

  return; 
}

void HQRates::readCoalesenceProb(
      string coal_prob_path)
{

    string filename;
    filename = coal_prob_path + "/Coalesence_probability.dat";
    ifstream file;
    file.open(filename.c_str(), ios::in);
    if(!file)
    {
        cerr << "[HQRates]: ERROR: Unable to open file " << filename << endl;
        exit(1);
    }

    double dummy, prob;
    for (int jj = 0; jj < 81 ; jj++)
    {
      file >> dummy >> prob;
      if (jj == 0) Coal_normalization = 1./prob;
      Coalesence_Prob_lattice[jj] = prob*Coal_normalization;
    }

    file.close();

    string filename_meson;
    filename_meson = coal_prob_path + "/Coalesence_Meson_probability.dat";
    ifstream file_meson;
    file_meson.open(filename_meson.c_str(), ios::in);
    if(!file_meson)
    {
        cerr << "[HQRates]: ERROR: Unable to open file " << filename_meson << endl;
        exit(1);
    }

    for (int jj = 0; jj < 81 ; jj++)
    {
      file_meson >> dummy >> prob;
      Coalesence_Meson_Prob_lattice[jj] = prob*Coal_normalization;
    }

    file_meson.close();
}

double HQRates::get_coalesence_prob(double pHQ)
{
    int ip_low, ip_high;
    double p_low, p_high, dp;
    double probability = 0.;
    if (pHQ > 20.) return probability; 

    ip_low  = pHQ/0.25;
    ip_high = ip_low + 1;
    p_low  = (double)ip_low * 0.25;
    p_high = (double)ip_high* 0.25;
    dp = 0.25;

    double prob_low  = Coalesence_Prob_lattice[ip_low];
    double prob_high = Coalesence_Prob_lattice[ip_high];

    //And now, the interpolation:
    probability = prob_low*(p_high - pHQ) + prob_high*(pHQ - p_low);
    probability /= dp;

    return probability; 
}

double HQRates::get_coalesence_meson_prob(double pHQ)
{
    int ip_low, ip_high;
    double p_low, p_high, dp;
    double probability = 0.;
    if (pHQ > 20.) return probability; 

    ip_low  = pHQ/0.25;
    ip_high = ip_low + 1;
    p_low  = (double)ip_low * 0.25;
    p_high = (double)ip_high* 0.25;
    dp = 0.25;

    double prob_low  = Coalesence_Meson_Prob_lattice[ip_low];
    double prob_high = Coalesence_Meson_Prob_lattice[ip_high];

    //And now, the interpolation:
    probability = prob_low*(p_high - pHQ) + prob_high*(pHQ - p_low);
    probability /= dp;

    return probability; 
}
