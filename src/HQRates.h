#ifndef HQRates_h
#define HQRates_h

#include<vector>
#include<cstring>
#include "HQratecell.h"
#include "Basics.h"

using namespace std;

class HQRates
{
    private:
        HQratecell** HQlattice;
	double* Coalesence_Prob_lattice; 
	double* Coalesence_Meson_Prob_lattice; 

    public:
        HQRates();   //constructor
        ~HQRates();  //destructor
	double Coal_normalization;
        void readHQRates(int scat_rad_flag, string rate_path);
        void getTandPdependent_coeff(double T, double pHQ, double& A0, double& B0, double& B1);
        void readCoalesenceProb(string coal_prob_path);
        double get_coalesence_prob(double pHQ);
        double get_coalesence_meson_prob(double pHQ);
	double get_Coal_normalization(){return Coal_normalization;};
};

#endif
