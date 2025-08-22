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

    public:
        HQRates();   //constructor
        ~HQRates();  //destructor
        void readHQRates(int scat_rad_flag, string rate_path);
        void getTandPdependent_coeff(double T, double pHQ, double& A0, double& B0, double& B1);
        
};

#endif
