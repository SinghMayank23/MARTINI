//This contains many stand-alone routines and should be renamed!

// Given a 4-vector p and a 4-velocity u, performs the Lorentz boost to determine 
// p'. It assumes u_mu*u^mu = 1.
#include <math.h>
#include <iostream>
#define N_SAMPLES 100000
//#define N_SAMPLES 100000
//#define E_JPSI 0.33
//#define E_CCB_BOUND 0.88
//#define E_Y -0.05
//#define E_BBB_BOUND 2.2
//#define E_BBC 0.18
//#define E_BBC_BOUND 1.52
#define E_JPSI 3.4
#define E_CCB_BOUND 3.78
#define E_Y 9.5
#define E_BBB_BOUND 11.7
#define E_BBC 6.3
#define E_BBC_BOUND 7.6
#define DMASS 1.87
#define BMASS 5.28
#define MJPsi 3.1
#define MExcitedJPsi 3.5
#define MUpsilon 9.46
#define MExcitedUpsilon 10.0
#define MB_c 6.1
#define MExcitedB_c 7.2
#define EPS_C 0.15
#define CHARM_MASS 1.27
//#define CHARM_MASS 1.45
#define BOTTOM_MASS 4.2

using namespace std;

void LorentzBooster(double *u, double *p, double *pnew){
  
  double pdotS = 0.0 ;
  for(int i=1; i<4 ; i++) pdotS +=  p[i] * u[i] ;

  for (int i=1; i< 4; i++) 
    pnew[i] = p[i] + u[i] * ( pdotS/( u[0] + 1.) - p[0] );  
  pnew[0] = (p[0] * u[0] - pdotS ) ;
  //cout << "M^2 at the beginning = " << p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3] << endl;
  //cout << "u^2 = " << u[0]*u[0]-u[1]*u[1]-u[2]*u[2]-u[3]*u[3] << endl;
  //cout << "M^2 at the end = " << pnew[0]*pnew[0]-pnew[1]*pnew[1]-pnew[2]*pnew[2]-pnew[3]*pnew[3] << endl;
  return ;

  //double u2 = sqrt(u[0]*u[0]-1.);
  //cout << "u2 = " << u2 << endl;
  //double udotp = u[1]*p[1]+u[2]*p[2]+u[3]*p[3];
  // 
  //if(u2 > 0.){
  //  pnew[0] = u[0]*p[0]-udotp;
  //  pnew[1] = p[1]-u[1]*p[0]+(u[0]-1.)*u[1]*udotp/u2;
  //  pnew[2] = p[2]-u[2]*p[0]+(u[0]-1.)*u[2]*udotp/u2;
  //  pnew[3] = p[3]-u[3]*p[0]+(u[0]-1.)*u[3]*udotp/u2;
  //}
  //else{
  //  pnew[0] = p[0];
  //  pnew[1] = p[1];
  //  pnew[2] = p[2];
  //  pnew[3] = p[3];
  //}
}  
  
//Determines the Cornell Potential, a phenomenological potential that describes well the 
//quarkonia spectra, which however also is in agreement with lattice results and well-motivated 
//physically:

//double CornellPotential(double r) {

////In GeV:
//if(r < 0.2) r = 0.2 ;
//double potential = -0.0985/r+0.812*r ;
////A cutoff for large separations:
//if(potential > 0.88){
//  potential = 0.88 ; }
//
//return potential ;
//}

//Here, the Cornell potential has been set to zero at r=0.4 fm, and increases continuously from 
//there where V_C(r=0.4 fm) = 0. The purpose is to make this exactly the color evaporation model 
//for all quarks below tau_0 = 0.4 fm:

//double CornellPotential(double r) {
//
// //In GeV:
//if(r < 0.4) r = 0.4 ;
//double potential = -0.0985/r+0.812*r
//  -(-0.0985/0.4+0.812*0.4);
////A cutoff for large separations:
//if(potential > 0.88){
//  potential = 0.88 ; }
//
//return potential ;
//}
