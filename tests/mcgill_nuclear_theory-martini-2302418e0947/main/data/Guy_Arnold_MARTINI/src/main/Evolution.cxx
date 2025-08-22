#include "Evolution.h"

using namespace std;

// This is a common variable used to keep the structure 
// of the program intact. I know. I shouldn't do that.

double COM_MD2_OVER_T2;

Evolution::Evolution(const FileNames* fileList, const Constants* constants, double idt, double idE, double iT, const double ialphas, int iNf, const double iinitialE, const double iwidth, const double imaxEnergy, const int iBetheHeitler, const int iBDMPS, const int idoRad, const int idoCol, const int inewRadGamma, const int size, const int rank)
{

  dt = idt;
  dE = idE;
  alphas = ialphas;
  Nf = iNf;
  T = iT;
  initialE = iinitialE;
  width = iwidth;
  maxEnergy = imaxEnergy; 
  BetheHeitler = iBetheHeitler;
  BDMPS = iBDMPS;
  doRad = idoRad;
  doCol = idoCol;
  newRadGamma = inewRadGamma;

  counter = 0;

  evolutionsFileNameList=fileList;
  
  PI = constants->PI;
  hbarc = constants->hbarc;
  EulerGamma = constants->EulerGamma;
  /// cb = -EulerGamma + Zeta'(2)/Zeta(2) + Log(2)
  cb = constants->cb; 
  /// cf = -EulerGamma + Zeta'(2)/Zeta(2)
  cf = constants->cf;
  cs = constants->cs;
   
  Nc = constants->Nc;


  LogEmax=constants->LogEmax;
  LogEmin=constants->LogEmin;
  LogStepE=constants->LogStepE;
  LogOmegaMin=constants->LogOmegaMin;
 LogStepOmega=constants->LogStepOmega;
  Tfile=constants->Tfile;
  stepT=constants->stepT;


  if(rank==0)
    {
      ofstream fout(evolutionsFileNameList->GetFN(1).c_str(),ios::app); 
      fout << " " << endl;
      fout << "[Evolution::Evolution]: " << endl;
      fout << "                Temperature = " << T << " GeV" << endl;
      fout << "                alpha_s = " << alphas << endl;
      fout << "                dt = " << dt << " GeV^(-1)" << endl;
      fout << "                dE = " << dE << " GeV" << endl;
      fout << "                Nf = " << Nf << endl;
      fout << "                initial jet energy = " << initialE << endl;
      fout << "                initial jet width = " << width << endl;
      fout << "                maximal energy = " << maxEnergy << " GeV" << endl;
      fout << "                Bethe-Heitler limit = " << BetheHeitler << endl;
      fout << "                BDMPS approximation = " << BDMPS << endl;
      fout << endl;
      fout.close();
    }

  if(doRad==1||newRadGamma==2) // prepare for radiative loss
    {
      if (newRadGamma==2) newRadGamma=1; // if the button to create new radiative transitiom rate list was hit 
                                         // - prepare to do it
      prep_equipment ( &dat , &Gam , &dGamma , &dGamma_gqq , &dGamma_ggg , 
						&dGamma_em );
    }
}

FILE *openfile ( char filename[] , char mode[] )
     /*  Gets a file's name and opens it in the mode given. */
{
  FILE *result;
  long  i;

  i = 0;
  while ( ( ( filename[i] = getchar() ) == '\n' ) || filename[i] == ' ' );
  while ( ( filename[++i] = getchar () ) != '\n' )
    if ( filename[i] == ' ' ) i--; /* read to end of line, throwing 
				      out empty spaces. */
  filename[i] = '\0';
  result = fopen ( filename , mode );

  return ( result );
}


//-----------------------------------------------------------------------------------//
// begin part for radiative energy loss - taken from Qin ----------------------------//
//-----------------------------------------------------------------------------------//

// This includes Peter's interpolation 
double Evolution::K(double z, int just_K0)
/* Sangyong's mod to include Peter Arnold's interpolation */
     /* BK(z) = K_0(z) + ln(z/2) + Gamma_E
     This is combined in the following way:
     N
     [   (1 + a^2/(s^2 - 1) ) BK(z, just_K0)
         - ( (a^2/s^2)/(s^2 - 1) ) BK(z s, just_K0) ]
     where a^2/s^2 = 0.3
     and s^2 = T^2/m_D^2
     The parameters are:
     N = 0.849346
     a^2 = 0.132827(T/m_D)^2
     s^2 = 0.82971(T/m_D)^2
     For the sake of keeping the structure of the program intact,
     I am going to use one global constant:
     COM_MD2_OVER_T2 (Debye mass squared over T2).

     */
{
double f;

 double a2, s2, s, normf;

 normf = 0.849346;
 s2 = 0.82971/COM_MD2_OVER_T2;
 s = sqrt(s2);
 a2 = 0.132827/COM_MD2_OVER_T2;

 f = (1.0 + a2/(s2-1.0))*BK(z, just_K0);
 f -= ((a2/s2)/(s2-1.0))*BK(s*z, just_K0);
 f *= normf; // overall normalization;

/*
 f = BK(z, just_K0);
*/

 return f;
}// end of K

// Sangyong's mod:
// This is the original K renamed to BK
// There are places in the code (keyword: Green) where K is changed to BK
// But that part is for tau = <1/\delta E> calculation
// and hence does not really concern us

double Evolution::BK(double z, int just_K0)
     /* Determines K(z) = K_0(z) + ln(z/2) + Gamma_E either by
	power series or by asymptotic series.  10 digit accurate */
     /* I = I_0(z) - 1.  (or I_0(z) to get just_K0).*/
     /* This is not very elegant.  It uses the power series except 
	for large z, where it uses an asymptotic series.  The power
	series suffers from big, cancelling terms, hence only 10 
	digit accuracy.  Is there a better way using recurrence relations? */
{
  double sofar , I , coeff , coeff2 , logplus;
  int    j;
  
  if ( z < 0 ) return ( BK ( -z , just_K0 ) );
  if ( z == 0 ) return ( 0 ); // Safety first.
  if ( z < 11 ) // Use power series 
    {
      logplus = log( 0.5 * z ) + EulerGamma;
      z *= 0.25 * z; // z*z/4 is what series is in. 
      I = ( just_K0 ? 1 : 0 );
      coeff2 = 1; // keeps track of 1+1/2+1/3+... which is added to logplus.
      j = 1;
      coeff = z;  // will become z^j / (j!)^2 
      sofar = 0;
      while ( coeff > 1.0e-12 )
	{
	  I += coeff;
	  sofar += coeff * coeff2;
	  j++;
	  coeff *= z / ( j * j ); /* z^j / (j!)^2 */
	  coeff2 += 1.0l/j;       /* 1 + 1/2 + 1/3 + .. 1/j */
	}
      sofar -= I * logplus;
      
      return ( sofar );
    }
  // That failing, use asymptotic series: 
  logplus = ( just_K0 ? 0 : log ( 0.5 * z ) + EulerGamma );
  I = sqrt(PI/(2*z)) * exp(-z); // (leading order asymptotic of K0) 
  z = 1.0l/z;
  sofar = ( ( ( -0.0732421875 * z + 0.0703125 ) * z - 0.125 ) * z + 1 );
  // the numbers are 1, -1/8, 1*9/8*8*2!, -1*9*25/8*8*8*3! 
  sofar = I * sofar + logplus;
  return ( sofar );
}

double Evolution::K_BDMPS ( double z )
     /* Evaluates the BDMPS approximation version of the above,
	that is, 0.5 * ( 1 - z K_1(z) ), K_1 the modified Bessel
	function. */
{
  double sofar , I , coeff1 , psi , the_log;
  int    j;

  if ( z < 0 ) return ( K_BDMPS ( -z ) );
  if ( z == 0 ) return ( 0 );
  if ( z < 9 )
    {
      // Note, I will really be (z/2) I_0.

      the_log = log ( 0.5 * z );
      z *= 0.25 * z; // Series is in z^2 / 4. 
      coeff1 = z;
      I = 0; sofar = 0;
      psi = - 2 * EulerGamma + 1; // psi(1) + psi(2) 
      j = 0;
      while ( coeff1 > 1.0e-13 )
	{
	  I += coeff1;
	  sofar += coeff1 * psi;
	  j++;
	  coeff1 *= z / ( j * (j+1) );
	  psi += 1.0l/j + 1.0l/(j+1);
	}
      sofar = 0.5 * sofar - the_log * I;
      // This is (1/2) (1 - z K_1(z))
    }
  else
    {
      sofar = sqrt ( 0.125 * z * PI ) * exp(-z);
      I = 1;
      coeff1 = 1;
      z = 0.125 / z;
      for ( j = 1 ; j < 12 ; j++ )
	{
	  coeff1 *= (4 - (2*j-1)*(2*j-1)) * ( z / j );
	  I += coeff1;
	}
      sofar = 0.5 - I * sofar;
    }
  return ( sofar );
}

/** ------------------------------------------------------------------------
    Now, function to evaluate derivative, for solving differential
    equation (for determining splitting kernel) by method of Runge-Kutta.  
    ------------------------------------------------------------------------
    
    f[0] = Re f
    f[1] = Im f
    f[2] = Re f'
    f[3] = Im f'
    Therefore (fpr =f'):
    fpr[0] = f[2]
    fpr[1] = f[3]
    fpr[2] = -3/b f[2] + A f[0] + DK f[1]
    fpr[3] = -3/b f[3] + A f[1] - DK f[0]                                     */

#define get_fpr(f,fpr,A,DK,b) {double b_coeff;\
        b_coeff=3/(b); fpr[0]=f[2];fpr[1]=f[3];\
        fpr[2]=-b_coeff*f[2]+A*f[0]+DK*f[1];\
        fpr[3]=-b_coeff*f[3]+A*f[1]-DK*f[0];}

/* -----------------------------------------------------------  */

void Evolution::RK_step ( double f[4] , double A , double D , double *b , double *k ,  
			  double C[3] , double r[3] ) 
     /* Performs a 4th order Runge-Kutta step.  
  	Internally determines appropriate step length. 
  	Note that I am stepping from larger to smaller b.   */
     /* C are the Casimirs for the three terms in the collision piece, 
  	r are the rescalings of b for those three terms.   */
{ 
  double fpr1[4] , fpr2[4] , fpr3[4] , fpr4[4] , ftmp[4]; 
  double step; 
  int    i; 
  
  /* First, determine how large a step is permitted.  */
  step = 3.0l / *b; 
  step = 0.03l / sqrt ( step * step + sqrt ( A*A + D*D* *k * *k ) ); 
  /*  The quantity in the square root is an estimate of the square of scale 
      of variation of the function f.  (scale)/30 is good enough with 
      Runge-Kutta to give 10^{-8} errors.  Tighter doesn't help (checked)  */
  
  /* Recall that RK update is: 
     Get fpr1 = f'[b,f[b]] 
     fpr2 = f'[b+step/2,f[b]+step fpr1 /2] 
     fpr3 = f'[b+step/2,f[b]+step fpr2 /2] 
     fpr4 = f'[b+step,f[b]+step fpr3] 
     fprime = (1/6) [fpr1 + fpr4 + 2 fpr2 + 2 fpr3] 
     The f[b]+... is an estimate of the value at b+step/2 or b+step 
     using the fprime's just derived.   */
  
  get_fpr ( f , fpr1 , A , D* *k , *b ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - 0.5 * fpr1[i] * step; 
  /* fpr1 = slope at b; ftmp = first midpt extrapolation.  */
  *b -= 0.5 * step; 
  if ( BDMPS )
    *k = C[0] * K_BDMPS (*b * r[0] ) + C[1] * K_BDMPS (*b * r[1] ) 
      + C[2] * K_BDMPS (*b * r[2] ); 
  else
    *k = C[0] * K (*b * r[0] , 0 ) + C[1] * K (*b * r[1] , 0 ) 
      + C[2] * K (*b * r[2] , 0 ); 
  get_fpr ( ftmp , fpr2 , A , D* *k , *b ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - 0.5 * fpr2[i] * step; 
  /* fpr2 = first guess slope at b-step/2; ftmp=second mid extrap.  */
  get_fpr ( ftmp , fpr3 , A , D* *k , *b ); 
  /* fpr3 = second guess slope at b-step/2.  */
  *b -= 0.5 * step; 
  if ( BDMPS )
    *k = C[0] * K_BDMPS (*b * r[0] ) + C[1] * K_BDMPS (*b * r[1] ) 
      + C[2] * K_BDMPS (*b * r[2] ); 
  else
    *k = C[0] * K (*b * r[0] , 0 ) + C[1] * K (*b * r[1] , 0 ) 
      + C[2] * K (*b * r[2] , 0 ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - fpr3[i] * step; 
  get_fpr ( ftmp , fpr4 , A , D* *k , *b ); 
  /* ftmp=endpt extrap and fpr4 is guess slope there.  */
  for ( i = 0 ; i < 4 ; i++ ) 
    f[i] -= step * ( fpr1[i] + 2*fpr2[i] + 2*fpr3[i] + fpr4[i] ) * 0.166666666666666666666667; 
} 

double Evolution::find_f_atzero ( double A , double D , double C[3] , double r[3] ,
				  double *one_over_dE ) 
  /* solves the differential equation with exponentially shrinking 
     large value boundary conditions to find the finite piece 
     orthogonal to the divergent piece at the origin. 
     Here C are the Casimirs and r the rescalings for the 3 terms 
     in the collision piece.   */
  /* This version also computes < 1 / delta E > by performing the
     integral of f over the Green function of 1 / delta E. */
{ 
  double f[4] , b , b_start , KD , b_min , result; 
  double bess_coeff[4] , j , jpr , y , ypr , k , shift; 
  double dE_inv_int[2] , dE_inv_increment[2] , delta_b , root_A , kk; 
  
  /* First task:  what b is big enough?  Based on rate of exponential  
     falloff in the tail--ignore 1/b term and b dependence of K.   */
  b = 2;  
  dE_inv_int[0] = 0;
  dE_inv_int[1] = 0;
  root_A = sqrt ( A );
  do 
    { 
      b_start = b; 
      if ( BDMPS )
	KD = D*( C[0] * K_BDMPS (b * r[0] ) + C[1] * K_BDMPS (b * r[1] ) +  
		 C[2] * K_BDMPS (b * r[2] ) ); 
      else
	KD = D*( C[0] * K (b * r[0] , 0 ) + C[1] * K (b * r[1] , 0 ) +  
		 C[2] * K (b * r[2] , 0 ) ); 
      b = 28 / sqrt ( A + sqrt ( A*A + KD * KD ) ); 
      /* denominator = real part of negative eigenvalue /sqrt(2)  */
      b_start /= b; 
    } 
  while ( b_start < 0.95 || b_start > 1.05 ); 
  f[0] = 1; 
  f[1] = 0; 
  f[2] = -sqrt ( 0.5 * ( sqrt ( A*A + KD*KD ) + A ) ); 
  f[3] =  sqrt ( 0.5 * ( sqrt ( A*A + KD*KD ) - A ) ); 
  /* Approximate shrinking solution--assumed 
     1) you can neglect f'/b term 
     2) you can take K constant.  Both valid far enough out in tail.  */
  /* Now find how small b must go before I accurately see b~=0 behavior  */
  b_min = 0.025 / sqrt ( 1 + sqrt ( A*A+D*D ) ); 
  if ( BDMPS )
    k = C[0] * K_BDMPS (b * r[0] ) + C[1] * K_BDMPS (b * r[1] ) 
      + C[2] * K_BDMPS (b * r[2] ); 
  else
    k = C[0] * K (b * r[0] , 0 ) + C[1] * K (b * r[1] , 0 ) 
      + C[2] * K (b * r[2] , 0 ); 
  kk = BK ( b*root_A , 1 ); /* Green function of 1 / delta E */
  /* Now run differential equation from b to b_min  */
  while ( b > b_min ) 
    { 
      delta_b = b;
      dE_inv_increment[0] = kk * ( f[0] + 0.5 * b * f[2] );
      dE_inv_increment[1] = kk * ( f[1] + 0.5 * b * f[3] );
      RK_step ( f , A , D , &b , &k , C , r ); 
      kk = BK ( b*root_A , 1 ); /* Green function of 1 / delta E */
      dE_inv_increment[0] += kk * ( f[0] + 0.5 * b * f[2] );
      dE_inv_increment[1] += kk * ( f[1] + 0.5 * b * f[3] );
      delta_b -= b;
      dE_inv_int[0] += dE_inv_increment[0] * delta_b;
      dE_inv_int[1] += dE_inv_increment[1] * delta_b;
      
      if ( f[0]*f[0] + f[1]*f[1] > 100 ) /* Apply harmless rescaling  */
  	{ 
	  f[0] *= 0.125; f[1] *= 0.125; 
	  f[2] *= 0.125; f[3] *= 0.125; 
	  dE_inv_int[0] *= 0.125; dE_inv_int[1] *= 0.125;
	} 
    } 
  /* Find projection of solution onto divergent and finite Bessel 
     type (D=0) solutions.  Then find the finite part, orthogonal to the 
     divergent one.  bess_coeff[0] and [1] are the Re, Im parts of 
     divergent, and bess_coeff[2] and [3] of finite, solutions, 
     normalized to to 1/b^2 and 1.   */
  
  /* fourth-order finite and divergent solutions with derivatives
     Note that "j" = 2 j_1/b and "y" = y_1/b plus some (A dependent) 
     multiple of j to make the logs simpler.  This doesn't matter 
     since I want "j" in the direction orthogonal to where "y" != 0.  */
  j = 1 + 0.125 * A*b*b * ( 1 + .04166666666*A*b*b ); 
  jpr = 0.25 * A * b * (1+.083333333333*A*b*b); 
  y = 1.0/(b*b) + 0.5 * A * ( log(b) - 0.5 + b*b*A*0.125* 
  			      (log(b)-1.25) ); 
  ypr = -2.0/(b*b*b) + 0.5 * A / b 
    + 0.125*A*A*b*(log(b)-0.75); 
  
  bess_coeff[0] = b * b * b * ( jpr * f[0] - j * f[2] ); 
  bess_coeff[1] = b * b * b * ( jpr * f[1] - j * f[3] ); 
  bess_coeff[2] = b * b * b * ( -ypr* f[0] + y * f[2] ); 
  bess_coeff[3] = b * b * b * ( -ypr* f[1] + y * f[3] ); 
  /* The b^3 is a Jacobian for the solutions . . .   */
  result = (bess_coeff[0]*bess_coeff[3]-bess_coeff[1]*bess_coeff[2])/ 
    (bess_coeff[0]*bess_coeff[0]+bess_coeff[1]*bess_coeff[1]); 
  shift = 0.0625 * D*b*b* ( C[0] * r[0]*r[0]*(1.5-EulerGamma-log(b*r[0]/2) ) 
  			    +C[1] * r[1]*r[1]*(1.5-EulerGamma-log(b*r[1]/2) ) 
  			    +C[2] * r[2]*r[2]*(1.5-EulerGamma-log(b*r[2]/2) ) ); 
  /* Correction due to D in interval from b to 0, leading 
     analytic form.  */
  result += shift; 
  
  *one_over_dE = (bess_coeff[0]*dE_inv_int[1]-bess_coeff[1]*dE_inv_int[0])/
    (bess_coeff[0]*bess_coeff[0]+bess_coeff[1]*bess_coeff[1]);
  *one_over_dE /= result;
  /* Ratio of (1 / delta E) operator acting on Re part / Re part */
  
  return ( result ); 
} 

double Evolution::solve_int_eqn ( double p , double k ,  
				  double df , double da , double cf , double ca , 
				  int  n_quark_flavors , int p_is_gluon , 
				  int all_are_gluons , int emit_photon , 
				  double *one_over_deltaE ) 
  /* Solves the integral equations to give the production  
     rate at p,k.  Contains all the annoying overall factors,
     except a factor of Nf in g->qq.  */
{ 
  double A , B , D , sofar , C[3] , r[3]; 
  double md2_over_g2T2 , kappa , tmp; 
  
  if ( k*k < 1.0e-9 || (k-p)*(k-p) < 1.0e-9 ) 
    { 
      /* Do average on either side of p, to avoid p=0 complications 
  	 in the following procedure.  */
      sofar = (solve_int_eqn( p , k-.001 , df , da , cf , ca , 
			      n_quark_flavors , p_is_gluon , all_are_gluons
			      , emit_photon, one_over_deltaE )
	       +solve_int_eqn ( p , k+.001 , df , da , cf , ca , 
				n_quark_flavors,p_is_gluon , all_are_gluons
				, emit_photon, one_over_deltaE )
		) * 0.5;
      return ( sofar );
    }
  md2_over_g2T2 = ( ca + n_quark_flavors * cf * df / da ) / 3.0l;
  kappa = 0.25 * cf / md2_over_g2T2;
  if ( p_is_gluon ) 
    { /* Then gluon goes either to 2 quarks or 2 gluons. */
      A = -1 / ( 4 * p );
      if ( all_are_gluons )
	A += 1 / ( 4.0 * k ) + 1 / ( 4.0 * (p-k) );
      else
	A += kappa * ( 1 / ( 2.0 * k ) + 1 / ( 2.0 * (p-k) ) );
    }
  else
    { /* p and (p-k) are quarks and k is a gluon. */
      A = kappa * ( 1 / ( 2.0 * ( p - k ) ) - 1 / ( 2 * p ) );
      if ( ! emit_photon )  /* Photon is approximately massless */
	A += 1 / ( 4 * k ); /* Gluon mass contribution. */
    }
  A *= md2_over_g2T2;
  B = md2_over_g2T2 * p / ( 2 * k * (p-k) );
  D = 1.0 / ( 2 * PI );
  
  /* Respectively, const and p^2 coefficients of delta E, and
     coefficient of collision term. */

  A /= B;
  D /= B;
  if ( BDMPS ) 
  A = 0; /* Because BDMPS forgot about thermal masses. */

  r[0] = 1;
  r[1] = ( ( k > 0 ) ? k / p : -k / p );
  r[2] = ( (p-k) > 0 ? (p-k)/p : (k-p)/p );
  if ( p_is_gluon )
    {
      if ( all_are_gluons )
	{
	  C[0] = ca / 2;
	  C[1] = ca / 2;
	  C[2] = ca / 2;
	}
      else
	{
	  C[0] = cf - ca / 2;
	  C[1] = ca / 2;
	  C[2] = ca / 2;
	}
    }
  else if ( emit_photon )
    {
      C[0] = 0;
      C[1] = cf;
      C[2] = 0;
    }
  else
    {
      C[0] = ca/2;
      C[1] = cf - ca/2;
      C[2] = ca/2;
    }

  /* Solving everything inside the p_parallel integration. */
  if ( BetheHeitler )
    {
      /* In this case we want to work perturbatively in D.  Do this
	 by determining things at D= (correct D)/100 and (correct D)/200
	 and use these results to extrapolate to zero D. */
      /* This is the "bad way."  I should fix this some day. */

      sofar = 400  * find_f_atzero ( A , D/200 , C , r , 
				     one_over_deltaE );
      tmp = *one_over_deltaE;
      sofar -= 100 * find_f_atzero ( A , D/100 , C , r , 
				     one_over_deltaE );
      *one_over_deltaE = 2*tmp - *one_over_deltaE;
      sofar *= 4 / ( PI * B );
    }
  else
    {
      sofar = 4 * find_f_atzero ( A , D , C , r , 
				  one_over_deltaE ) / ( PI * B );
    }  
  *one_over_deltaE /= B;
  sofar *= p*p / (16 * PI ) * md2_over_g2T2 * md2_over_g2T2;
  if ( p_is_gluon )
    {
      if ( all_are_gluons )
	{
	  sofar *= ca * ( p*p*p*p + k*k*k*k + (p-k)*(p-k)*(p-k)*(p-k) )
	    / ( p*p*p*k*k*k*(p-k)*(p-k)*(p-k) )
	    * BoseStim(k) * BoseStim((p-k));
	}
      else
	{
	  sofar *= cf * ( k*k + (p-k)*(p-k) ) / ( k*k*(p-k)*(p-k)*p*p*p )
	    * PauliBlock(k) * PauliBlock((p-k));
	}
    }
  else
    {
      sofar *= ( p*p + (p-k)*(p-k) ) / ( p*p*(p-k)*(p-k)*k*k*k )
	* PauliBlock((p-k));
      if ( emit_photon ) 
	{
//////////////////////////////////////////////////////////////////////////
//	The old program double counts this factor
//	(4n/2+(n+1)/2)/(9n)=6/27 for n=3, n/2=1
//////////////////////////////////////////////////////////////////////////
/*	  sofar *= ( 4*(n_quark_flavors/2) + ( (n_quark_flavors+1)/2 ) )
	    /((double) (9*n_quark_flavors));	*/
//////////////////////////////////////////////////////////////////////////
	/* Average charge squared, assuming number of up-type is equal or one
	     less than number of down type fermions. */
	  if ( k < 0 )
	    sofar = 0; /* Because there are no photons to "pick up". */
	}
      else
	{
	  sofar *= cf * BoseStim(k);
	}
    }
  return ( sofar );
  /* Rate per unit k, except for g^4 T factor (and Nf for g->qqbar)! */
}


/** ------------------------------------------------------------------------
    Now, equipment for building a table for interpolation, and for
    storing, reading, and using that table. 
    ------------------------------------------------------------------------

    NORMALIZATION:
    
    In
     build_table
     write_table
     read_table
     use_table
     prepare_table
     prep_dGamma
     and in the stored tables
    
    The normalization is, that the dGamma(p,k)/dk dx given, is
    in units of 1/g^4 except for the photon one, which is in units
    of 1/g^2 e^2.  That means, that to convert dGamma(p,k)/dk dx
    to the rate of emission per momentum range $k$ and distance or
    time $x$, one should multiply by g^4 or g^2 e^2.
    In other words, the dGamma, integrated \int dk/T, gives the
    likelihood to make an emission in a length 1/(g^4 T) of plasma.

    In
     find_dP
     evolve_one_step
   
    the normalization of dx or dt is, that it is 1/g^4 times an inverse
    energy (in GeV).

    It is in evolve_in_medium that the conversion to Fermis is performed,
    meaning the inclusion of the 1/g^4 and of the relation between a
    Fermi and a GeV^-1.  
*/

void Evolution::build_table ( Gamma_info * dat , dGammas *Gam )
     /* Determines the ranges of p,k and spacing dp to consider,
	and loads up an array of values of dGamma/dkdx needed for 
	updating. */
{
  double dtau , k , p , rate, md2_over_g2T2;
  int    i_p , i_k , i_g;
  char   c;
  
  dat->ddf = Nc;
  dat->dda = Nc * Nc - 1;
  dat->dcf = dat->dda / ( 2 * dat->ddf );
  dat->dca = Nc;

 printf("Nf = %d\n", Nf);
 printf("alphas = %e\n", alphas);
 printf("PI = %e\n", PI);

  md2_over_g2T2 = ( dat->dca + Nf * (dat->dcf) * (dat->ddf) /(dat->dda) ) / 3.0l;
  printf("md2/g2T2 = %e\n", md2_over_g2T2);
  
  COM_MD2_OVER_T2 = md2_over_g2T2*(alphas*4.0*PI);
  printf("md2/T2 = %e\n", COM_MD2_OVER_T2);

  dat->include_gluons = 1;  /* Gluon inclusion will be the default, since we are including photons as well! */

  for ( i_p = 0 ; i_p < NP ; i_p++ )
    for ( i_k = 0 ; i_k < NK ; i_k++ )
      {
	Gam->dGamma[i_p][i_k] = 0; 
	Gam->dGamma_gqq[i_p][i_k] = 0;
	Gam->dGamma_ggg[i_p][i_k] = 0; /* Initialize them all to zero. */
	Gam->dGamma_em[i_p][i_k] = 0;
	Gam->tau[i_p][i_k] = 0;
	Gam->tau_gqq[i_p][i_k] = 0;
	Gam->tau_ggg[i_p][i_k] = 0;
	Gam->tau_em[i_p][i_k] = 0;
      }
  for ( i_g = 0 ; i_g < 4 ; i_g++ )
    for ( i_p = 0 ; i_p < NP ; i_p++ )
      {
	if ( i_p == 0 )
	  p = 4.01;
	else
	  p = p * 1.04119; /* spaced so 6---1000 is 0--127 */
	printf ( "Starting p %e\n" , p );
	for ( i_k = 0 ; i_k < ( ( (i_g % 3) == 0 ) ? NK : NK - 160 ) ; i_k++ )
	  {
	    if ( i_k < 50 )        /* spaced by 0.2  from -12 to -2 */
	      k = -12 + i_k * 0.2;
	    else if ( i_k < 60 )   /* spaced by 0.1  from -2  to -1 */
	      k = -2 + (i_k-50) * 0.1;
	    else if ( i_k < 100 )  /* spaced by 0.05 from -1  to +1 */
	      k = -1 + (i_k-60) * 0.05;
	    else if ( i_k < 110 )  /* spaced by 0.1  from +1  to +2 */
	      k = 1 + (i_k-100) * 0.1;
	    else if ( i_k < 270 )  /* spaced complicated, +2 to p-2 */
	      {
		k = 0.1 * (i_k-190);
		k = 2 + (p-4) * ( -0.0003353501304664781l
				  + 1.000670700260932956l / (1+exp(-k)) );
	      }
	    else if ( i_k < 280 )  /* spaced by 0.1  from p-2 to p-1 */
	      k = p - 2 + 0.1 * (i_k-270);
	    else if ( i_k < 320 )  /* spaced by 0.05 from p-1 to p+1 */
	      k = p + 0.05 * (i_k - 300);
	    else if ( i_k < 330 )  /* spaced by 0.1  from p+1 to p+2 */
	      k = p + 0.1 * (i_k - 310);
	    else                   /* spaced by 0.2  from p+2 to p+12 */
	      k = p + 0.2 * (i_k - 320);
	    
	    if ( ABS(k) > 0.001 )
	      {
		rate = 
		  solve_int_eqn ( p , k , dat->ddf , dat->dda , dat->dcf , 
				  dat->dca , Nf , i_g%3 , (i_g == 2) , 
				  (i_g == 3) , &dtau );
		switch ( i_g )
		  {
		  case 0:
		    rate *=  k;
		    if ( k < 20 ) rate *= 1 - exp(-k);
		    if ( k > p - 20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma[i_p][i_k] = rate;
		    Gam->tau[i_p][i_k] = dtau;
		    break;
		  case 1:
		    rate *= p;
		    if ( k < 20 ) rate *= 1 + exp(-k);
		    if ( k > p - 20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma_gqq[i_p][i_k] = rate;
		    Gam->tau_gqq[i_p][i_k] = dtau;
		    break;
		  case 2:
		    rate *= k * (p-k) / p;
		    if ( k < 20 ) rate *= 1 - exp(-k);
		    if ( k > p - 20 ) rate *= 1 - exp(k-p);
		    Gam->dGamma_ggg[i_p][i_k] = rate;
		    Gam->tau_ggg[i_p][i_k] = dtau;
		    break;
		  case 3:
		    rate *= k;
		    if ( k < 0 ) rate = 0; /* No photon pickup allowed */
		    if ( k > p-20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma_em[i_p][i_k] = rate;
		    Gam->tau_em[i_p][i_k] = dtau;
		    break;
		  }
	      }
	  }
	Gam->tau[i_p][80] = 0; /* it is! */
	Gam->tau_gqq[i_p][80] = 0;
	Gam->tau_ggg[i_p][80] = 0;
	Gam->tau_em[i_p][80] = 0;
	switch ( i_g )
	  {
	  case (0):
	    Gam->dGamma[i_p][80] = 0.5 * ( Gam->dGamma[i_p][79] 
					   + Gam->dGamma[i_p][81] );
	    break;
	  case(1):
	    Gam->dGamma_gqq[i_p][80] = 0.5 * ( Gam->dGamma_gqq[i_p][79] 
					       + Gam->dGamma_gqq[i_p][81] );
	    break;
	  case(2):
	    Gam->dGamma_ggg[i_p][80] = 0.5 * ( Gam->dGamma_ggg[i_p][79] 
					       + Gam->dGamma_ggg[i_p][81] );
	  case(3):
	    Gam->dGamma_em[i_p][80] = 0.5 * ( Gam->dGamma_em[i_p][79] 
					      + Gam->dGamma_em[i_p][81] );
	  }
      }
}

void Evolution::write_table ( Gamma_info *dat , dGammas Gam )
     /* Writes the table, in binary, to a file. */
{
  char fname[100];
  FILE *wfile;

  printf ( "Enter the file name to write Gamma information into.\n" );
  wfile = openfile ( fname , "wb" );
  fwrite ( (char *) (&dat->ddf) , sizeof ( double ) , 1 , wfile );
  fwrite ( (char *) (&dat->dda) , sizeof ( double ) , 1 , wfile );
  fwrite ( (char *) (&dat->dcf) , sizeof ( double ) , 1 , wfile );
  fwrite ( (char *) (&dat->dca) , sizeof ( double ) , 1 , wfile );
  fwrite ( (char *) (&dat->Nc) , sizeof ( int ) , 1 , wfile );
  fwrite ( (char *) (&dat->Nf) , sizeof ( int ) , 1 , wfile );
  fwrite ( (char *) (&dat->BetheHeitler) , sizeof ( int ) , 1 , wfile );
  fwrite ( (char *) (&dat->BDMPS) , sizeof ( int ) , 1 , wfile );
  fwrite ( (char *) (&dat->include_gluons) , sizeof ( int ) , 1 , wfile );
  fwrite ( (char *) Gam.dGamma , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.tau , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.dGamma_gqq , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.tau_gqq , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.dGamma_ggg , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.tau_ggg , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.dGamma_em , sizeof ( double ) , NP * NK , wfile );
  fwrite ( (char *) Gam.tau_em , sizeof ( double ) , NP * NK , wfile );
  fclose ( wfile );
}


double Evolution::use_table ( double p , double k , double dGamma[NP][NK] , 
			      double tau[NP][NK] , double *inv_dE , int which_kind )
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
  *inv_dE = (1-a) * ( (1-b) * tau[n_p][n_k] + b * tau[n_p][n_k+1] )
    +         a * ( (1-b) * tau[n_p+1][n_k] + b * tau[n_p+1][n_k+1] );
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

void Evolution::prep_dGamma ( Gamma_info *dat , dGammas *gam , 
			      double T , double beta , double cos_phi , 
			      double *dGamma , double *dGamma_gqq , 
			      double *dGamma_ggg , double *dGamma_em )
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
	    use_table ( p , k , gam->dGamma , gam->tau , &junk , 0 );
	  dGamma_em[here] = jacobian *
	    use_table ( p , k , gam->dGamma_em , gam->tau_em , &junk , 0 );

	  tot_dGamma += dGamma[here];
	  if ( 4 * i_k - 2 * dat->n_kmin <= 3 + i_p + dat->n_pmin )
	    {
	      dGamma_gqq[here] = Nf * jacobian *
		use_table ( p , k , gam->dGamma_gqq , gam->tau , 
			    &junk , 1 );
	      dGamma_ggg[here] = jacobian *
		use_table ( p , k , gam->dGamma_ggg , gam->tau , &junk , 2 );
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

double Evolution::find_dP ( double *P , double *Pg , double *Pem , double *dGamma , 
			    double *dGamma_gqq , double *dGamma_ggg , 
			    double *dGamma_em , double *dP , 
			    double *dPg , double *dPem , Gamma_info dat )
     /* Evaluates dP/dt and computes, and returns, the "safety margin"
	for the evolution algorithm.  

	At this point, dt is in units of 1/g^4 E, with E in GeV.  */
{
  register int    i_p , i_k , i_ppr , i_kpr , n_kmin , n_p , n_pmin , inc_glue;
  double          tot_change , tot_prob , big_change , ratio_alphas;
  register double chg_here , scale , weight;

  n_kmin = dat.n_kmin;
  n_p = dat.n_p;
  inc_glue = dat.include_gluons;
  n_pmin = dat.n_pmin;
  ratio_alphas = dat.alpha / dat.alpha_s;
  scale = 2 * dat.dp;
  /* The rescaling of the step, to turn dGamma/dk into a finite difference
     in k.  Remember that the k values are spaced with twice the 
     interval of the p values. */
  for ( i_p = 0 ;  i_p < n_p ; i_p ++ )
    {
      dP[i_p] = 0;
      dPg[i_p] = 0;
      dPem[i_p] = 0;
    }
  for ( i_k = 0 ; i_k < dat.n_k ; i_k ++ )
    for ( i_p = 0 ; i_p < n_p ; i_p ++ )
      { /* Update probability due to emission by particle, momentum p. */
	i_ppr = i_p - 2*i_k + n_kmin; /* location of final momentum. */
	i_kpr = 2 * i_k - n_kmin - n_pmin; /* location of k momentum. */

	/* First, quark emitting a gluon. */
	chg_here = P[i_p] * dGamma[i_p + n_p * i_k] * scale;	
	dP[i_p] -= chg_here;
	if ( (i_ppr >= 0) & (i_ppr < n_p) )
	  dP[i_ppr] += chg_here; /* Keep track of scattered state. */
	if ( ( i_kpr >= 0 ) & ( i_kpr < n_p ) )
	  dPg[i_kpr] += 0.5*chg_here;
	if ( ( i_kpr+1 >= 0 ) & ( i_kpr+1 < n_p ) )
	  dPg[i_kpr+1] += 0.25*chg_here;
	if ( ( i_kpr-1 >= 0 ) & ( i_kpr-1 < n_p ) )
	  dPg[i_kpr-1] += 0.25*chg_here;
	/* Accounts for the emitted gluon from the quark. */

	/* Next, quark emitting a photon. */
	chg_here = P[i_p] * dGamma_em[i_p + n_p * i_k] * scale
	  * ratio_alphas;
	if ( (i_kpr >= 0) & (i_kpr < n_p) )
	  dPem[i_kpr] += 0.5*chg_here;
	if ( (i_kpr+1 >= 0) & (i_kpr+1 < n_p) )
	  dPem[i_kpr+1] += 0.25*chg_here;
	if ( (i_kpr-1 >= 0) & (i_kpr-1 < n_p) )
	  dPem[i_kpr-1] += 0.25*chg_here;
	/* I ONLY keep track of accumulated photon, ignoring 
	   reduction in quark number, which is small. */
	/* NOTE:  the factor of q^2 of the quark is NOT included here!
	   You have to remember to multiply by 4/9 (up) or 1/9 (down
	   or strange) at the end.  */

	/* Next, gluon breakup.  Assume k<p/2 by symmetry of final states */
	if ( i_ppr >= i_kpr - 1 )
	  {
	    /* First, get edge effects right! */
	    if ( i_ppr - i_kpr < 2 )
	      weight = 0.5 + 0.25 * ( i_ppr - i_kpr );
	    else
	      weight = 1;
	    
	    /* First, gluon to two quarks */
	    chg_here = Pg[i_p] * weight * 
	      dGamma_gqq[i_p + n_p * i_k] * scale;

	    dPg[i_p] -= chg_here;
	    if ( ( i_ppr >= 0 ) & ( i_ppr < n_p ) )
	      dP[i_ppr] += chg_here; /* Gain is to quarks. */
	    
	    if ( ( i_kpr >= 0 ) & ( i_kpr < n_p ) )
	      dP[i_kpr] += 0.5 * chg_here;
	    if ( ( i_kpr-1 >= 0 ) & ( i_kpr-1 < n_p ) )
	      dP[i_kpr-1] += 0.25 * chg_here;
	    if ( ( i_kpr+1 >= 0 ) & ( i_kpr+1 < n_p ) )
	      dP[i_kpr+1] += 0.25 * chg_here;
	    
	    /* Finally, gluon to gluons. */
	    chg_here = Pg[i_p] * weight * 
	      dGamma_ggg[i_p + n_p * i_k] * scale;

	    dPg[i_p] -= chg_here;
	    if ( ( i_ppr >= 0 ) & ( i_ppr < n_p ) )
	      dPg[i_ppr] += chg_here;
	    
	    if ( ( i_kpr >= 0 ) & ( i_kpr < n_p ) )
	      dPg[i_kpr] += 0.5*chg_here;
	    if ( ( i_kpr+1 >= 0 ) & ( i_kpr+1 < n_p ) )
	      dPg[i_kpr+1] += 0.25*chg_here;
	    if ( ( i_kpr-1 >= 0 ) & ( i_kpr-1 < n_p ) )
	      dPg[i_kpr-1] += 0.25*chg_here;
	  }
      }

  /* Determines dP/dx.  Now, decide on a delta x. */
  
  tot_prob = 0;
  tot_change = 0;
  big_change = 0;
  for ( i_p = 0 ; i_p < n_p ; i_p ++ )
    {
	  tot_prob += P[i_p];
      if ( P[i_p] < 0 )
	printf ( "Severe problem, negative probability! (timestep too large ?) \n" );
      /* Does not seem to happen. */
      tot_change += ABS(dP[i_p]);
      tot_prob += Pg[i_p];
      if ( Pg[i_p] < 0 )
	{
	  printf ( "Severe problem, negative gluon probability! (timestep too large ?) \n" );
	  cout << "for i_p=" << i_p << endl;
	}
      if ( dPg[i_p] < 0 )
	{
	  chg_here = - dPg[i_p] / Pg[i_p];
	  if ( big_change < chg_here )
	    big_change = chg_here;
	}
      if ( dP[i_p] < 0 )
	{
	  chg_here = - dP[i_p] / P[i_p]; /* can't be negative if P[i_p]=0 */
	  if ( big_change < chg_here )
	    big_change = chg_here;
	}
      
    }
  tot_change /= tot_prob;
  if ( big_change > 2 * tot_change ) tot_change = 0.5 * big_change;
  /* Whichever effect is larger--total change, or downward change in
     one location--will decide how big the step size is allowed to be. */
  
  return ( tot_change );
}

double Evolution::evolve_one_step ( double *P , double *Pg , double *Pem , 
				    double *dGamma , double *dGamma_gqq , 
				    double *dGamma_ggg , double *dGamma_em ,  
				    Gamma_info dat , double xmax )
{
  double *dP1 , *dPg1 , *dP2 , *dPg2 , *dPem1 , *dPem2;
  double dx , one_over_dP;
  long   i_p;

  dP1 = (double *) malloc ( sizeof(double) * dat.n_p );
  dPg1 = (double *) malloc ( sizeof(double) * dat.n_p );
  dP2 = (double *) malloc ( sizeof(double) * dat.n_p );
  dPg2 = (double *) malloc ( sizeof(double) * dat.n_p );
  dPem1 = (double *) malloc ( sizeof(double) * dat.n_p );
  dPem2 = (double *) malloc ( sizeof(double) * dat.n_p );

  one_over_dP = find_dP ( P , Pg , Pem , dGamma , dGamma_gqq , dGamma_ggg , 
			  dGamma_em , dP1 , dPg1 , dPem1 , dat );
  dx = dat.delta_x / one_over_dP;  /* step to use. */
  if ( dx > dat.dx_max )
    dx = dat.dx_max;
  if ( dx > xmax ) 
    dx = xmax; /* Don't change too much. */

  for ( i_p = 0 ; i_p < dat.n_p ; i_p ++ )
    {
      P[i_p] += dx * dP1[i_p];
      Pg[i_p] += dx * dPg1[i_p];
      Pem[i_p] += dx * dPem1[i_p];
    }

  /* Find linear order future values. */
  /* Use linear order future values to find future derivatives. */
  find_dP ( P , Pg , Pem , dGamma , dGamma_gqq , dGamma_ggg , 
	    dGamma_em , dP2 , dPg2 , dPem2 , dat );
  /* Now, update; replace dx * dP1, as now performed, with
     dx (dP1 + dP2)/2, by adding dx (dP2-dP1)/2. */
  for ( i_p = 0 ; i_p < dat.n_p ; i_p ++ )
    {
      P[i_p] += 0.5 * dx * ( dP2[i_p] - dP1[i_p]);
      Pg[i_p] += 0.5 * dx * ( dPg2[i_p] - dPg1[i_p]);
      Pem[i_p] += 0.5 * dx * ( dPem2[i_p] - dPem1[i_p]);
    }

  free ( dP1 );
  free ( dPg1 );
  free ( dP2 );
  free ( dPg2 );
  free ( dPem1 );
  free ( dPem2 );

  return ( dx );
}

void Evolution::evolve_in_medium ( double *P , double *Pg , double *Pem , 
				   double *dGamma , double *dGamma_gqq , 
				   double *dGamma_ggg , double *dGamma_em ,  
				   Gamma_info *dat , dGammas *gam , double T , 
				   double beta , double cos_phi , double x_inGeVtoMinusOne )
     /* Evolves in a medium of temperature T, velocity beta, particles
	with angle phi WRT flow velocity vector.  The evolution
	is by an amount, IN GeV^-1, of x.

	First, x in 1/(g^4 GeV).  Then, previous equipment is brought to bear. */
{
  //cout << "T=" << T << endl;
  //cout << "alphas=" << alphas << endl;

  double x_times_GeVg4;
  x_times_GeVg4 = 16 * PI * PI * alphas * alphas * 
     x_inGeVtoMinusOne;
  //cout << " x_inGeVtoMinusOne=" << x_inGeVtoMinusOne << endl;
  //cout << " x_times_GeVg4=" << x_times_GeVg4 << endl;
 
  prep_dGamma ( dat , gam , T , beta , cos_phi , dGamma , dGamma_gqq , 
		dGamma_ggg , dGamma_em );

  while ( x_times_GeVg4 > 0.0000000001 )
    {
      x_times_GeVg4 -= evolve_one_step ( P , Pg , Pem , dGamma , dGamma_gqq , 
				       dGamma_ggg , dGamma_em , 
				       *dat , x_times_GeVg4 );

    }
}

void Evolution::read_table ( Gamma_info *dat , dGammas *Gam )
     /* Reads in the binary stored file of dGamma values. */
{
  FILE *rfile;
  string path="./rates/";
  string pathandfile = path+evolutionsFileNameList->GetFN(9).c_str();
  const char* filename = pathandfile.c_str();
  cout << "opening file " << path+evolutionsFileNameList->GetFN(9).c_str() << " for reading." << endl;
  
  rfile = fopen ( filename , "rb" ); 
  fread ( (char *) (&dat->ddf) , sizeof ( double ) , 1 , rfile );
  fread ( (char *) (&dat->dda) , sizeof ( double ) , 1 , rfile );
  fread ( (char *) (&dat->dcf) , sizeof ( double ) , 1 , rfile );
  fread ( (char *) (&dat->dca) , sizeof ( double ) , 1 , rfile );
  fread ( (char *) (&dat->Nc) , sizeof ( int ) , 1 , rfile );
  fread ( (char *) (&dat->Nf) , sizeof ( int ) , 1 , rfile );
  fread ( (char *) (&dat->BetheHeitler) , sizeof ( int ) , 1 , rfile );
  fread ( (char *) (&dat->BDMPS) , sizeof ( int ) , 1 , rfile );
  fread ( (char *) (&dat->include_gluons) , sizeof ( int ) , 1 , rfile );
  fread ( (char *) Gam->dGamma , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->tau , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->dGamma_gqq , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->tau_gqq , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->dGamma_ggg , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->tau_ggg , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->dGamma_em , sizeof ( double ) , NP * NK , rfile );
  fread ( (char *) Gam->tau_em , sizeof ( double ) , NP * NK , rfile );
  fclose ( rfile );
}

void Evolution::prepare_table ( Gamma_info *dat , dGammas *Gam )
     /* Either reads in or generates gamma table.  Simple driver. */
{
  char c;

  if ( newRadGamma==1 )
    {
      /*
	printf ( "Do you really want to generate a new file for the radiative dGamma/dkdt? (enter y or n) \n" );
	READ_LETTER ( c , stdin );
	if ( c == 'y' )
	{
	build_table ( dat , Gam );
	write_table ( dat , *Gam );
	}
	else
	{
	read_table ( dat , Gam );
	} */
      build_table ( dat , Gam );
      write_table ( dat , *Gam );
    }
  else if ( newRadGamma==0 )
    {
      read_table ( dat , Gam );
    }
}

void Evolution::prep_equipment ( Gamma_info *dat , dGammas *Gam , 
		      double ** dGam1 , double ** dGam2 , 
		      double ** dGam3 , double ** dGam4 )
  // Learns the basic information needed for the radiative part . . . . 
{
  prepare_table ( dat , Gam );
  dat->Nc=Nc;
  dat->Nf=Nf;
  dat->BetheHeitler=BetheHeitler;
  dat->BDMPS=BDMPS;
  dat->dp=dE;
  dat->p_max=maxEnergy;
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

  /* Arranges the k's to be in the range from about -4 GeV
     to about p_max + 4 GeV.  For T=400 MeV, this is enough that
     population functions are only cut off at 10 T. */
  *dGam1 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
  *dGam2 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
  *dGam3 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
  *dGam4 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
  /* Allocate space for the required records. */
  dat->delta_x=dt;
  dat->alpha_s=alphas;
  dat->alpha=1./137.;
}


//-----------------------------------------------------------------------------------//
//-- begin part for collisional energy loss -----------------------------------------//
//-----------------------------------------------------------------------------------//

// here I define the energy loss rates dE/dt for use with the diffusion method

double Evolution::dEdtqq(double T, double En)
{
  return 2./9.*Nf*PI*alphas*alphas*T*T*(log(En*T/mgs(T))+cf+23./12.+cs);
}

double Evolution::dEdtqg(double T, double En)
{
  return 4./3.*PI*alphas*alphas*T*T*(log(En*T/mgs(T))+cb+13./6.+cs);
}

double Evolution::dEdtgq(double T, double En)
{
  return 1./2.*Nf*PI*alphas*alphas*T*T*(log(En*T/mgs(T))+cf+13./6.+cs);
}

double Evolution::dEdtgg(double T, double En)
{
  return 3.*PI*alphas*alphas*T*T*(log(En*T/mgs(T))+cb+131./48.+cs);
}

/// q->g conversion rate
double Evolution::Gammaqg(double T, double En)
{
  return 4./3.*2.*PI*alphas*alphas*T*T/(3.*En)*(0.5*log(En*T/mqs(T))-0.36149);
}

/// g->q conversion rate
double Evolution::Gammagq(double T, double En)
{
  return Nf*3./8.*4./3.*2.*PI*alphas*alphas*T*T/(3.*En)*(0.5*log(En*T/mqs(T))-0.36149);
}

// this is the variation of the quark distribution function per time step for the diffusion method:
double Evolution::dPqdt(double T, double En, double* PqV, double* PgV)
{
  int i = static_cast<int>(En/dE+0.0001);
  return 1./dE*(PqV[i+1]*Gammaqq1(T,En+dE)+PqV[i-1]*Gammaqq2(T,En-dE)-PqV[i]*(Gammaqq1(T,En)+Gammaqq2(T,En)))
    -PqV[i]*Gammaqg(T,En)
    +0.5*(PgV[i+1]*Gammagq(T,En+dE)+PgV[i-1]*Gammagq(T,En-dE));
}

// this is the variation of the gluon distribution function per time step for the diffusion method:
double Evolution::dPgdt(double T, double En, double* PqV, double* PgV)
{
  int i = static_cast<int>(En/dE+0.0001);
  return 1./dE*(PgV[i+1]*Gammagg1(T,En+dE)+PgV[i-1]*Gammagg2(T,En-dE)-PgV[i]*(Gammagg1(T,En)+Gammagg2(T,En)))
    -PgV[i]*Gammagq(T,En)
    +0.5*(PqV[i+1]*Gammaqg(T,En+dE)+PqV[i-1]*Gammaqg(T,En-dE));
;
}

// smooth transition rates (not the diffusion method) for hard q-> hard q
double Evolution::GammaqqS(double T, double En, double omega, double* trqq, double* trqg)
{
  En=En/T*Tfile;// Tfile is the temperature that was used to generate the tabulated transition rate
                // I do not use Tfile=1 because q* in method A depends on T 
                // and may cause too large deviations from the ideal q*
                // for a physically reasonable T.
  omega=omega/T*Tfile;
  
  // the following is used to read from the table with exponentially increasing step sizes 
  int Esize=static_cast<int>((LogEmax-LogEmin)/LogStepE+1);
  int omegaSize=static_cast<int>((2*(LogEmax-LogOmegaMin))/LogStepOmega+3);
  
  int iE=floor((log(En)-LogEmin)/LogStepE);
  if(iE<0) iE=0;
  
  int iOmega;
  if(omega>0) iOmega=floor((log(omega)-LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega+1);
  if(omega<0) iOmega=floor((-log(-omega)+LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega);

  double fracO;
  if(omega>0) 
    fracO=(omega-(exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))))
      /(exp(((iOmega)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))
	-exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin)));
  
  if(omega<0)
    fracO=(omega-(-exp((-(iOmega)*LogStepOmega+LogEmax))))
      /((exp((-(iOmega)*LogStepOmega+LogEmax)))-exp((-(iOmega+1)*LogStepOmega+LogEmax)));
  
  double fracE=(En-exp((LogEmin+iE*LogStepE)))/(exp((LogEmin+(iE+1)*LogStepE))-exp((LogEmin+iE*LogStepE)));

  return ((1-fracE)*((1-fracO)*trqq[iE*(omegaSize)+iOmega]+fracO*trqq[iE*(omegaSize)+(iOmega+1)])
    +fracE*((1-fracO)*trqq[(iE+1)*(omegaSize)+iOmega]+fracO*trqq[(iE+1)*(omegaSize)+(iOmega+1)])+
    (1-fracE)*((1-fracO)*trqg[iE*(omegaSize)+iOmega]+fracO*trqg[iE*(omegaSize)+(iOmega+1)])
    +fracE*((1-fracO)*trqg[(iE+1)*(omegaSize)+iOmega]+fracO*trqg[(iE+1)*(omegaSize)+(iOmega+1)]));
}

// smooth transition rates (not the diffusion method) for hard g-> hard g
double Evolution::GammaggS(double T, double En, double omega, double* trgq, double* trgg)
{
  En=En/T*Tfile;// Tfile is the temperature that was used to generate the tabulated transition rate
                // I do not use Tfile=1 because q* in method A depends on T 
                // and may cause too large deviations from the ideal q*
                // for a physically reasonable T.
  omega=omega/T*Tfile;

  // the following is used to read from the table with exponentially increasing step sizes 
  int Esize=static_cast<int>((LogEmax-LogEmin)/LogStepE+1);
  int omegaSize=static_cast<int>((2*(LogEmax-LogOmegaMin))/LogStepOmega+3);

  int iE=floor((log(En)-LogEmin)/LogStepE);
  if(iE<0) iE=0;
  
  int iOmega;
  if(omega>0) iOmega=floor((log(omega)-LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega+1);
  if(omega<0) iOmega=floor((-log(-omega)+LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega);
  
  double fracO;
  if(omega>0) 
    fracO=(omega-(exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))))
      /(exp(((iOmega)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))
	-exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin)));
  
  if(omega<0)
    fracO=(omega-(-exp((-(iOmega)*LogStepOmega+LogEmax))))
      /((exp((-(iOmega)*LogStepOmega+LogEmax)))-exp((-(iOmega+1)*LogStepOmega+LogEmax)));
  
  double fracE=(En-exp((LogEmin+iE*LogStepE)))/(exp((LogEmin+(iE+1)*LogStepE))-exp((LogEmin+iE*LogStepE)));
  
  return ((1-fracE)*((1-fracO)*trgq[iE*(omegaSize)+iOmega]+fracO*trgq[iE*(omegaSize)+(iOmega+1)])
    +fracE*((1-fracO)*trgq[(iE+1)*(omegaSize)+iOmega]+fracO*trgq[(iE+1)*(omegaSize)+(iOmega+1)])+
    (1-fracE)*((1-fracO)*trgg[iE*(omegaSize)+iOmega]+fracO*trgg[iE*(omegaSize)+(iOmega+1)])
    +fracE*((1-fracO)*trgg[(iE+1)*(omegaSize)+iOmega]+fracO*trgg[(iE+1)*(omegaSize)+(iOmega+1)]));
}

// this is the variation of the quark distribution function per time step for the full method A or B:
double Evolution::dPqdtS(double T, double En, double* PqV, double* PgV, double* trqq, double* trqg, double* trgq, double* trgg)
{
  int i = static_cast<int>(En/dE+0.0001);
  int imaxE = static_cast<int>(maxEnergy/dE+0.0001);
  double value=0.;
  double omega;

  int omegaSize=static_cast<int>((2*(LogEmax-LogOmegaMin))/LogStepOmega+2);
  

  for(int j=-i;j<imaxE-i;j++)
    {
      omega=j*dE;
      if(En+omega>exp(LogEmin) && En>exp(LogEmin) && omega!=0 && omega<=exp(LogEmax)) 
	{
	  value+=dE*(PqV[i+j]*GammaqqS(T,En+omega,omega,trqq,trqg)-PqV[i]*GammaqqS(T,En,omega,trqq,trqg));
	}
    }
  value+=-PqV[i]*Gammaqg(T,En)+0.5*(PgV[i+1]*Gammagq(T,En+dE)+PgV[i-1]*Gammagq(T,En-dE));// conversion
  
  return value; 
}

// this is the variation of the gluon distribution function per time step for the full method A or B:
double Evolution::dPgdtS(double T, double En, double* PqV, double* PgV, double* trqq, double* trqg, double* trgq, double* trgg)
{
  int i = static_cast<int>(En/dE+0.0001);
  int imaxE = static_cast<int>(maxEnergy/dE+0.0001);
  double value=0.;
  double omega;

  int omegaSize=static_cast<int>((2*(LogEmax-LogOmegaMin))/LogStepOmega+2);


  for(int j=-i;j<imaxE-i;j++)
    {
      omega=j*dE;
      if(En+omega>exp(LogEmin) && En>exp(LogEmin) && omega!=0 && omega<=exp(LogEmax))
	{
	  value+=dE*(PgV[i+j]*GammaggS(T,En+omega,omega,trgq,trgg)-PgV[i]*GammaggS(T,En,omega,trgq,trgg));
	}    
    }
  value+=-PgV[i]*Gammagq(T,En)+0.5*(PqV[i+1]*Gammaqg(T,En+dE)+PqV[i-1]*Gammaqg(T,En-dE)); //conversion

  return value; 
}

// write dE/dt of the four processes into file. not really needed - just as a test.
void Evolution::dEdtOutput()
{

  ofstream fout(evolutionsFileNameList->GetFN(2).c_str(),ios::out); 
 
  int start = static_cast<int>(4/dE);
  int stop = static_cast<int>(50/dE+1);

  fout << "E " << " dEdtqq " << " dEdtqg " << " dEdtgq " << " dEdtgg " << endl;
  
  for(int i=start; i<stop;i++)
    {
      fout << dE*i 
	   << " " << dEdtqq(T,dE*i)/hbarc 
	   << " " << dEdtqg(T,dE*i)/hbarc 
	   << " " << dEdtgq(T,dE*i)/hbarc 
	   << " " << dEdtgg(T,dE*i)/hbarc 
	   << endl;
    }
}

// initialize a single quark jet with energy initialE:
void Evolution::initP(double* Pq, double* Pg, double* Pem)
{
  int esize = static_cast<int>(maxEnergy/dE+0.0001);

  for(int i=0;i<=esize;i++)
    {
      Pq[i]=exp(-(i*dE-initialE)*(i*dE-initialE)/(2*width*width))/(sqrt(2*PI)*width);
      Pg[i]=0;
      Pem[i]=0;
    }
}

//here we go one step in time in the diffusion method
void Evolution::run(double* Pq, double* Pg, double* Pem, double* PqA, double* PgA)
{
  int esize = static_cast<int>(maxEnergy/dE+0.0001);
  int estart = static_cast<int>(2/dE+1);
  //cout << "dPqdt(T,16)=" << dPqdt(T,16,Pq,Pg) << " Pq(16)=" << Pq[static_cast<int>(16/dE+0.0001)]<< endl;

  for(int i=0;i<static_cast<int>(2.5/dE);i++)
    {
      Pq[i]=0;
      Pg[i]=0;
    }   

  if(doRad) 
    {
      if(counter==0) cout << "do radiative" << endl;
      evolve_in_medium(Pq, Pg, Pem, dGamma, dGamma_gqq, dGamma_ggg, dGamma_em, &dat, &Gam, T, 0, 1, dt);
    }
  
  for(int i=0;i<static_cast<int>(2.5/dE);i++)
    {
      Pq[i]=0;
      Pg[i]=0;
    }   

  counter++; //counter just for output of line "do ... " - nothing important
    
  if(doCol)
    {

      //determine time step:

      double dtmax;
      double mydt;
      double myt=0.;
      double totchange=0.;
      double totprob=0.;
      double chg_here;
      double big_change=0;
      double tot_gamma=0;

      if(counter<2) 
	{
	  cout << "do collisional" << endl;
	}
      
      for(int i=estart;i<esize;i++)
	{
	  double En = i*dE;
	  tot_gamma += Gammaqq1(T,En)+Gammaqq2(T,En)+Gammagg1(T,En)+Gammagg2(T,En);
	}
      dtmax=0.8/(tot_gamma*2*dE);

      
      //cout << "dtmax=" << dtmax << endl;
      //cout << "dt=" << dt << endl;
      
      for(int i=estart;i<esize;i++)
	{
	  double En = i*dE;
	  double add=dPqdt(T,En,Pq,Pg);
	  totprob+=Pq[i]+Pg[i];	  
	  if(add<0) add=-add;
	  totchange+=add;
	  if ( dPqdt(T,En,Pq,Pg) < 0 )
	    {
	      chg_here = - dPqdt(T,En,Pq,Pg) / Pq[i];
	      if ( big_change < chg_here )
		big_change = chg_here;
	    } 
	  if ( dPgdt(T,En,Pq,Pg) < 0 )
	    {
	      chg_here = - dPgdt(T,En,Pq,Pg) / Pg[i];
	      if ( big_change < chg_here )
		big_change = chg_here;
	    } 
	}

      if(totprob!=0) totchange/=totprob;
      if(big_change>2*totchange) totchange=0.5*big_change;
      
      
      mydt=dtmax*1/dE/totchange;
      
      if(mydt>dt) mydt=dt;
      if(mydt>0.08*dE) mydt=0.08*dE;
      //cout << " dt=" << dt << endl;
      //cout << " mydt=" << mydt << endl;
      int imax= static_cast<int>(dt/mydt);
      //cout << "imax=" << imax << endl;

      for(int i=0;i<=esize;i++)
	{
	  PqA[i]=0;
	  PgA[i]=0;
	}

      // time step determined. now integrate / solve DE:

      ///second order Runge Kutta:
      for(int it=0;it<imax;it++)
	{
	  for(int i=estart;i<esize;i++)
	    {
	      double En = i*dE;
	      PqA[i]=Pq[i]+0.5*(mydt*dPqdt(T,En,Pq,Pg));
	      PgA[i]=Pg[i]+0.5*(mydt*dPgdt(T,En,Pq,Pg));
	    }
	  
	  for(int i=estart;i<esize;i++)
	    {
	      double En=i*dE;     
	      Pq[i]=Pq[i]+mydt*dPqdt(T,En,PqA,PgA);
	      Pg[i]=Pg[i]+mydt*dPgdt(T,En,PqA,PgA);
	    }
	}
      double diff=dt-(imax)*mydt; 
      //cout << "difference dt-(imax)*mydt=" << diff << endl;
	  for(int i=estart;i<esize;i++)
	    {
	      double En = i*dE;
	      PqA[i]=Pq[i]+0.5*(diff*dPqdt(T,En,Pq,Pg));
	      PgA[i]=Pg[i]+0.5*(diff*dPgdt(T,En,Pq,Pg));
	    }
	  
	  for(int i=estart;i<esize;i++)
	    {
	      double En=i*dE;     
	      Pq[i]=Pq[i]+diff*dPqdt(T,En,PqA,PgA);
	      Pg[i]=Pg[i]+diff*dPgdt(T,En,PqA,PgA);
	    }
    }
}

//here we go one step in time using the full collisional rates (method A or B)
void Evolution::runsmooth(double* Pq, double* Pg, double* Pem, double* PqA, double* PgA, double* trqq, double* trqg, double* trgq, double* trgg)
{
  int esize = static_cast<int>(maxEnergy/dE+0.0001);
  int estart = static_cast<int>(exp(LogEmin)/dE+1+0.0001);
  //cout << "dPqdt(T,16)=" << dPqdt(T,16,Pq,Pg) << " Pq(16)=" << Pq[static_cast<int>(16/dE+0.0001)]<< endl;

  for(int i=0;i<static_cast<int>(2.5/dE);i++)
    {
      Pq[i]=0;
      Pg[i]=0;
    }   
  

  if(doRad) 
    {
      if(counter==0) cout << "do radiative" << endl;
      evolve_in_medium(Pq, Pg, Pem, dGamma, dGamma_gqq, dGamma_ggg, dGamma_em, &dat, &Gam, T, 0, 1, dt);
    }
  counter++;

  for(int i=0;i<static_cast<int>(2.5/dE);i++)
    {
      Pq[i]=0;
      Pg[i]=0;
    }   
  
  if(doCol)
    {
      if(counter<2) cout << "do collisional" << endl;


      //determine time step:
      double dtmax;
      double mydt;
      /* (dont bother for now... works fine with the given dt)
      double myt=0.;
      double totchange=0.;
      double totprob=0.;
      double chg_here;
      double big_change=0;
      double tot_gamma=0;
      double value=0;
      double omega;
      int imaxE = maxEnergy/dE+0.0001;


      if(counter<2) 
	{
	  cout << "do collisional" << endl;
	}
      
      for(int i=estart;i<esize;i++)
	{
	  double En = i*dE;
	  for(int j=-i;j<imaxE-i;j++)
	    {
	      omega=j*dE;
	      if(En>exp(LogEmin) && omega!=0 && omega<=exp(LogEmax)) 
		{
		  value+=dE*(GammaqqS(T,En,omega,trqq,trqg)+GammaggS(T,En,omega,trgq,trgg));
		}
	    }
	  	  
	  tot_gamma += value;
	}
      dtmax=100/(tot_gamma*2*dE);

      cout << "dtmax=" << dtmax << endl;
      cout << "dt=" << dt << endl;
      
      for(int i=estart;i<esize;i++)
	{
	  double En = i*dE;
	  double add=dPqdt(T,En,Pq,Pg);
	  totprob+=Pq[i]+Pg[i];	  
	  if(add<0) add=-add;
	  totchange+=add;
	  if ( dPqdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg) < 0 )
	    {
	      chg_here = - dPqdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg) / Pq[i];
	      if ( big_change < chg_here )
		big_change = chg_here;
	    } 
	  if ( dPgdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg) < 0 )
	    {
	      chg_here = - dPgdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg) / Pg[i];
	      if ( big_change < chg_here )
		big_change = chg_here;
	    } 
	}

      //cout << "totchange=" << totchange << endl;
      //cout << "totprob=" << totprob << endl;

      if(totprob!=0) totchange/=totprob;
      if(big_change>2*totchange) 
	{
	  totchange=0.5*big_change;
	  cout << "using big change!" << endl;
	}
      
      //cout << "totchange=" << totchange << endl;

      mydt=dtmax*1/dE/totchange;

      if(mydt>dt) mydt=dt;
      cout << " dt=" << dt << endl;
      cout << " mydt=" << mydt << endl;
      */

      //cout << "imax=" << imax << endl;
      
      mydt=dt;
      int imax= static_cast<int>(dt/mydt);

      for(int i=0;i<=esize;i++)
	{
	  PqA[i]=0;
	  PgA[i]=0;
	}

   
      ///second order Runge Kutta:
      for(int it=0;it<imax;it++)
	{
	  for(int i=estart;i<esize;i++)
	    {
	      double En=i*dE;     
	      PqA[i]=Pq[i]+0.5*(mydt*dPqdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg));
	      PgA[i]=Pg[i]+0.5*(mydt*dPgdtS(T,En,Pq,Pg,trqq,trqg,trgq,trgg));
	    }
	  
	  for(int i=estart;i<esize;i++)
	    {
	      double En=i*dE;     
	      Pq[i]=Pq[i]+mydt*dPqdtS(T,En,PqA,PgA,trqq,trqg,trgq,trgg);
	      Pg[i]=Pg[i]+mydt*dPgdtS(T,En,PqA,PgA,trqq,trqg,trgq,trgg);
	    }
	}
      double diff=dt-(imax)*mydt; 
      //cout << "difference dt-(imax)*mydt=" << diff << endl;
      for(int i=estart;i<esize;i++)
	{
	  double En = i*dE;
	  PqA[i]=Pq[i]+0.5*(diff*dPqdt(T,En,Pq,Pg));
	  PgA[i]=Pg[i]+0.5*(diff*dPgdt(T,En,Pq,Pg));
	}
      
      for(int i=estart;i<esize;i++)
	{
	  double En=i*dE;     
	  Pq[i]=Pq[i]+diff*dPqdt(T,En,PqA,PgA);
	  Pg[i]=Pg[i]+diff*dPgdt(T,En,PqA,PgA);
	}
    }
}
