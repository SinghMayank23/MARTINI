// Random.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the random generator, which is

/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

#include "Random.h"

Random::Random()
{
}

Random::~Random()
{
}

/* initializes mt[NNM] with a seed */
void Random::init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NNM; mti++) 
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void Random::init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length)
{
    unsigned long long i, j, k;
    init_genrand64(19650218ULL);
    i=1; j=0;
    k = (NNM>key_length ? NNM : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NNM) { mt[0] = mt[NNM-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NNM-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
          - i; /* non linear */
        i++;
        if (i>=NNM) { mt[0] = mt[NNM-1]; i=1; }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long Random::genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NNM) { /* generate MM2 words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NNM+1) 
            init_genrand64(5489ULL); 

        for (i=0;i<NNM-MM2;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM2] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NNM-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM2-NNM)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NNM-1]&UM)|(mt[0]&LM);
        mt[NNM-1] = mt[MM2-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }
  
    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long Random::genrand64_int63(void)
{
    return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double Random::genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double Random::genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double Random::genrand64_real3(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}

Vec4 Random::thermal(double temp, int kind)
{
  // this algorithm uses 5.5/(1+px2) as envelope function and then uses the rejection method
  // tan (Pi/2*x) is the inverse of the integral of the envelope function
  // this can be improved
  // kind: -1 = Boson
  //        1 = Fermion
  Vec4 p;
  double gm = 5.5;
  int repeat = 1;
  double gx;
  double px;
  double cost, sint;
  double phi;
  double px2;
  double ex;
  do 
    {
      cost = 2*genrand64_real1() - 1.0;
      sint = sqrt(1.0 - cost*cost);
      phi = 2*M_PI*genrand64_real1();
      px = tan (M_PI * genrand64_real1() / 2.0);
      px2 = px*px;
      ex = sqrt (px2);
      if (kind == 0)  gx = px2 * (1.0 + px2) * exp (- ex);
      else if (kind == -1) gx = px2 * (1.0 + px2) * 1/(exp(ex)-1);
      else if (kind == +1) gx = px2 * (1.0 + px2) * 1/(exp(ex)+1);
      if ( genrand64_real1() < gx / gm) 
	{
	  repeat = 0;
	  p.px(px * sint * cos(phi));
	  p.py(px * sint * sin(phi));
	  p.pz(px * cost );
	  p.e(ex);
      	} 
      else if (gx > gm) 
	{
	}
    } while (repeat);   
  p = p*temp;
  return p;
}

// sampling thermal parton that makes a recoil parton on-shell.
double Random::thermal2(double k_min, double temp, int kind)
{
  // this algorithm uses 5.5/(1+px2) as envelope function and then uses the rejection method
  // tan (Pi/2*x) is the inverse of the integral of the envelope function
  // this can be improved
  // kind: -1 = Boson
  //        1 = Fermion
  double p;
  double gm = 5.5;
  int repeat = 1;
  double gx;
  double px;
  double phi;
  double px2;
  double ex;

  int i = 0;
  do 
    {
      i++;
      phi = 2*M_PI*genrand64_real1();
      px = k_min/temp + tan (M_PI * genrand64_real1() / 2.0);
      px2 = px*px;
      ex = sqrt (px2);
      if (kind == 0) gx = px2 * (1.0 + px2) * exp (- ex);
      else if (kind == -1) gx = px2 * (1.0 + px2) * 1/(exp(ex)-1);
      else if (kind == +1) gx = px2 * (1.0 + px2) * 1/(exp(ex)+1);
      if ( genrand64_real1() < gx / gm)
        {
          repeat = 0;
          p = ex;
      	} 
      else if (gx > gm) 
	{
	}
    } while (repeat);   
  
  p = p*temp;
  return p;
}
