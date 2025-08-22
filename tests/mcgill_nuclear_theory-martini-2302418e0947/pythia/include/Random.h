
#ifndef Random_h
#define Random_h

#include <stdio.h>
#include <iostream>

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

class Random
{
 private:
  /* The array for the state vector */
  unsigned long long mt[NN]; 
  /* mti==NN+1 means mt[NN] is not initialized */
  int mti;
 public:
  Random();//constructor
  ~Random();//destructor
  void init_genrand64(unsigned long long seed);
  void init_by_array64(unsigned long long init_key[],
		       unsigned long long key_length);
  unsigned long long genrand64_int64(void);
  long long genrand64_int63(void);
  double genrand64_real1(void);
  double genrand64_real2(void);
  double genrand64_real3(void);
};

#endif
