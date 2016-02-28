// Random number generator


#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
namespace nmlib{


#include "mt64.h"  // 64bit Mersenne Twister


// RNG (random number generator)
inline void init_rand(unsigned long long seed){ init_genrand64(seed); }         // initialization
inline unsigned long long irand(void){ return genrand64_int64(); }              // [0,2^64) uniform
inline double urand(void){ return genrand64_real3(); }                          // U(0,1) uniform
inline double nrand(void){ return sqrt(-2*log(urand()))*cos(2*M_PI*urand()); }  // N(0,1) gaussian
inline double erand(void){ return -log(urand()); }                              // Ex(1)  exponential


// LDS(low-discrepancy sequence) - van der Corput
inline double lds(ulong n, ulong b){
  ulong q=0, p=1;
  while(p<=ulong(-1)/b){ q=q*b+(n%b); p*=b; n/=b; }
  return (q+0.5)/p;  // n=(ABCDE)b --> (0.EDCBA)b + 1/2max
}


}  //namespace nmlib
#endif //RANDOM_H
