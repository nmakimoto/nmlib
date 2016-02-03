// Random number generator


#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
namespace nmlib{


#include "mt64.h"  // 64bit Mersenne Twister


inline void init_rand(unsigned long long seed){ init_genrand64(seed); }         // initialization
inline unsigned long long irand(void){ return genrand64_int64(); }              // [0,2^64) uniform
inline double urand(void){ return genrand64_real3(); }                          // U(0,1) uniform
inline double nrand(void){ return sqrt(-2*log(urand()))*cos(2*M_PI*urand()); }  // N(0,1) gaussian
inline double erand(void){ return -log(urand()); }                              // Ex(1)  exponential

inline double lds(unsigned n, unsigned b=2){  // low-discrepancy sequence (van der Corput)
  unsigned g=0, bl=1;
  while(n){ g*=b; g+=(n%b); bl*=b; n/=b; }  // revert digits of b-ary n
  return (g+0.5)/bl;  // normalize to (0,1)
}


}  //namespace nmlib
#endif //RANDOM_H
