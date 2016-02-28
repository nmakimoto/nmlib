// Random number generator


#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
#include "matrix.H"
namespace nmlib{


#include "mt64.h"  // 64bit Mersenne Twister


// RNG (random number generator)
inline void init_rand(unsigned long long seed){ init_genrand64(seed); }         // initialization
inline unsigned long long irand(void){ return genrand64_int64(); }              // [0,2^64) uniform
inline double urand(void){ return genrand64_real3(); }                          // U(0,1) uniform
inline double nrand(void){ return sqrt(-2*log(urand()))*cos(2*M_PI*urand()); }  // N(0,1) gaussian
inline double erand(void){ return -log(urand()); }                              // Ex(1)  exponential
inline Matrix urand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=urand(); return x; }  // multidim uniform
inline Matrix nrand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=nrand(); return x; }  // multidim normal


// LDS(low-discrepancy sequence) - van der Corput
inline double corput(ulong n, ulong b){
  ulong q=0, p=1;
  while(p<=ulong(-1)/b){ q=q*b+(n%b); p*=b; n/=b; }
  return (q+0.5)/p;  // n=(ABCDE)b --> (0.EDCBA)b + 1/2max
}


class Corput{
public:
  Corput(ulong base=2, ulong seqno=1): b(base),s(seqno) {}
  double operator()(void) const { return corput(s++,b); }
private:
  ulong b;
  mutable ulong s;
};


class Halton{
public:
  Halton(size_t dim, ulong seqno=1): s(seqno) {
    if(dim>5) throw std::runtime_error("Halton(dim): currently dim must be <=5");
    ulong pp[]={2,3,5,7,11};
    bb=std::vector<ulong>(pp,pp+dim);
  }
  Matrix operator()(void) const {
    Matrix v(bb.size());
    for(size_t k=0; k<bb.size(); k++) v(k)=corput(s,bb[k]);
    s++;
    return v;
  }
private:
  std::vector<ulong> bb;
  mutable ulong s;
};


}  //namespace nmlib
#endif //RANDOM_H
