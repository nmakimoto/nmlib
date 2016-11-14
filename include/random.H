// Random number generator


#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
#include "matrix.H"
namespace nmlib{


#include "mt64.h"  // 64bit Mersenne Twister


/******** Utilitiy I/F ********/

void init_rand(unsigned long long seed);  // initialization
unsigned long long irand(void);  // [0,2^64) uniform
double urand(void);  // U(0,1) uniform
double nrand(void);  // N(0,1) gaussian
double erand(void);  // Ex(1)  exponential
Matrix urand_m(size_t n);  // multidim uniform U(0,1)^n
Matrix nrand_m(size_t n);  // multidim gaussian N(0,I)

double corput(ulong seq, ulong base=2);  // LDS(low-discrepancy sequence) - van der Corput
Matrix halton(ulong seq, size_t dim);    // multidim LDS


/******** Implementation ********/

inline void init_rand(unsigned long long seed){ init_genrand64(seed); }
inline unsigned long long irand(void){ return genrand64_int64(); }
inline double urand(void){ return genrand64_real3(); }
inline double nrand(void){ return sqrt(-2*log(urand()))*cos(2*M_PI*urand()); }
inline double erand(void){ return -log(urand()); }
inline Matrix urand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=urand(); return x; }
inline Matrix nrand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=nrand(); return x; }

inline double corput(ulong n, ulong b){
  ulong q=0, p=1;
  while(p<=ulong(-1)/b){ q=q*b+(n%b); p*=b; n/=b; }
  return (q+0.5)/p;  // n=(ABCDE)b --> (0.EDCBA)b + 1/2max
}
inline Matrix halton(ulong n, size_t dim){
  const ulong pp[]={2,3,5,7,11,13};
  if(dim>sizeof(pp)/sizeof(pp[0])) throw std::runtime_error("halton(): dimension too large (>6)");
  Matrix x(dim);
  for(size_t k=0; k<dim; k++) x(k)=corput(n,pp[k]);
  return x;
}


}  //namespace nmlib
#endif //RANDOM_H
