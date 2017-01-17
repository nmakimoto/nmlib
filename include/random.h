// Random number generator


#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
#include <cstdint>
#include "matrix.h"
namespace nmlib{


/******** Utilitiy I/F ********/

void init_rand(uint_fast64_t seed=0);  // initialization
uint_fast64_t irand(void);  // [0,2^64) uniform
double urand(void);  // U(0,1) uniform
double nrand(void);  // N(0,1) gaussian
double erand(void);  // Ex(1)  exponential
Matrix urand_m(size_t n);  // multidim uniform U(0,1)^n
Matrix nrand_m(size_t n);  // multidim gaussian N(0,I)

double corput(uint_fast64_t seq, uint_fast64_t base=2);  // LDS(low-discrepancy sequence) - van der Corput
Matrix halton(uint_fast64_t seq, size_t dim);    // multidim LDS


/******** Implementation ********/

inline Matrix urand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=urand(); return x; }
inline Matrix nrand_m(size_t n){ Matrix x(n); for(size_t k=0; k<n; k++) x(k)=nrand(); return x; }

inline double corput(uint_fast64_t n, uint_fast64_t b){
  uint_fast64_t q=0, p=1;
  while(p<=uint_fast64_t(-1)/b){ q=q*b+(n%b); p*=b; n/=b; }
  return double(q)/p;  // n=(ABCDE)b --> (0.EDCBA)b + 1/2max
}
inline Matrix halton(uint_fast64_t n, size_t dim){
  const uint_fast64_t pp[]={2,3,5,7,11,13};
  if(dim>sizeof(pp)/sizeof(pp[0])) throw std::runtime_error("halton(): dimension too large (>6)");
  Matrix x(dim);
  for(size_t k=0; k<dim; k++) x(k)=corput(n,pp[k]);
  return x;
}


}  //namespace nmlib
#endif //RANDOM_H
