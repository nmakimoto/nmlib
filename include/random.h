// Random number generators and low-discrepancy sequences


#ifndef RANDOM_H
#define RANDOM_H


#include <random>
#include <vector>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include "matrix.h"
namespace nmlib{


/******** Class I/F ********/

// Random number generator
class Rng {
public:
  typedef uint_fast64_t Uint;

  Rng      (Uint seed=0);
  void init(Uint seed=0);
  double operator()(void)     const;  // U(0,1)
  Matrix operator()(size_t d) const;  // U(0,1)^d

  double u(double a=0, double b=1)         const;  // U(a,b)
  double n(double m=0, double s=1)         const;  // N(m,s^2)
  double e(double m=1)                     const;  // Ex(m)
  Matrix u(double a,   double b, size_t d) const;  // U(a,b)^d
  Matrix n(double m,   double s, size_t d) const;  // N(m,s)^d
  template<class Func>
    double rand_pdf(const Func& pdf, double pmax, double a, double b) const;  // given PDF on (a,b)

private:
  mutable std::mt19937_64 rng;
};


// Low-discrepancy sequence
class Lds {
public:
  typedef uint_fast64_t     Uint;
  typedef std::vector<Uint> Uvec;

  Lds      (Uint seqno=0);
  void init(Uint seqno=0);
  double operator()(      Uint  base=2) const;  // U(0,1)
  Matrix operator()(const Uvec& base  ) const;  // U(0,1)^d

  double u(double a=0, double b=1,       Uint  base=2)     const;  // U(a,b)
  double n(double m=0, double s=1, const Uvec& base={2,3}) const;  // N(m,s)
  double e(double m=1,                   Uint  base=2)     const;  // Ex(m)
  Matrix u(double a,   double b,   const Uvec& base)       const;  // U(a,b)^d
  //Matrix n(double m=0, double s=1, const Uvec& base={2,3}) const;  // N(m,s)^d - not supported
  template<class Func>
    double rand_pdf(const Func& pdf, double pmax, double a, double b, const Lds::Uvec& base) const;  // given PDF on (a,b)

private:
  void incr(void) const { seq++; if(seq==0) seq++; }
  mutable uint_fast64_t seq;
};


/******** Utility I/F ********/

// Low-discrepancy sequence
double corput(Lds::Uint seq,       Lds::Uint  base=2);  // 1dim LDS - van der Corput sequence
Matrix halton(Lds::Uint seq, const Lds::Uvec& base  );  // multidim LDS - Halton sequence


// Conversion utils
double box_muller(const Matrix& u);  // U(I^2) -> N(0,1) converter
Matrix conv_i2s2 (const Matrix& u);  // U(I^2) -> U(S^2) converter (returns a unit 3-vector)


// The following functions are not supported. Use global Rng/Lds objects instead.
/*
// RNG functions
extern Rng RNG0;  // must be defined and linked in the user code
inline void   init_rand(Rng::Uint seed=0){ RNG0.init(seed); }
inline double urand(double a=0, double b=1){ return RNG0.u(a,b); }  // U(a,b)
inline double nrand(double m=0, double s=1){ return RNG0.n(m,s); }  // N(m,s^2)
inline double erand(double m=1)            { return RNG0.e(m); }    // Ex(m)
inline Matrix urand(double a,   double b, size_t d){ return RNG0.u(a,b,d); }  // U(a,b)^d
inline Matrix nrand(double m,   double s, size_t d){ return RNG0.n(m,s,d); }  // N(m,s^2)^d
template<class Func> double rand_pdf(const Func& pdf, double pmax, double a, double b){ return RNG0.rand_pdf(pdf,pmax,a,b); }  // PDF on (a,b)

// LDS functions
extern Lds LDS0;  // must be defined and linked in the user code
inline void   init_lds(Lds::Uint seqno=0){ LDS0.init(seqno); }
inline double u_lds(double a=0, double b=1,       Lds::Uint  base=2    ){ return LDS0.u(a,b,base); }  // U(a,b)
inline double n_lds(double m=0, double s=1, const Lds::Uvec& base={2,3}){ return LDS0.n(m,s,base); }  // N(m,s^2)
inline double e_lds(double m=1,                   Lds::Uint  base=2    ){ return LDS0.u(m,  base); }  // Ex(m)
inline Matrix u_lds(double a,   double b,   const Lds::Uvec& base      ){ return LDS0.u(a,b,base); }  // U(a,b)^d
//inline Matrix n_lds(double m, double s, const Lds::Uvec& base){ return LDS0.n(m,s,base); }  // N(m,s^2)
template<class Func> double rand_pdf(const Func& pdf, double pmax, double a, double b, const Lds::Uvec& base){ return LDS0.rand_pdf(pdf,pmax,a,b,base); }  // PDF on (a,b)
*/


/******** Implementation ********/

// RNG class methods
inline        Rng::Rng (Rng::Uint seed){ init(seed); }
inline void   Rng::init(Rng::Uint seed){ rng=std::mt19937_64(seed ? seed : std::random_device()()); }
inline double Rng::operator()(void)      const { return u(0,1); }
inline Matrix Rng::operator()(size_t d)  const { return u(0,1,d); }
inline double Rng::u(double a, double b) const { std::uniform_real_distribution<double> dist(a,b); return dist(rng); }
inline double Rng::n(double m, double s) const { std::normal_distribution      <double> dist(m,s); return dist(rng); }
inline double Rng::e(double m)           const { std::exponential_distribution <double> dist(m);   return dist(rng); }
inline Matrix Rng::u(double a, double b, size_t d) const { Matrix x(d); for(size_t i=0; i<d; i++) x(i)=u(a,b); return x; }
inline Matrix Rng::n(double m, double s, size_t d) const { Matrix x(d); for(size_t i=0; i<d; i++) x(i)=n(m,s); return x; }
template<class Func> double Rng::rand_pdf(const Func& pdf, double pmax, double a, double b) const {
  double x;
  while(1){
    x=u(a,b);
    if( u(0,pmax)<pdf(x) ) break;  // rejection sampling
  }
  return x;
}


// LDS class methods
inline        Lds::Lds (Lds::Uint seqno){ init(seqno); }
inline void   Lds::init(Lds::Uint seqno){ seq=seqno; }
inline double Lds::operator()(      Lds::Uint  base)            const { incr(); return corput(seq,base); }
inline Matrix Lds::operator()(const Lds::Uvec& base)            const { incr(); return halton(seq,base); }
inline double Lds::u(double a, double b,       Lds::Uint  base) const { incr(); return corput(seq,base)*(b-a)+a; }
inline double Lds::n(double m, double s, const Lds::Uvec& base) const { incr(); return box_muller(halton(seq,base))*s+m; }
inline double Lds::e(double m,                 Lds::Uint  base) const { incr(); return -log(corput(seq,base))*m; }
inline Matrix Lds::u(double a, double b, const Lds::Uvec& base) const {
  incr();
  Matrix x=halton(seq,base);
  for(size_t i=0; i<x.dim(); i++) x(i)=x(i)*(b-a)+a;
  return x;
}
template<class Func> double Lds::rand_pdf(const Func& pdf, double pmax, double a, double b, const Lds::Uvec& base) const {
  double x;
  while(1){
    incr();
    x=corput(seq,base[0])*(b-a)+a;
    if( corput(seq,base[1])*pmax < pdf(x) ) break;  // rejection sampling
  }
  return x;
}


// Utilities
inline double corput(Lds::Uint n, Lds::Uint b){
  Lds::Uint q=0, p=1;
  while(n){ q=q*b+(n%b); p*=b; n/=b; }
  return double(q)/p;  // n=(ABCDE)b --> (0.EDCBA)b
}
inline Matrix halton(Lds::Uint n, const Lds::Uvec& base){
  size_t dim=base.size();
  Matrix x(dim);
  for(size_t k=0; k<dim; k++) x(k)=corput(n,base[k]);
  return x;
}
inline double box_muller(const Matrix& u){
  if( u.dim()<2 || u(0)<0 ) throw std::domain_error("box_muller: domain error");
  return sqrt(-2*log(u(0)))*cos(2*M_PI*u(1));
}
inline Matrix conv_i2s2(const Matrix& u){
  if( u.dim()<2 || u(0)<0 || 1<u(0) ) throw std::domain_error("conv_i2s2: domain error");
  double th=asin(u(0)*2-1), ph=2*M_PI*u(1);
  return Matrix({cos(th)*cos(ph), cos(th)*sin(ph), sin(th)});
}


}  //namespace nmlib
#endif //RANDOM_H
