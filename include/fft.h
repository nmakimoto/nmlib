// FFT(Fast Fourier Transformation) template library
// Calculates XX=fft(xx) etc. for a vector<complex<T>> xx of 2^nu dim.


#ifndef FFT_H
#define FFT_H


#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace nmlib{


/******** Prototype ********/

template<class Vec> Vec fft    (const Vec& xx, bool inverse=false);  // FFT
template<class Vec> Vec ifft   (const Vec& ff);                      // inverse FFT
template<class Vec> Vec fft2d  (const Vec& xx, size_t r, size_t c, bool inverse=false);
template<class Vec> Vec ifft2d (const Vec& ff, size_t r, size_t c);

template<class Vec> Vec powersp(const Vec& xx);                      // power spectrum
template<class Vec> Vec lowpass(const Vec& xx, unsigned f_cut);      // lowpass filter
template<class Vec> Vec convolv(const Vec& xx, const Vec& yy);       // convolution

template<class Vec> Vec ft_slow (const Vec& xx, bool inverse=false);
template<class Vec> Vec ift_slow(const Vec& xx);


/******** Implementation ********/

inline unsigned shuffle(unsigned i,unsigned n){
  unsigned j=0;
  for(n>>=1; n; i>>=1,n>>=1) j=(j<<1)|(i&1);  // bit reversal (0..n-1)
  return j;
}
template<class Cpx> void butterfly(Cpx& x1, Cpx& x2, const Cpx& w){
  Cpx x3=w*x2;
  x2 =x1-x3;  // x1-w*x2
  x1+=   x3;  // x1+w*x2
}
template<class Vec> Vec fft(const Vec& xx0, bool inverse){
  typedef typename Vec::value_type Cpx;
  unsigned n=xx0.size(), i, j, p, q;
  Vec xx(xx0);

  for(i=1; i<n; i<<=1)
    if(n&i) throw std::domain_error("fft: size is not a power of 2");

  for(i=0; i<n; i++)
    if(i<shuffle(i,n)) std::swap(xx[i],xx[shuffle(i,n)]);

  for(p=n/2,q=2; p; p>>=1,q<<=1){
    Cpx w1=exp(Cpx(0,-2*M_PI/q)), wj;
    if(inverse) w1=conj(w1);

    for(i=0; i<p; i++)
      for(j=0,wj=1; j<q/2; j++,wj*=w1)
	butterfly(xx[i*q+j], xx[i*q+q/2+j], wj);
  }

  if(inverse)
    for(i=0; i<n; i++) xx[i]/=n;
  return xx;
}


template<class Vec> Vec fft2d(const Vec& xx, size_t r, size_t c, bool finv){
  Vec yy(r*c),y1(c);
  for(size_t i=0; i<r; i++){
    for(size_t j=0; j<c; j++) y1[j]=xx[i*c+j];
    y1=fft(y1,finv);
    for(size_t j=0; j<c; j++) yy[i*c+j]=y1[j];
  }
  Vec zz(r*c),z1(r);
  for(size_t j=0; j<c; j++){
    for(size_t i=0; i<r; i++) z1[i]=yy[i*c+j];
    z1=fft(z1,finv);
    for(size_t i=0; i<r; i++) zz[i*c+j]=z1[i];
  }
  return zz;
}


template<class Vec> Vec ft_slow(const Vec& xx0, bool inverse){
  typedef typename Vec::value_type Cpx;
  unsigned n=xx0.size(), i, j;
  Vec xx(n);

  for(i=0; i<n; i++){
    Cpx w1=exp(Cpx(0,-2*M_PI*i/n)), wj;
    if(inverse) w1=conj(w1);
    for(j=0,wj=1; j<n; j++,wj*=w1) xx[i]+=wj*xx0[j];
  }

  if(inverse)
    for(i=0; i<n; i++) xx[i]/=n;
  return xx;
}


template<class Vec> Vec ifft    (const Vec& ff){  return fft    (ff,true);  }
template<class Vec> Vec ift_slow(const Vec& ff){  return ft_slow(ff,true);  }
template<class Vec> Vec ifft2d  (const Vec& ff, size_t r, size_t c){ return fft2d(ff,r,c,true); }


template<class Vec> Vec powersp(const Vec& xx){
  Vec ff=fft(xx);
  for(unsigned i=0; i<ff.size(); i++) ff[i]=norm(ff[i]);  // note: norm(C)=|C|^2 (!=|C|)
  return ff;
}
template<class Vec> Vec lowpass(const Vec& xx, unsigned f_cut){
  if(f_cut==0) return Vec(xx.size());
  Vec ff=fft(xx);
  for(unsigned i=f_cut; i<=ff.size()/2; i++) ff[i]=ff[ff.size()-i]=0;
  return ifft(ff);
}
template<class Vec> Vec convolv(const Vec& xx, const Vec& yy){
  if(xx.size()!=yy.size())
    throw std::domain_error("convolv: sizes mismatch");
  Vec ff=fft(xx), gg=fft(yy);
  for(unsigned i=0; i<ff.size(); i++) ff[i]*=gg[i];
  return ifft(ff);
}


}  //namespace nmlib
#endif  //FFT_H
