// Usage of FFT code
// - Include "fft.h" and use "nmlib" namespace. No linking is required.
// - Features include FFT, inverse FFT and their by-products.


#include <iostream>
#include <vector>
#include <complex>
#include "fft.h"
using namespace nmlib;


int main(void){
  size_t n=1024, ncut=20;
  std::vector<std::complex<double> > x(n),y, z;

  for(size_t k=0; k<n; k++) x[k]=k%64;  // sample wave - triangle

  y= fft(x);  // y[f] = \sum_k exp(-2\pi ikf/n) x[k]
  z=ifft(y);  // x[k] = \sum_f exp( 2\pi ikf/n) y[f] / n
  y=powersp(x);
  z=lowpass(x,ncut);  // cut ncut-th frequency (ncut/T) or higher
  //...

  for(size_t k=0; k<n; k++)
    std::cout << k << "\t" << x[k].real() << "\t" << y[k].real() << "\t" << z[k].real() << std::endl;

  return 0;
}
