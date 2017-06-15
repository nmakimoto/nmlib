// Unit test (FFT)


#include <gtest/gtest.h>

#include <iostream>
#include <stdexcept>
#include <cstdlib>  // random
#include "fft.h"
#include "matrix.h"
using namespace nmlib;


typedef std::complex<double> Cpx;
typedef std::vector<Cpx>  VecC;
typedef matrix<Cpx> MatC;

static double urand1(void){ return (random()%RAND_MAX+0.5)/RAND_MAX; }
static VecC random_sample(size_t n){
  VecC vc(n);
  for(size_t k=0; k<n; k++){ vc[k]=Cpx(urand1()*2-1, urand1()*2-1); }
  return vc;
}


TEST(fft,core){
  size_t n=512;
  VecC xx=random_sample(n);
  EXPECT_NEAR(norm(MatC(ift_slow( ft_slow(xx)))-MatC(xx)), 0, 1.e-8);  // slowIFT slowFT  = ID
  EXPECT_NEAR(norm(MatC( ft_slow(ift_slow(xx)))-MatC(xx)), 0, 1.e-8);  // slowFT  slowIFT = ID
  EXPECT_NEAR(norm(MatC(ifft(fft(xx)))-MatC(xx)), 0, 1.e-8);  // IFFT  FFT = ID
  EXPECT_NEAR(norm(MatC(fft(ifft(xx)))-MatC(xx)), 0, 1.e-8);  //  FFT IFFT = ID
  EXPECT_NEAR(norm(MatC( fft(xx))-MatC( ft_slow(xx))), 0, 1.e-8);  //  FFT = slow FT
  EXPECT_NEAR(norm(MatC(ifft(xx))-MatC(ift_slow(xx))), 0, 1.e-8);  // IFFT = slow IFT
}


TEST(fft,size){
  VecC xx;
  for(size_t n=1; n<=2048; n*=2){
    xx=random_sample(n);
    EXPECT_NEAR(norm(MatC( fft(xx))-MatC( ft_slow(xx)))/sqrt(n), 0, 1.e-8);
    EXPECT_NEAR(norm(MatC(ifft(xx))-MatC(ift_slow(xx)))/sqrt(n), 0, 1.e-8);
    EXPECT_NEAR(norm(MatC(ifft( fft(xx)))-MatC(xx))/sqrt(n), 0, 1.e-8);
    EXPECT_NEAR(norm(MatC( fft(ifft(xx)))-MatC(xx))/sqrt(n), 0, 1.e-8);
  }
  xx=VecC();   EXPECT_EQ(fft(xx).size(), 0);
  xx=VecC(20); EXPECT_THROW(fft(xx), std::domain_error);
  xx=VecC(31); EXPECT_THROW(fft(xx), std::domain_error);
}


TEST(fft,twodim){
  size_t r=64,c=128;
  VecC xx=random_sample(r*c),yy(r*c);
  EXPECT_NEAR(norm(MatC(ifft2d(fft2d(xx,r,c),r,c))-MatC(xx)), 0, 1.e-8);
  EXPECT_NEAR(norm(MatC(fft2d(ifft2d(xx,r,c),r,c))-MatC(xx)), 0, 1.e-8);
}
