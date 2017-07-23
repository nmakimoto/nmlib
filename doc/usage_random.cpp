// Usage of random number generators
// - Include "random.h" and use namespace "nmlib".
// - Compile and link "random.C" if you use *rand() functions. Nothing to link otherwise.
// - Features: RNG(random number generator), LDS(low-discrepancy sequence)


#include <iostream>
#include "random.h"
using namespace nmlib;


int main(void){
  Rng rng(12345);  // seed is optional
  Lds lds;

  auto pdf = [](double x){ return sqrt(1-x*x); };  // example probability density function
  double pdf_max=1;                                // and its upper bound (used by rejection)

  for(int seq=0; seq<100; seq++){
    int dim=3;
    double x, a=0, b=1;
    Matrix y;

    x=rng();       // U(0,1) uniform
    x=rng.u(a,b);  // U(a,b) uniform
    x=rng.n(a,b);  // N(a,b^2) gaussian
    x=rng.e(a);    // Ex(a)  exponential
    x=rng.rand_pdf(pdf,pdf_max,a,b);  // sample of given PDF
    y=rng.u(a,b,dim);  // U(a,b)^dim uniform
    y=rng.n(a,b,dim);  // N(a,b^2)^dim gaussian
    //...some other methods...

    x=lds();      // LDS on (0,1)
    //...methods analogous to RNG...
  }
  return 0;
}
