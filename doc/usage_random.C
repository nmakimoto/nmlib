// Usage of random number generators
// - Include "random.H" and use namespace "nmlib".
// - Compile and link "random.C" if you use *rand() functions. Nothing to link otherwise.
// - Features: RNG(random number generator), LDS(low-discrepancy sequence)


#include <iostream>
#include "random.H"  // *rand(), lsd()
using namespace nmlib;


int main(void){
  init_rand(12345);  // initialize RNG with seed

  for(int seq=0; seq<100; seq++){
    unsigned long long i;
    double xu,xn,xe,xl;

    i =irand();  // [0,2^64-1] uniform RNG
    xu=urand();  // U(0,1) uniform
    xn=nrand();  // N(0,1) gaussian
    xe=erand();  // Ex(1)  exponential

    xl=lds(seq,2);  // LDS (base=2 vander Corput)

    std::cout<<i<<' '<<xu<<' '<<xn<<' '<<xe<<' '<<xl<<'\n';
  }
  return 0;
}
