// Usage of Chrono class
// - Include "chrono.H" and use namespace "nmlib". Nothing to compile/link.
// - Features: time measurement (resolution=1.0[usec])


#include "chrono.H"
using namespace nmlib;


int main(void){
  Chrono sw;  // new object - lap() is reset() to zero
  double t;

  //...task...

  t=sw.lap();  // time[sec] since last reset()
  sw.reset();  // reset lap() to zero for subsequent measurements

  return 0;
}
