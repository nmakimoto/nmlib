// Usage of Spline
// - Include "spline.h" and use namespace "nmlib". No linking is required.
// - Features: interpolation(splint/linear/nearest)


#include <iostream>
#include <map>
#include "spline.h"
using namespace nmlib;


int main(void){
  std::map<double,double> x2y;  // control points (x0,y0)...(xn,yn)
  x2y[10]=4;
  x2y[20]=7;
  x2y[35]=3;
  x2y[40]=5;
  //...

  Spline f(x2y);  // specify control points
  //Spline f;
  //f.set(x2y);  // specify control points later
  //f.get(x2y);  // retrieve control points later

  for(double x=0; x<50.1; x+=1){
    double y,y1,y0,dy;
    y =f(x);       // deg=default(=3): natural spline
    y1=f(x,1);     // deg=1: peicewise linear
    y0=f(x,0);     // deg=0: nearest neighbour
    dy=f.grad(x);  // d/dx of spline
    std::cout<<x<<' '<<y<<' '<<y1<<' '<<y0<<' '<<dy<<'\n';
  }

  return 0;
}
