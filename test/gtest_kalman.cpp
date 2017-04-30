// Unit Test (random number generator)


#include <gtest/gtest.h>
#include <stdexcept>

#include <map>
#include "kalman.h"
#include "matrix.h"
using namespace nmlib;


// target model: linear motion with constant speed
inline Matrix f(const Matrix& x, double dt){
  int d=x.dim()/2;
  Matrix p=getsub(x,0,0,d,1);  // position
  Matrix v=getsub(x,d,0,d,1);  // velocity
  return vcat(p+v*dt,v);
}

// sensor model: observe position
inline Matrix h(const Matrix& x){
  int d=x.dim()/2;
  Matrix p=getsub(x,0,0,d,1);  // position
  return p;
}


TEST(kalman,linear){
  Kalman kf;
  Matrix Q,R,p0,v0,x,y;
  double t,dt;

  // initialize tracker
  kf.t =0;
  kf.x =Matrix(6).fill(10.0);
  kf.dx=Matrix(6).fill(1.e-4);
  kf.P =Matrix(6,6)+10000.0;

  // initialize models
  Q =Matrix(6,6)+0.1*0.1;
  R =Matrix(3,3)+0.1*0.1;
  p0=Matrix({0,0,0});
  v0=Matrix({1,2,3});
  dt=2;

  // track target
  for(int i=0; i<=10; i++){
    t=i*dt;
    x=p0+t*v0;  // + gaussian noise;
    y=h(vcat(p0+t*v0,v0));
    kf.predict(t,f,Q);
    kf.update (y,h,R);
  }

  EXPECT_NEAR(norm(kf.x-vcat(p0+t*v0,v0)), 0, 1.e-2);
  EXPECT_NEAR(norm(kf.P), 0, 1.e-1);
}
