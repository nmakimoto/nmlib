// Unit Test (robot)


#include <gtest/gtest.h>
#include <stdexcept>

#include <fstream>
#include "robot.h"
using namespace nmlib;


// sample data
static void init_robot(Robot& rbt, Matrix& th, Matrix& hom){
  Matrix y90=rotabout(1,M_PI/2), z90=rotabout(2,M_PI/2);

  rbt = Robot();
  rbt.hom0(0) = homtrsf(Matrix({1000,2000,1000}), z90);
  rbt.hom0(7) = homtrsf(Matrix({   0,1000,   0}), z90*y90);

  th = vcat(Matrix({20,-20,20}),Matrix({-20,20,-20})) * (M_PI/180);

  hom = homtrsf(Matrix({1000,3000,5000}), tp(z90));
}


TEST(robot,homtrsf){
  Robot rbt;
  Matrix th, hom, hom2;
  init_robot(rbt,th,hom);

  EXPECT_NEAR(hom(3,3),1, 1.e-8);
  EXPECT_NEAR(norm(getsub(hom,3,0,1,3)), 0, 1.e-8);
  EXPECT_NEAR(norm(tp(getsub(hom,0,0,3,3))*getsub(hom,0,0,3,3)-1.0), 0, 1.e-8);

  hom2=hom; hom2(0,0)+=0.1; EXPECT_FALSE(is_homtrsf(hom2));
  hom2=hom; hom2(3,3)+=0.1; EXPECT_FALSE(is_homtrsf(hom2));
  hom2=hom; hom2(3,0)+=0.1; EXPECT_FALSE(is_homtrsf(hom2));
  hom2=hom; hom2(0,3)+=0.1; EXPECT_TRUE(is_homtrsf(hom2));

  EXPECT_TRUE(is_homtrsf(hom));
  EXPECT_TRUE(is_homtrsf(inv_homtrsf(hom)));
  EXPECT_NEAR(norm(inv_homtrsf(hom)*hom-1.0), 0, 1.e-8);
  EXPECT_NEAR(norm(inv_homtrsf(inv_homtrsf(hom))-hom), 0, 1.e-8);

  EXPECT_NEAR(norm(homtrsf(hom2pos(hom),hom2rot(hom))-hom), 0, 1.e-8);
}


TEST(robot,fk){
  Robot rbt;
  Matrix th0, hom;
  init_robot(rbt,th0,hom);

  EXPECT_THROW(rbt.fk(th0,-1), std::domain_error);
  EXPECT_THROW(rbt.fk(th0, 8), std::domain_error);

  for(int j=0; j<6; j++){
    // initial position of end tool relative to (j+1)-st part
    Matrix pos0 = getsub(inv_homtrsf(rbt.fk(th0,j+1))*rbt.fk(th0), 0,3,3,1);

    for(int i=0; i<10; i++){
      // motion of end tool relative to the initial position of (j+1)-st part...
      Matrix pos1, th1;
      th1  = th0;
      th1(j)+=2*M_PI*(i+1.0)/10;
      pos1 = getsub(inv_homtrsf(rbt.fk(th0,j+1))*rbt.fk(th1), 0,3,3,1);
      EXPECT_TRUE(is_homtrsf(rbt.fk(th1,j+1)));
      EXPECT_NEAR(norm(rbt.fk(th1)-rbt.fk(th1,7)), 0, 1.e-8);

      // will be a rotation around Z-axis of (j+1)-st part
      double t=2*M_PI*(i+1.0)/10, c=cos(t), s=sin(t);
      EXPECT_NEAR(pos1(0)-(c*pos0(0)-s*pos0(1)), 0, 1.e-8);
      EXPECT_NEAR(pos1(1)-(s*pos0(0)+c*pos0(1)), 0, 1.e-8);
      EXPECT_NEAR(pos1(2)-pos0(2), 0, 1.e-8);

      // parts 0..j will not move
      EXPECT_NEAR(norm(rbt.fk(th1,j)-rbt.fk(th0,j)), 0, 1.e-8);
      EXPECT_NEAR(norm(rbt.fk(th1,0)-rbt.fk(th0,0)), 0, 1.e-8);
    }
  }
}


TEST(robot,jacobian){
  Robot rbt;
  Matrix th0, hom;
  init_robot(rbt,th0,hom);

  EXPECT_THROW(rbt.jacobian(th0,-1), std::domain_error);
  EXPECT_THROW(rbt.jacobian(th0, 8), std::domain_error);
  EXPECT_NEAR(norm(rbt.jacobian(th0,0)), 0, 1.e-8);  // the world will not move

  for(int k=0; k<6+2; k++){
    for(int j=0; j<6; j++){
      Matrix dth, df1, df2;
      dth=Matrix(6);
      dth(j)=1.e-8;

      // upper 3 rows: compare velocity with difference quotients
      df1 = (hom2pos(rbt.fk(th0+dth/2.0,k)) - hom2pos(rbt.fk(th0-dth/2.0,k))) / dth(j);
      df2 = getsub(rbt.jacobian(th0,k), 0,j,3,1);
      if(1.e-8 < norm(df1)) EXPECT_NEAR(norm(df2-df1), 0, norm(df1)*1.e-4);
      if(1.e-8 < norm(df2)) EXPECT_NEAR(norm(df2-df1), 0, norm(df2)*1.e-4);

      // lower 3 rows: compare angular velocity with difference quotients
      df1 = rot2vec(hom2rot(rbt.fk(th0+dth/2.0,k)) * tp(hom2rot(rbt.fk(th0-dth/2.0,k)))) / dth(j);
      df2 = getsub(rbt.jacobian(th0,k), 3,j,3,1);
      if(1.e-8 < norm(df1)) EXPECT_NEAR(norm(df2-df1), 0, norm(df1)*1.e-4);
      if(1.e-8 < norm(df2)) EXPECT_NEAR(norm(df2-df1), 0, norm(df2)*1.e-4);

      if(k==7) EXPECT_NEAR(norm(rbt.jacobian(th0)-rbt.jacobian(th0,k)), 0, 1.e-8);  // k=7 by default
    }
  }
}


TEST(robot,ik){
  Robot rbt;
  Matrix th0, pos0, rot0, hom0, th1, pos1, rot1, hom1;
  init_robot(rbt,th0,hom1);

  // start and goal
  hom0=rbt.fk(th0);
  pos0=hom2pos(hom0);
  rot0=hom2rot(hom0);
  th1 =th0 + vcat(Matrix({11,12,13}),Matrix({14,15,16})) * (M_PI/180);
  hom1=rbt.fk(th1);
  pos1=hom2pos(hom1);
  rot1=hom2rot(hom1);

  Matrix th=th0, pos, rot, hom;
  int n=10;
  for(int i=0; i<n; i++){
    // set next checkpoint by interpolation
    double a=(i+1.0)/n;
    pos = a*pos1 + (1-a)*pos0;
    rot = vec2rot(a*rot2vec(rot1*tp(rot0))) * rot0;
    hom = homtrsf(pos,rot);

    // checkpoint
    th = rbt.ik(hom,th);
    EXPECT_NEAR(norm(rbt.fk(th)-hom), 0, 1.e-8);
  }

  // goal
  EXPECT_NEAR(norm(rbt.fk(th)-hom1), 0, 1.e-8);
}
