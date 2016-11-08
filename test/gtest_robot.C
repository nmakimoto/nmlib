// Unit Test (robot)


#include <gtest/gtest.h>
#include <stdexcept>

#include <fstream>
#include "robot.H"
#include "io.H"
using namespace nmlib;


// common data
inline void init_robot(Robot& rbt, Matrix& th, Matrix& hom){
  std::ifstream str("gtest_robot.dat");
  str >> rbt   // Yaskawa Motoman ES280D-230 (6-axis manipulator)
      >> th    // initial joint angles [deg] (6-vector)
      >> hom;  // target pose of end tool (homogeneous transformation matrix)
  th*=M_PI/180;
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

  for(int k=0; k<6; k++){
    // initial position of end tool relative to k-th part
    Matrix hom00,pos0;
    hom00 = (k>0 ? rbt.fk(th0,k-1) : Matrix(4,4)+1.0);
    pos0 = getsub(inv_homtrsf(rbt.fk(th0,k))*rbt.fk(th0), 0,3,3,1);

    for(int i=0; i<10; i++){
      // rotation around k-th joint
      Matrix pos1, th1;
      th1  = th0;
      th1(k)+=2*M_PI*(i+1.0)/10;
      pos1 = getsub(inv_homtrsf(rbt.fk(th0,k))*rbt.fk(th1), 0,3,3,1);
      EXPECT_TRUE(is_homtrsf(rbt.fk(th1,k)));
      EXPECT_NEAR(norm(rbt.fk(th1)-rbt.fk(th1,6)), 0, 1.e-8);

      // must be a rotation around Z-axis of k-th part
      double t=2*M_PI*(i+1.0)/10, c=cos(t), s=sin(t);
      EXPECT_NEAR(pos1(0)-(c*pos0(0)-s*pos0(1)), 0, 1.e-8);
      EXPECT_NEAR(pos1(1)-(s*pos0(0)+c*pos0(1)), 0, 1.e-8);
      EXPECT_NEAR(pos1(2)-pos0(2), 0, 1.e-8);
    }
  }
}


TEST(robot,jacobian){
  Robot rbt;
  Matrix th0, hom;
  init_robot(rbt,th0,hom);

  for(int k=0; k<6+1; k++){
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

      if(k==6) EXPECT_NEAR(norm(rbt.jacobian(th0)-rbt.jacobian(th0,k)), 0, 1.e-8);  // k=6 by default
    }
  }
}


TEST(robot,ik){
  Robot rbt;
  Matrix th0, pos0, rot0, hom0, pos1, rot1, hom1;
  init_robot(rbt,th0,hom1);

  // start and goal
  hom0=rbt.fk(th0);
  pos0=hom2pos(hom0);
  rot0=hom2rot(hom0);
  pos1=hom2pos(hom1);
  rot1=hom2rot(hom1);

  Matrix th=th0, pos, rot, hom;
  int n=10;
  for(int i=0; i<n; i++){
    // set next checkpoint by interpolation
    double a=(i+1.0)/n;
    pos = a*pos1 + (1-a)*pos0;
    rot = vec2rot(a*rot2vec(rot1) + (1-a)*rot2vec(rot0));
    hom = homtrsf(pos,rot);

    // checkpoint
    th = rbt.ik(hom,th);
    EXPECT_NEAR(norm(rbt.fk(th)-hom), 0, 1.e-8);
  }

  // goal
  EXPECT_NEAR(norm(rbt.fk(th)-hom1), 0, 1.e-8);
}
