// Unit Test (route)


#include <gtest/gtest.h>
#include <stdexcept>

#include "route.h"
using namespace nmlib;


TEST(route,convention){
  Route route;
  Route::Path path;
  path={99,99};  EXPECT_EQ(route.length(path), 0);
  path={11,99};  EXPECT_EQ(route.length(path),-1);
}


TEST(route,dijkstra){
  Route route;

  // node 11-88
  for(int i=1; i<8; i++){
    for(int j=1; j<=8; j++) route.dist[i*10+j][(i+1)*10+j]=route.dist[(i+1)*10+j][i*10+j]=10;  // horizontal
    for(int j=1; j<=8; j++) route.dist[j*10+i][j*10+(i+1)]=route.dist[j*10+(i+1)][j*10+i]= 1;  // vertical
  }
  // node 91-99
  for(int k=91; k<99; k++) route.dist[k][k+1]=route.dist[k+1][k]=route.dist[k+1][k]=10;

  Route::Path path;
  path=route(11,88);  EXPECT_NEAR(route.length(path), 10*7+1*7, 1.e-10);
  path=route(91,99);  EXPECT_NEAR(route.length(path), 10*8    , 1.e-10);
  path=route(11,99);  EXPECT_EQ(path.size(), 0U);
}
