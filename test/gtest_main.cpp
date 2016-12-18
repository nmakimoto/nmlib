// Google test for nmlib
// Compile: g++ -o gtest_main gtest_*.C -I/usr/src/gtest -lpthread
// Usage: $0 --help


#include <gtest/gtest.h>
#include <stdexcept>
#include <iostream>


int main(int argc, char* argv[]){
  ::testing::InitGoogleTest(&argc,argv);
  if(argc>1){ std::cerr<<"E unknown options <"<<argv[1]<<"...> (ABORTED)\n"; return -1; }
  return RUN_ALL_TESTS();
}
