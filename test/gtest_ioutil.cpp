// Unit Test (I/O utility)


#include <gtest/gtest.h>
#include <stdexcept>

#include <vector>
#include <map>
#include <sstream>
#include "ioutil.h"
using namespace nmlib;


class Dummy{
public:
  double r,c;
};
std::istream& operator>>(std::istream& str,       Dummy& x){ str >> x.r >> x.c; return str; }
std::ostream& operator<<(std::ostream& str, const Dummy& x){ str << x.r << '\t' << x.c << '\t'; return str; }


TEST(ioutil,vector){
  std::stringstream str1, str2;
  std::vector<double> xx;

  str1 << "3  +11.1  2.22e1  3.33e+1";
  str1 >> xx;
  EXPECT_EQ(xx.size(), 3U);
  for(size_t i=0; i<xx.size(); i++) EXPECT_NEAR(xx[i], (i+1)*11.1, 1.e-8);

  xx.push_back(44.4);
  str2 << xx;
  EXPECT_EQ(str2.str(),"4\n11.1\n22.2\n33.3\n44.4\n");
}


TEST(ioutil,map){
  std::stringstream str1, str2;
  std::map<int,double> xx;

  str1 << "3  11  +11.1  22  222e-1  33  3.33e+1";
  str1 >> xx;
  EXPECT_EQ(xx.size(), 3U);
  int k=0;
  for(std::map<int,double>::const_iterator i=xx.begin(); i!=xx.end(); i++){
    k++;
    EXPECT_EQ(i->first, k*11);
    EXPECT_NEAR(i->second, k*11.1, 1.e-8);
  }

  xx[44]=44.4;
  str2 << xx;
  EXPECT_EQ(str2.str(),"4\n11\t11.1\n22\t22.2\n33\t33.3\n44\t44.4\n");
}


TEST(ioutil,string){
  std::string s="12.3 45.6";
  double x;
  x=0;    str2any(s,x);        EXPECT_NEAR(x,12.3,1.e-8);
  x=0;  x=str2any<double>(s);  EXPECT_NEAR(x,12.3,1.e-8);

  char s1[80]="45.6";
  x=0;    str2any(s1,x);        EXPECT_NEAR(x,45.6,1.e-8);
  x=0;  x=str2any<double>(s1);  EXPECT_NEAR(x,45.6,1.e-8);

  s = "3 11 11.1 22 22.2 33 33.3";
  std::map<int,double> xx = str2any<std::map<int,double> >(s);
  EXPECT_EQ(xx.size(), 3U);
  int k=0;
  for(std::map<int,double>::const_iterator i=xx.begin(); i!=xx.end(); i++){
    k++;
    EXPECT_EQ(i->first, k*11);
    EXPECT_NEAR(i->second, k*11.1, 1.e-8);
  }

  EXPECT_EQ(any2str(xx), "3\n11\t11.1\n22\t22.2\n33\t33.3\n");

  s="!@#$";
  EXPECT_THROW(str2any(s,x), std::runtime_error);
  s = "3 11 11.1 22 22.2 33 ";
  EXPECT_THROW(str2any(s,xx), std::runtime_error);
  s = "3 11 11.1a 22 22.2 33 33.3";
  EXPECT_THROW(str2any(s,xx), std::runtime_error);
}


TEST(ioutil,userclass){
  std::string s="11 22";
  Dummy x;

  x=str2any<Dummy>(s);
  EXPECT_EQ(any2str(x), "11\t22\t");
  EXPECT_EQ(any2str(x,5).find("1.10000e+"), 0U);  // precision=5
  EXPECT_NE(any2str(x,3).find("\t2.200e+"), std::string::npos);
}
