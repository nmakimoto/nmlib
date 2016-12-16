// Random number generator


#include <random>
#include "random.H"
namespace nmlib{


static std::mt19937_64 MT;
static std::uniform_real_distribution<double> URAND(0.0,1.0);
static std::normal_distribution      <double> NRAND(0.0,1.0);
static std::exponential_distribution <double> ERAND(1.0);

void init_rand(ulong seed){
  if(seed==0) seed=std::random_device()();
  MT=std::mt19937_64(seed);
}

ulong  irand(void){ return ulong(MT()); }
double urand(void){ return URAND(MT); }
double nrand(void){ return NRAND(MT); }
double erand(void){ return ERAND(MT); }


}  //namespace
