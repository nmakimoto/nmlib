// Random number generator


#include <random>
#include "random.H"
namespace nmlib{


static std::mt19937_64 MT;
static std::uniform_real_distribution<double> URAND(0.0,1.0);
static std::normal_distribution      <double> NRAND(0.0,1.0);
static std::exponential_distribution <double> ERAND(1.0);

void init_rand(unsigned long long seed){
  if(seed==0) seed=std::random_device()();
  MT=std::mt19937_64(seed);
}

unsigned long long irand(void){ return MT(); }
double urand(void){ return URAND(MT); }
double nrand(void){ return NRAND(MT); }
double erand(void){ return ERAND(MT); }


}  //namespace
