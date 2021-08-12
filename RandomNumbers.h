// C++ implementation of the ran2.c routine in NR

#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

#include <cmath>
#include <cstdlib>
using namespace std;

class random2 {

 private:

  long unsigned IM1, IM2, IMM1;
  long IA1, IA2, IQ1, IQ2, IR1, IR2;
  double AM, EPS, RNMX, NDIV;
  const static int NTAB=32;
  long idum, idum2;
  long iy;
  long iv[NTAB];

 public:

  random2(const long lseed) {
    IM1=2147483563;
    IM2=2147483399;
    AM=(1.0/IM1);
    IMM1=(IM1-1);
    IA1=40014;
    IA2=40692;
    IQ1=53668;
    IQ2=52774;
    IR1=12211;
    IR2=3791;
    NDIV=(1+IMM1/NTAB);
    EPS=1.2e-7;
    RNMX=(1.0-EPS);

    initiate(lseed);
  }


  void initiate(long seed) {
    long k;
    idum=seed;
    if (idum <= 0) {
    if (-idum < 1) idum=1;
    else idum = -idum;
    idum2=idum;
    for (int j=NTAB+7; j>=0; j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] =  idum;
    }
    iy=iv[0];
    }
    
  }

  float newvalue() {

    long k;
    int j; 
    float temp;

    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
  }

  float operator() () {return newvalue();} 

};

#endif  
