#include "general.h"

using namespace std;

int main(int argc, char *argv[])

{

 if(argc != 2) {
    cerr<<"input JK number"<<endl; //varying from 0 to Lsub/DeltaJK-1
    exit(10);
  }

 int JKsample=atoi(argv[1]);
 int NJK=10;
 int Jix, Jiy;
 //JKsample=Jix*NJK+Jiy
 Jiy=(JKsample-1) % NJK;
 Jix=(JKsample-1-Jiy)/NJK;

 cerr<<"ix="<<Jix<<" iy="<<Jiy<<endl;

 return 0;

}
