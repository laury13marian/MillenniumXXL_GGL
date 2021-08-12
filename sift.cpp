#include "general.h"

using namespace std;

int main()

{

  vector<long long> VN;
  vector<int> Vsub;
  const int NN=216;
  string file = "/u/lmarian/LauraTreeCode/LauraMPI/Log_check_54_v2";
  ifstream in(file.c_str(), ios::in);
  check_stream(in, file);
  for(int i=0; i<NN; i++) {
    string s[6]; int i; long long l;
    in>>s[0]>>s[1]>>s[2]>>s[3]>>s[4]>>s[5]>>i>>l;
    if(l==0 && i==0)
      continue;
    VN.push_back(l); Vsub.push_back(i);
  }
  in.close();
  size_t sum=0;
  cerr<<VN.size()<<' '<<Vsub.size()<<endl;
  for(int i=0; i<VN.size(); i++) {
    cerr<<Vsub[i]<<' '<<VN[i]<<endl;
    sum+=VN[i];
  }

  size_t orig_number=303464448000;
  cerr<<"total number of particles: "<<sum<<" original number: "<<orig_number<<" difference: "<<sum-orig_number<<endl;


  return 0;

}
