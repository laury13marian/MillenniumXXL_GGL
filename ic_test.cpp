#include "read_data.h"
#include "tree_processing_projected.h"


using namespace std;

int main(int argc, char *argv[])

{
  if(argc!=2) {
    cerr<<"introduce task number"<<endl;
    exit(1);
  }

  const double rbox=500.; //500.;//3000.;
  const float rMAX=150.;//110.;//30.;
  const float rMIN=0.2;
  const float zMIN=0.;
  const float zMAX=50.;
  const int NB=30; 
  const int NZ=10;
  float rb[3][NB]; //radial bin vector
  float zb[3][NZ]; //z-bin vector
  const size_t NRAN=2000000;
  const signed long rng_seed_gen = -1381718569;
  float SIDE=0.; //chosen subcube for replication  
  vector<float> binsRU(NB), binsRL(NB), binsRM(NB), binsZM(NZ), binsZL(NZ), binsZU(NZ);
  ostringstream osNB, osNZ, osZMAX, osNRAN; 
  osNB<<NB; osNRAN<<NRAN;
  osNZ<<NZ; osZMAX<<zMAX;
  double RR[NB][NZ];
  const string path_RAN =
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/PROJECTED_CORR_FUNC/RANDOM_CATALOGUE/";
  const string file_RAN = path_RAN + "pair_counts_NRAN_" + osNRAN.str() + "_binsR_" + osNB.str()
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_64" + "_task_" + argv[1] + ".dat" ; 
  ifstream in(file_RAN.c_str(), ios::in);                                                                              
  check_stream(in, file_RAN);                                                                                          
  double Norm =0.;
  double NormR=0.;
  const double PI=4.*atan(1.);
  const double Vtot=pow(rbox,3);
  double Vol=0.;
  for(int ir=0; ir<NB; ir++) 
    for(int iz=0; iz<NZ; iz++) {                                                                             
      double d1, d2;                                                                                                            
      float r1, r2, r3, z1, z2, z3;                                                                                             
      in>>r1>>r2>>r3>>z1>>z2>>z3>>d1>>d2; 
      binsRL[ir]=r1; binsRM[ir]=r2; binsRU[ir]=r3;
      binsZL[iz]=z1; binsZM[iz]=z2; binsZU[iz]=z3;
      RR[ir][iz]=d1;
      NormR+=RR[ir][iz];
      Norm+=(RR[ir][iz] * PI * (pow(binsRU[ir], 2)-pow(binsRL[ir], 2)) * (binsZU[iz]-binsZL[iz]));
      Vol+=(PI * (pow(binsRU[ir], 2)-pow(binsRL[ir], 2)) * (binsZU[iz]-binsZL[iz]));
    }                                                                                                                           
  in.close();                                                                                                                 

  cerr<<NormR<<' '<<Norm<<' '<<Vol<<' '<<Vol/Vtot<<' '<<pow(Vol/Vtot, 2)<<endl;

  return 0;

}
