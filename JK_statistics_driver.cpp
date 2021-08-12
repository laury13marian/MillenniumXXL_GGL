#include "JK_statistics.h"

using namespace std;

int main(int argc, char *argv[])

{

  if(argc!=2) {
    cerr<<"task number"<<endl;
    exit(10);
  }

  const int TASK=atoi(argv[2]);
  const int Nmags=5;
  const int NM = Nmags-1;
  const int NJK = 64;
  const int SUB = 100; //this is the subcube that we jack-knifed
  ostringstream ossNJK;
  ossNJK<<NJK;
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.;
 
  string TYPE = "gggg";
  vector<string> file_cov, file_corr;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/JACK_KNIFE/";
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/STATISTICS/JACK_KNIFE/SUBCUBE_100/";

  
  for(int im=0; im<NM; im++) {
    ostringstream oss;
    oss<<im;
    string fcov, fcorr;
    if(TASK==1000) {
      fcov = path_out + "covariance_matrix_" + TYPE + "_NJK_" + ossNJK.str() + 
	"_magnitude_bin_" + oss.str() + ".dat";
      fcorr = path_out + "correlation_matrix_" + TYPE + "_NJK_" + ossNJK.str() + 
	"_magnitude_bin_" + oss.str() + ".dat";
    }
    else {
      fcov = path_out + "ONE_TASK/covariance_matrix_" + TYPE + "_JK_block_" + argv[1] + 
	"_magnitude_bin_" + oss.str() + "_one_task_" + argv[2] + ".dat";
      fcorr = path_out + "ONE_TASK/correlation_matrix_" + TYPE + "_JK_block_" + argv[1] + 
	"_magnitude_bin_" + oss.str() + "_one_task_" + argv[2] + ".dat";
    }
    
    file_cov.push_back(fcov);
    file_corr.push_back(fcorr);
  }
  
  JK_statistis OB(path_gen, mags, SUB, TASK);
  OB.signal_average();
  //OB.write_covariance_matrix(file_cov, TYPE);
  //OB.write_correlation_matrix(file_corr, TYPE);


  return 0;

}

