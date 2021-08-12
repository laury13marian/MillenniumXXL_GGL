#include "jack_knife.h"

using namespace std;

int main(int argc, char *argv[])

{

  if(argc!=3) {
    cerr<<"introduce JK block size and task number"<<endl;
    exit(10);
  }
  const int block_size=atoi(argv[1]);
  const int TASK=atoi(argv[2]);
  const int Nmags=10;
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;
  
  //const int block_size=8;
  //ostringstream oblock;
  //oblock<<block_size;
  string TYPE = "ggsgm";
  vector<string> file_cov, file_corr;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/";
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/STATISTICS/JACK_KNIFE/";

  
  for(int im=0; im<4; im++) {
    ostringstream oss;
    oss<<im;
    string fcov, fcorr;
    if(TASK==1000) {
      fcov = path_out + "covariance_matrix_" + TYPE + "_JK_block_" + argv[1] + 
	"_magnitude_bin_" + oss.str() + ".dat";
      fcorr = path_out + "correlation_matrix_" + TYPE + "_JK_block_" + argv[1] + 
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
  
  jack_knife OB(path_gen, mags, block_size);
  //OB.signal_average();
  OB.write_covariance_matrix(file_cov, TYPE);
  OB.write_correlation_matrix(file_corr, TYPE);


  return 0;

}

