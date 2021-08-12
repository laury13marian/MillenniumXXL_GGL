#include "statistics.h"

using namespace std;

int main(int argc, char *argv[])
//int main()

{
  if(argc!=2) {
    cerr<<"introduce TASK number"<<endl;
    exit(10);
  }
  int TASK=atoi(argv[1]);
  const int Nmags=10;
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;
  
  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/";
  const int ic=2; //0 for noIC; 1 for DD; 2 for W
  const string xx="IC_W"; //"noIC"; //"IC_W"; //"IC_DD";
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/STATISTICS/";

  //output file for the average correlation functions
  //string type_ave = "gm"; 
  //string file_wp = path_out + "ONE_TASK/ave_" + type_ave + '_' + xx + "_one_task_" + argv[1] + ".dat";
  //string file_wp = path_out + "ave_" + type_ave + '_' + xx + ".dat";
  
  /*string file_bias[2]; //output files for the average bias 
  file_bias[0] = path_out + "ave_bias_gg_" + xx + ".dat";
  file_bias[1] = path_out + "ave_bias_gm_" + xx + ".abs.dat";

  string file_coeff; //output for the correlation coefficient
  file_coeff = path_out + "correlation_coefficient_" + xx + ".dat";*/

  vector<string> file_cov;
  string type_return = "covariance";//"covariance";"correlation";
  string type = "ggsgm"; //"sgmsgm";
  //string type = "smmsmm";
  //string fs = path_out + "COVARIANCE_MATRIX/SHEAR_" + xx + "/covariance_matrix_" + type + '_' + xx + ".dat";
  //string fs = path_out + "COVARIANCE_MATRIX/SHEAR_" + xx + "/ONE_TASK/covariance_matrix_" + type + '_' + xx 
  //+ "_one_task_" + argv[1] + ".dat";
  //file_cov.push_back(fs);

  for(int im=0; im<Nmags-1; im++) {
    ostringstream oss;
    oss<<im;
    //string fs = path_out + "COVARIANCE_MATRIX/SHEAR_" + xx + "/ONE_TASK/V2/covariance_matrix_" + type + '_' + xx 
    //+ "_magnitude_bin_" + oss.str() + "_one_task_" + argv[1] + ".dat";
    string fs = path_out + "COVARIANCE_MATRIX/CROSS_GG_GGL_" + xx + "/ONE_TASK/covariance_matrix_" + type + '_' + xx 
      + "_magnitude_bin_" + oss.str() + "_one_task_" + argv[1] + ".dat"; 
    //string fs = path_out + "COVARIANCE_MATRIX/CROSS_GG_GGL_" + xx + "/covariance_matrix_" + type + '_' + xx
    //+ "_magnitude_bin_" + oss.str() + ".dat";
    //string fs = path_out + "COVARIANCE_MATRIX/CROSS_GG_GGL_" + xx + "/correlation_matrix_" + type + '_' + xx
    //+ "_magnitude_bin_" + oss.str() + ".dat";
    file_cov.push_back(fs);
    }
  statistics OB(path_gen, ic, mags, TASK);
  //OB.write_average_functions(file_wp, type_ave);
  //statistics OB(path_gen, ic, mags); 
  //OB.cov_corr_function(type, type_return);
  OB.write_cov_corr(file_cov, type, type_return);

  //OB.write_ave_corr(file_corr);
  //OB.write_bias(file_bias);
  //OB.write_correlation_coefficient(file_coeff);
  
  return 0;

}
