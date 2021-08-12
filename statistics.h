//code to compute statistics for the correlation functions

#ifndef STATISTICS_H
#define STATISTICS_H
#include "general.h"
#include "/u/lmarian/ZURICH_CODES/Cosmology/utilities/Matrix.h"

class statistics {

 private:

  static const int Nsub=216;
  static const int snap=54;
  static const int NZ=10;
  static const int NR=30;
  static const int nsampM=4000;
  static const int ntasks=32;
  static const long long NRAN=1000000;
  static const float zmax=100.;
  const string path_gen;
  const int FLAG_IC;
  int NM;
  vector<float> binsR;
  vector<double> mags;
  vector<double> Vmm, Vgg, Vgm, SHmm, SHgm; //shear vectors
  vector<double> Amm, Agg, Agm, ASHmm, ASHgm; //shear average
  ostringstream osSNAP, osNZ, osNR, osNM, osNSAMP, osT, osNRAN, osZMAX;

 public:

  statistics(const string path_gen_, const int ic_, vector<double> mags_, int TASK=1000);
  void read_data(int sub, int TASK);
  void average_functions();
  vector<double> average_bias_gg();
  vector<double> average_bias_gm();
  vector<double> correlation_coefficient();
  vector<SqDMatrix> cov_corr_function(string type, string type_return);
  void write_average_functions(string file, string type);
  void write_bias(string file[2]);
  void write_correlation_coefficient(string file);
  void write_cov_corr(vector<string> file, string type, string type_return);
 
};

#endif
