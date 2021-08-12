// code to compute statistics for the correlation functions

#ifndef JACK_KNIFE_H
#define JACK_KNIFE_H
#include "general.h"
#include "/u/lmarian/ZURICH_CODES/Cosmology/utilities/Matrix.h"


class jack_knife {

 private:

  static const int Nsub=216;
  static const int snap=54;
  static const int NZ=10;
  static const int NR=30;
  static const int nsampM=4000;
  static const int ntasks=32;
  static const long long NRAN=1000000;
  static const float zmax=100.;
  static const int NM2=4; // only 4 magnitude bins are actually used here
  const string path_gen;
  const int block_size;
 
  int NM, NJK;
  vector<float> binsR;
  vector<double> mags;
  vector<double> Vgg, SHgm, VJKgg, VJKshgm; //shear vectors
  vector<double> Agg, Ashgm; //shear average
  ostringstream osSNAP, osNZ, osNR, osNM, osNSAMP, osT, osNRAN, osZMAX;

 public:

  jack_knife(const string path_gen_, vector<double> mags_, const int block_size_, int TASK=1000);
  void read_data(int sub, int TASK);
  void jack_knife_data(int ijack);
  void jack_knife_data();
  void signal_average();
  vector<SqDMatrix> signal_covariance(string type);
  void write_covariance_matrix(vector<string> file, string type);
  void write_correlation_matrix(vector<string> file, string type);


};

#endif
