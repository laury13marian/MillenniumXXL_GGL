//code to compute statistics for the correlation functions

#ifndef JK_STATISTICS_H
#define JK_STATISTICS_H
#include "general.h"
#include "/u/lmarian/ZURICH_CODES/Cosmology/utilities/Matrix.h"


class JK_statistics {

 private:

  static const int snap=54;
  static const int NR=30;
  static const int ntasks=32;
  static const int NM=4; // only 4 magnitude bins are actually used here
  static const int NJK=64;
  const string path_gen;
  const int sub; //subcube that has been jack-knifed, for me number 100
  vector<float> binsR;
  vector<double> mags;
  vector<double> Vgg, SHgm; //initial vectors
  vector<double> Agg, Ashgm; //shear average
  ostringstream osSNAP, osNR, osNM, osT, osSC;

 public:

  JK_statistics(const string path_gen_, vector<double> mags_, const int sub_, int TASK=1000);
  void read_data(int JKsample, int TASK);
  void signal_average();
  vector<SqDMatrix> signal_covariance(string type);
  void write_covariance_matrix(vector<string> file, string type);
  void write_correlation_matrix(vector<string> file, string type);


};

#endif
