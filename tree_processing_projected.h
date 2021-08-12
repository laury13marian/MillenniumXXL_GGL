#ifndef TREE_PROCESSING_PROJECTED_H
#define TREE_PROCESSING_PROJECTED_H

#include "tree_code.h"


class projected_correlations {
  
 private:
  

  build_tree *Ptree1, *Ptree2; 
  const int NBIN, ZBIN;
  const float rmin, rmax, zmin, zmax;
  const vector<float> vecPER; //to be initialized in the constructor
  static const int iPER=0; //same
  int dim;
  size_t Npart1, Npart2;
  double rbox;
  bool FLAG_PER;
  vector<float> rbinsL, rbinsU, rbinsM; //lower, upper, middle r_perp bins
  vector<float> zbinsL, zbinsU, zbinsM; //r_parallel=z bins
  float delta_bin_rperp, delta_bin_z; //log, log10, lin
  float maxD_rperp, minD_rperp, maxD_z, minD_z, dist_rperp, dist_z, amin;
  double weight;
  //functions

  inline void dual_tree_count(tree_node *node1, tree_node *node2, float rLO, float rHI, 
			      float zLO, float zHI, double &num);
  inline void dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI, 
				 float zLO, float zHI, vector<float> vecP, int iP, double &num);
  inline void get_pairs_dual_tree(tree_node *node1, tree_node *node2, float rLO, float rHI,
				  float zLO, float zHI, vector<float> vecP, int iP, float minD_rperp,
				  float maxD_rperp, float minD_z, float maxD_z, double &num);
  void set_bins_rperp();  
  void set_bins_z();
 public:
  
  enum BinType { LOG, LOG10, LIN };
  BinType bin_type;
 
  projected_correlations(build_tree *Ptree1_, build_tree *Ptree2_, projected_correlations::BinType bt_, 
			 const float rmin_, const float rmax_, const int NBIN_, const float zmin_, 
			 const float zmax_, const int ZBIN_, bool flag_per_, const vector<float> VP_); 
  vector<double> dual_tree_2p(string ss);
  vector<double> brute_force_XX();
  vector<double> brute_force_XY();
  vector<float> return_bins_rperp();
  vector<float> return_bins_z();
  vector<double> generate_RR(bool write_file); 
  vector<double> volume_bins();
  double IntegralConstraint;
};

#endif
