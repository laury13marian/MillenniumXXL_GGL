#ifndef TREE_PROCESSING_H
#define TREE_PROCESSING_H

#include "tree_code.h"


class correlations {

 private:
  

  build_tree *Ptree1, *Ptree2; 
  const int NBIN;
  const float rmin, rmax;
  const vector<float> vecPER; //to be initialized in the constructor
  static const int iPER=0; //same
  //int iPER;
  //static const float amin=0.;
  int dim;
  size_t Npart1, Npart2;
  double rbox;
  bool FLAG_PER;
  vector<float> rbinsL, rbinsU, rbinsM; //lower, upper, middle 
  float delta_bin; //log, log10, lin
  //float maxDIST, minDIST, dist, amin;
  float amin;
  double weight;
  //functions

  //double range_count(tree_node *node, vector<float> qv, float rHI, vector<float> vecP, int iP);
  inline void dual_tree_count(tree_node *node1, tree_node *node2, float rLO, float rHI, double &num);
  inline void dual_tree_count_all_bins(tree_node *node1, tree_node *node2, vector<double> &num);
  inline void dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI, 
				 vector<float> vecP, int iP, double &num);
  inline void dual_tree_count_PB_all_bins(tree_node *node1, tree_node *node2, 
					  vector<float> vecP, int iP, vector<double> &num);
  inline void get_pairs_dual_tree_all_bins(tree_node *node1, tree_node *node2, vector<float> vecP,
					   int iP, float minDIST, float maxDIST, vector<double> &num);
  inline void get_pairs_dual_tree(tree_node *node1, tree_node *node2, float rLO, float rHI, vector<float> vecP,
				  int iP, float minDIST, float maxDIST, double &num);

  void set_bins();  

 public:


  enum BinType { LOG, LOG10, LIN };
  BinType bin_type;
  correlations(build_tree *Ptree1_, build_tree *Ptree2_, correlations::BinType bt_, 
	       const float rmin_, const float rmax_, const int NBIN_, bool flag_per_, 
	       const vector<float> VP_); 
  vector<double> dual_tree_2p(string ss, bool VERBOSE=false);
  vector<double> dual_tree_2p_all_bins(string ss);
  vector<double> brute_force_XX(string ss);
  vector<double> brute_force_XY();
  vector<float> bins_return();
  vector<double> generate_RR(bool write_file, const string file_name);  
  vector<double> volume_bins();
 
};

#endif
