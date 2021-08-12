//header file for tree code

#ifndef TREE_CODE_H
#define TREE_CODE_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <string>
#include <new>
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace std;
//first define some necessary types.

class particle {

 public:

  vector<float> pos;
  //float x, y, z;
  //particle(float x_, float y_, float z_) {x=x_; y=y_; z=z_;}
  particle() {};
  particle(vector<float> pos_) {pos=pos_;} //0 for x, 1 for y, 2 for z 
  //but in principle the code can be used for 2D data also
  int dim_data(){ return pos.size(); } //what is the dimension of the objects to be correlated
 
};


class tree_node {
  
 public:
  //made cm a vector to avoid giving it a dimension, and for quicker initialization
  tree_node(){};
  tree_node(int dnum_, int level_, int l_, int u_, float val_, vector<float> cm_, float rad_) { 
    dnum=dnum_; level=level_, l=l_; u=u_; val=val_; cm=cm_; rad=rad_; 
  }
  int dnum, level; //dimension to cut// depth of the tree
  int l, u; //indices of particles in this node (lower and upper)
  float val, rad; //where to cut the box side, radius of node wrt to CM
  vector<float> cm; // position CM
  tree_node *left, *right; //child pointers 
  void set_tree_node(int dnum_, int level_, int l_, int u_, float val_, 
		     vector<float> cm_, float rad_, tree_node *pl, tree_node *pr) {
    dnum=dnum_; level=level_, l=l_; u=u_; val=val_; cm=cm_; rad=rad_; left=pl; right=pr;
    //cerr<<"number of particles in this node with l "<<l<<", u "<<u<<", level "<<level<<" is "<<u-l+1<<endl;
  }
  ~tree_node() { if(left != 0) left=0; if(right != 0) right=0; }
};




//combine Robert's create_tree, build_tree, build_tree_for range and tree_master_record

class build_tree { 

 private:

  static const int bucket_size=50;
  static const float MinNodeSize=1.; //Mpc/h
  const double rbox; //size of simulation box
  const int dim; //number of spatial dimensions
  const size_t Npart;
  vector<particle> &Vpart; //particle data vector as reference
  int *indexes;
  //tree_node *root; //pointer to the root of the tree; WHAT DO WE DO HERE?????
  int level;
  int count_levels;
  vector<long> count_buckets;
  //functions

  void select_on_coord(int c, int k, int l_arg, int u_arg);
  float spread_in_coord(int c, int l, int u); 
  int most_spread_coord(int l, int u);
  void NodePhenom(vector<float> &cm, float &rad, int l, int u);

 public:
  tree_node *root;
  build_tree(const double rbox_, vector<particle> &Vpart_);
  tree_node* build_tree_for_range(const int l, const int u, const int lev);
  ~build_tree();
  friend class correlations;
  friend class projected_correlations;

};


  
#endif 

