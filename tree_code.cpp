//code for tree correlation function

#include "tree_code.h"


//PRIVATE FUNCTIONS OF THE BUILD_TREE CLASS


float build_tree::spread_in_coord(int c, int l, int u) //c=dimension to be cut, c is 0,1,2 or 0,1
//l and u are both used as Vpart[l] and Vpart[u]--which means the elements with indexes l+1, u+1
{
  //my version of this function, shorter than Rob's, performs equally well
  vector<float> VTemp;
  for(int i=l; i<=u; i++)
    VTemp.push_back(Vpart[indexes[i]].pos[c]);
  sort(VTemp.begin(), VTemp.end());
  
  return (VTemp[VTemp.size()-1]-VTemp[0]);
  
}



void build_tree::select_on_coord(int c, int k, int l_arg, int u_arg) 

{

  //k is a boundary, below which all indexes will correspond to elements of Vpart
  //with coord in c smaller than Vpart[indexes[k]][c], and above will be bigger
  //but the elements are not ordered beyond that
  //again, for us V[l_arg] and V{u_arg] refer to l_arg+1 and u_arg+1 elements, so CHECK later 
  //this function changes the indexes array pointer !!!
  //cerr<<"inside select_on_coord for dim spread "<<c<<", boundary "<<k<<", lower "<<l_arg
  //<<", upper "<<u_arg<<endl;
  //cerr<<"boundary value: "<<Vpart[indexes[k]].pos[c]<<endl<<endl;
  int l=l_arg;
  int u=u_arg; //initialize u & l, which will be updated within the next loop

  while(l<u) { //used to have a while here
    int t=indexes[l];
    int m=l; 
    int s;
    for(int j=l+1; j<=u; j++) {
      if(Vpart[indexes[j]].pos[c]<Vpart[t].pos[c]) {
	m++; 
	s=indexes[m]; //this piece here seems to be redundant, but affects greatly the speed of the code
	indexes[m]=indexes[j];
	indexes[j]=s;
      }
    }
    s=indexes[l];
    indexes[l]=indexes[m];
    indexes[m]=s;
    if(m<=k) l=m+1;
    if(m>=k) u=m-1;
  }
  
}


int build_tree::most_spread_coord(int l, int u)

{

  float bsp=-1.;
  int res;
  for(int i=0; i<dim; i++) {
    float sp=spread_in_coord(i, l, u);
    if( sp > bsp ) {
      res=i;
      bsp=sp;
    }
  }

  return res;
  
}


void build_tree::NodePhenom(vector<float> &cm, float &rad, int l, int u)

{
  //we consider Vpart[l], Vpart[u], which refer to the elements l+1, u+1--CHECK later 
  //compute the center of mass
  if(cm.size()!=dim) {
    cerr<<"incompatibility between dimension of cm given as arg";
    cerr<<" in NodePhenom and space dimension of the build_tree class"<<endl;
    exit(1);
  }
  for(int i=0; i<dim; i++) {
    cm[i]=0.;
    for(int j=l; j<=u; j++)
      cm[i]+=Vpart[indexes[j]].pos[i];
    cm[i]/=(u-l+1.);
  }
  
  //compute the radius
  
  float radMAX=0.;
  for(int j=l; j<=u; j++) {
    float rr=0.;
    for(int i=0; i<dim; i++) 
      rr+=pow((Vpart[indexes[j]].pos[i]-cm[i]), 2);
    if( rr > radMAX ) radMAX=rr;
  }
  rad=sqrt(radMAX);
  
}



//PUBLIC FUNCTIONS OF THE BUILD_TREE CLASS


build_tree::build_tree(const double rbox_, vector<particle> &Vpart_):
  rbox(rbox_), Vpart(Vpart_), Npart(Vpart_.size()), dim(Vpart_[0].dim_data()) 

{
  cerr<<"number of particles building the tree: "<<Npart<<endl;
  cerr<<"dimensions of the particle data: "<<dim<<endl;
  if(Npart>0) {
  indexes = new int[Npart]; //do some exception handling here
  for(size_t i=0; i<Npart; i++) 
    indexes[i]=i;
  level=1;
  count_levels=level;
  root = new tree_node;
  root = build_tree_for_range(0, Npart-1, level);
  cerr<<"TOTAL NUMBER OF LEVELS: "<<count_levels<<endl;
  }
}



build_tree::~build_tree()

{

  if(Npart>0) {  
    delete[]indexes; cerr<<"deleting indexes of the tree"<<endl;
    delete root; cerr<<"deleting root of the tree"<<endl;
    cerr<<"TOTAL NUMBER OF LEVELS: "<<count_levels<<endl;
  }

}


tree_node* build_tree::build_tree_for_range(const int l, const int u, const int lev)

{

  if( l<0 || l>=Npart ) { //changed here from Rob's code, which had l<1,l>Npart
    cerr<<"illegal lower index for build_tree_for_range: "<<l<<endl;
    exit(1);
  }

  if( u<0 || u>=Npart ) {//same change
    cerr<<"illegal upper index for build_tree_for_range: "<<u<<endl;
    exit(1);
  }

  if( u<l ) {
    cerr<<"upper index smaller than lower one: upper "<<u<<" lower "<<l<<endl;
    exit(1);
  }

  int dnum; //dnum the coord to be split
  vector<float> cm(dim, 0.);
  double val; float rad(0.);
  tree_node *res;

  res=new tree_node; //not sure this is the best way to do this, 
  //maybe simple alloc via object would be better

  if( (u-l+1) <= bucket_size ) {
    NodePhenom(cm, rad, l, u); //this should compute cm and rad, and return them by reference
    dnum=0; val=0.; 
    tree_node *pL, *pR;
    pL=0; pR=0;
    res->set_tree_node(dnum, lev, l, u, val, cm, rad, pL, pR);
    count_buckets.push_back(u-l+1);
  }
  
  else {
    int c=most_spread_coord(l, u); //this settles the dimension to be cut
    int m=(l+u)/2; //this settles the boundary of the cut
    select_on_coord(c, m, l, u); //this modifies the indexes pointer array
    NodePhenom(cm, rad, l, u); //this should compute cm and rad, and return them by reference
    tree_node *pL, *pR;
   
    if(2.*rad < MinNodeSize) {
      dnum=0; val=0.; pL=0; pR=0;
      res->set_tree_node(dnum, lev, l, u, val, cm, rad, pL, pR);
      count_buckets.push_back(u-l+1);
    }
    else{
      dnum=c; val=Vpart[indexes[m]].pos[c]; 
      int level_left=lev+1; 
      int level_right=lev+1;
      count_levels++;
      pL=build_tree_for_range(l, m, level_left);
      pR=build_tree_for_range(m+1, u, level_right);
      res->set_tree_node(dnum, lev, l, u, val, cm, rad, pL, pR);
      //count_buckets.push_back(u-l+1);
    }
    
  }
  /*cerr<<"EXITING build_tree_for_range with args: "<<endl; 
  cerr<<"level: "<<lev<<endl;
  cerr<<"cofm: "<<cm[0]<<' '<<cm[1]<<' '<<cm[2]<<endl;
  cerr<<"rad: "<<rad<<endl;
  cerr<<"left & right: "<<l<<' '<<u<<endl;
  cerr<<"dnum & val: "<<dnum<<' '<<val<<endl;
  cerr<<"most spread coord: "<<dnum<<endl;
  cerr<<"particles in this level: "<<endl;
  for(int j=l; j<=u; j++) {
    for(int d=0; d<dim; d++)
      cerr<<Vpart[indexes[j]].pos[d]<<"     ";
    cerr<<endl;
    }*/

  return res;
  
}
