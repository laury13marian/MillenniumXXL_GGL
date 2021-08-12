
#include "tree_processing.h"


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

correlations::correlations(build_tree *Ptree1_, build_tree *Ptree2_, correlations::BinType bt_,
			   const float rmin_, const float rmax_, const int NBIN_, 
			   bool flag_per_, const vector<float> VP_):
  Ptree1(Ptree1_), Ptree2(Ptree2_), rmin(rmin_), rmax(rmax_), NBIN(NBIN_),
  FLAG_PER(flag_per_), vecPER(VP_)

{

  amin=0.;
  bin_type=bt_; 
  set_bins();
  dim=Ptree1->dim;
  rbox=Ptree1->rbox;
  Npart1=Ptree1->Npart;
  Npart2=Ptree2->Npart;
  if(Npart1!=Npart2) weight = Npart1 * Npart2;
  else weight = Npart1*(Npart1-1.)/2.;
  cerr<<"weight is: "<<weight<<endl;
  //if pointers Ptree1 and Ptree2 point to different catalogues, must be careful here
  cerr<<"start correlation function code for tree 1 with "<<Npart1<<" particles."
      <<" and tree 2 with "<<Npart2<<endl;
  //Npart1 is used only for the single_tree calc, so 2 catalogues are no problem
  cerr<<"number of levels in tree 1 (root=level 1)"<<Ptree1->count_levels
      <<" and number of levels in tree 2: "<<Ptree2->count_levels<<endl;
  cerr<<"size count_buckets tree1:  "<<Ptree1->count_buckets.size()
      <<" and tree 2: "<<Ptree2->count_buckets.size()<<endl<<endl;
 
}


void correlations::set_bins()

{

  if(bin_type==LOG) {

    double dd=log(rmax/rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(exp(log(rmin) + dd*i));
      rbinsU.push_back(exp(log(rmin) + dd*(i+1)));
      rbinsM.push_back(exp(log(rmin) + dd*(i+0.5)));
    }
    delta_bin=dd;

  }
 
  else if(bin_type==LOG10) {
  
    double dd=log10(rmax/rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(pow(10., (log10(rmin) + dd*i)));
      rbinsU.push_back(pow(10., (log10(rmin) + dd*(i+1))));
      rbinsM.push_back(pow(10., (log10(rmin) + dd*(i+0.5))));
    }
    delta_bin=dd;
  }
  else if(bin_type==LIN) {

    double dd=(rmax-rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(rmin + dd*i);
      rbinsU.push_back(rmin + dd*(i+1));
      rbinsM.push_back(rmin + dd*(i+0.5));
    }
    delta_bin=dd;
  }
  else {
    cerr<<"wrong bin_type in set_bins()"<<endl;
    exit(1);
  }
  cerr<<"setting up the radial bins for the correlation function: "<<endl;
  //for(int i=0; i<NBIN; i++)
  //cerr<<i<<' '<<rbinsL[i]<<' '<<rbinsM[i]<<' '<<rbinsU[i]<<endl;

  cerr<<endl;
}


vector<float> correlations::bins_return() 

{

  vector<float> V;
  for(int i=0; i<NBIN; i++)
    V.push_back(rbinsL[i]);
  for(int i=0; i<NBIN; i++)
    V.push_back(rbinsU[i]);
  for(int i=0; i<NBIN; i++)
    V.push_back(rbinsM[i]);

  return V;

}



void correlations::dual_tree_count(tree_node *node1, tree_node *node2, float rLO, float rHI, double &num)
				     
{
 float dist, minDIST, maxDIST;
  //get the max and min distances between the 2 nodes
  dist=0.; //distance bet. the 2 nodes
  for(int i=0; i<dim; i++)
    dist+=pow(node1->cm[i]-node2->cm[i], 2);
  dist=sqrt(dist);
  //cerr<<"DISTANCE: "<<dist<<" COUNTS: "<<num<<endl;
  maxDIST = dist + node1->rad + node2->rad;
  minDIST = max(dist-node1->rad-node2->rad, amin);

   
  if( (minDIST > rHI) || (maxDIST < rLO) ) {
    //no counts from these nodes
    return; //exits the function
  }
  
  else if( (minDIST >= rLO) && (maxDIST < rHI) ) {
    //include both nodes in the count
    num+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
    return;
  }
  //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.;
	for(int j=0; j<dim; j++)
	  rr+=pow((Ptree1->Vpart[Ptree1->indexes[i1]].pos[j])-(Ptree2->Vpart[Ptree2->indexes[i2]].pos[j]), 2);
	rr=sqrt(rr);
	if( (rr < rHI) && (rr >= rLO) ) num+=1.;
      }
    return;
  }
  
  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count(node1, node2->left, rLO, rHI, num);
    dual_tree_count(node1, node2->right, rLO, rHI, num);
    return;
  }
  
  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count(node1->left, node2, rLO, rHI, num); 
    dual_tree_count(node1->right, node2, rLO, rHI, num);
    return;
  }
  else {
    //split both nodes
    dual_tree_count(node1->left, node2->left, rLO, rHI, num);
    dual_tree_count(node1->left, node2->right, rLO, rHI, num);
    dual_tree_count(node1->right, node2->left, rLO, rHI, num);
    dual_tree_count(node1->right, node2->right, rLO, rHI, num);
    return;
  }
}




/*void correlations::dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI, 
				      vector<float> vecP, int iP, double &num)
				    
{
  float dist, minDIST, maxDIST;
  //get the max and min distances between the 2 nodes
  dist=0.; //distance bet. the 2 nodes
 
  for(int i=0; i<dim; i++) 
    dist+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  dist=sqrt(dist);
  maxDIST = dist + node1->rad + node2->rad;
  minDIST = max(dist-node1->rad-node2->rad, amin);
  
  //if( FLAG_PER && (iP == 0) && (minDIST > rHI) && (maxDIST > rbox/2.) ) {
    if( FLAG_PER && (iP == 0) && (maxDIST > rbox/2.) ) {
    //if( FLAG_PER && (iP == 0) && (minDIST>0.) && (maxDIST > rbox/2.) ) {
    iP=1; //do periodic boundary conditions
    for(int i1=-1; i1<=1; i1++)
      for(int i2=-1; i2<=1; i2++)
	for(int i3=-1; i3<=1; i3++) {
	  vecP[0]=i1*rbox; vecP[1]=i2*rbox; vecP[2]=i3*rbox;
	  dist=0.;
	  for(int j=0; j<dim; j++)
	    dist+=pow(node1->cm[j]-node2->cm[j]+vecP[j], 2);
	  dist=sqrt(dist);
	  maxDIST = dist + node1->rad + node2->rad;
	  minDIST = max(dist-node1->rad-node2->rad, amin);
	  //inject function here
	  if( (minDIST > rHI) || (maxDIST < rLO) ) {
	    //no counts from these nodes
	    continue; //exits the function
	  }

	  else if( (minDIST >= rLO) && (maxDIST < rHI) ) {
	    //include both nodes in the count
	    num+=((node1->u) - (node1->l) + 1.)*((node2->u) - (node2->l) + 1.);
	    continue;
	  }
	  else if( (node1->left == 0) && (node2->left == 0) ) {
	    //leaf nodes
	    for(register int i1=node1->l; i1<=node1->u; i1++)
	      for(register int i2=node2->l; i2<=node2->u; i2++) {
		float rr=0.;
		for(int j=0; j<dim; j++) 
	  	  rr+=pow(Ptree1->Vpart[Ptree1->indexes[i1]].pos[j]-
			  Ptree2->Vpart[Ptree2->indexes[i2]].pos[j]+vecP[j], 2);
		rr=sqrt(rr);
		if( (rr < rHI) && (rr >= rLO) ) num+=1.;
	      }
	    continue;
	  }
	  else if( (node1->left != 0) && (node2->left == 0) ) {
	    //split node 1 only
	    dual_tree_count_PB(node1->left, node2, rLO, rHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2, rLO, rHI, vecP, iP, num);
	    continue;
	  }
	  else if( (node1->left == 0) && (node2->left != 0) ) {
	    //split node 2 only
	    dual_tree_count_PB(node1, node2->left, rLO, rHI, vecP, iP, num);
	    dual_tree_count_PB(node1, node2->right, rLO, rHI, vecP, iP, num);
	    continue;
	  }
	  else if( (node1->left != 0) && (node2->left != 0) )  {
	    //split both nodes
	    dual_tree_count_PB(node1->left, node2->left, rLO, rHI, vecP, iP, num);
	    dual_tree_count_PB(node1->left, node2->right, rLO, rHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2->left, rLO, rHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2->right, rLO, rHI, vecP, iP, num);
	    continue;
	  }
	  
	  //get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minDIST, maxDIST, num) ;
	}
    }
  

  else {
    //get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minDIST, maxDIST, num);
    if( (minDIST > rHI) || (maxDIST < rLO) ) {
      //no counts from these nodes
      return; //exits the function
    }
    
    else if( (minDIST >= rLO) && (maxDIST < rHI) ) {
      //include both nodes in the count
      num+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
      return;
    }
    //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
    else if( (node1->left == 0) && (node2->left == 0) ) {
      //leaf nodes
      for(register int i1=node1->l; i1<=node1->u; i1++)
	for(register int i2=node2->l; i2<=node2->u; i2++) {
	  float rr=0.;
	  for(int j=0; j<dim; j++) 
	    rr+=pow((Ptree1->Vpart[Ptree1->indexes[i1]].pos[j])-
		    (Ptree2->Vpart[Ptree2->indexes[i2]].pos[j])+vecP[j], 2);
	  rr=sqrt(rr);
	  if( (rr < rHI) && (rr >= rLO) ) num+=1.;
	}
      return;
    }
    else if( (node1->left == 0) && (node2->left != 0) ) {
      //split node 2 only
      dual_tree_count_PB(node1, node2->left, rLO, rHI, vecP, iP, num);
      dual_tree_count_PB(node1, node2->right, rLO, rHI, vecP, iP, num);
      return;
    }
    else if( (node1->left != 0) && (node2->left == 0) ) {
      //split node 1 only
      dual_tree_count_PB(node1->left, node2, rLO, rHI, vecP, iP, num); 
      dual_tree_count_PB(node1->right, node2, rLO, rHI, vecP, iP, num);
      return;
    }
    else if( (node1->left != 0) && (node2->left != 0) )  {
      //split both nodes
      dual_tree_count_PB(node1->left, node2->left, rLO, rHI, vecP, iP, num);
      dual_tree_count_PB(node1->left, node2->right, rLO, rHI, vecP, iP, num);
      dual_tree_count_PB(node1->right, node2->left, rLO, rHI, vecP, iP, num);
      dual_tree_count_PB(node1->right, node2->right, rLO, rHI, vecP, iP, num);
      return;
      }
       
  }
} */



void correlations::dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI, 
				      vector<float> vecP, int iP, double &num)
				    
{

  float dist, minDIST, maxDIST;
  //get the max and min distances between the 2 nodes
  dist=0.; //distance bet. the 2 nodes
 
  for(int i=0; i<dim; i++) 
    dist+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
    dist=sqrt(dist);
  maxDIST = dist + node1->rad + node2->rad;
  minDIST = max(dist-node1->rad-node2->rad, amin);
  
  //if( FLAG_PER && (iP == 0) && (minDIST > rHI) && (maxDIST > rbox/2.) ) {
  if( FLAG_PER && (iP == 0) && (maxDIST > rbox/2.) ) {
  //if( FLAG_PER && (iP == 0) && (minDIST>0.) && (maxDIST > rbox/2.) ) {
    iP=1; //do periodic boundary conditions
    for(int i1=-1; i1<=1; i1++)
      for(int i2=-1; i2<=1; i2++)
	for(int i3=-1; i3<=1; i3++) {
	  vecP[0]=i1*rbox; vecP[1]=i2*rbox; vecP[2]=i3*rbox;
	  dist=0.;
	  for(int j=0; j<dim; j++)
	    dist+=pow(node1->cm[j]-node2->cm[j]+vecP[j], 2);
	  dist=sqrt(dist);
	  maxDIST = dist + node1->rad + node2->rad;
	  minDIST = max(dist-node1->rad-node2->rad, amin);
	  get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minDIST, maxDIST, num) ;
	}
    }
  

  else 
    get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minDIST, maxDIST, num);

    }
 

void correlations::get_pairs_dual_tree(tree_node *node1, tree_node *node2, float rLO, float rHI,
				       vector<float> vecP, int iP, float minDIST, float maxDIST, double &num)

{

  if( (minDIST > rHI) || (maxDIST < rLO) ) {
    //no counts from these nodes
    return; //exits the function
  }

  else if( (minDIST >= rLO) && (maxDIST < rHI) ) {
    //include both nodes in the count
    num+=((node1->u) - (node1->l) + 1.)*((node2->u) - (node2->l) + 1.);
    return;
  }
  //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.;
	for(int j=0; j<dim; j++)
	  rr+=pow(Ptree1->Vpart[Ptree1->indexes[i1]].pos[j]-
		  Ptree2->Vpart[Ptree2->indexes[i2]].pos[j]+vecP[j], 2);
      rr=sqrt(rr);
      if( (rr < rHI) && (rr >= rLO) ) num+=1.;
      }
    return;
  }

  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count_PB(node1, node2->left, rLO, rHI, vecP, iP, num);
    dual_tree_count_PB(node1, node2->right, rLO, rHI, vecP, iP, num);
    return;
  }

  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count_PB(node1->left, node2, rLO, rHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2, rLO, rHI, vecP, iP, num);
    return;
  }

  else {
    //split both nodes
    dual_tree_count_PB(node1->left, node2->left, rLO, rHI, vecP, iP, num);
    dual_tree_count_PB(node1->left, node2->right, rLO, rHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2->left, rLO, rHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2->right, rLO, rHI, vecP, iP, num);
    return;
  }
 
}



void correlations::dual_tree_count_all_bins(tree_node *node1, tree_node *node2, vector<double> &num)
				    
{
  float dist, minDIST, maxDIST;
  //get the max and min distances between the 2 nodes
  dist=0.; //distance bet. the 2 nodes
  for(int i=0; i<dim; i++)
    dist+=pow(node1->cm[i]-node2->cm[i], 2);
  dist=sqrt(dist);
  maxDIST = dist + node1->rad + node2->rad;
  minDIST = max(dist-node1->rad-node2->rad, amin);
  
  if( (node1->left != 0) && (node2->left != 0) ) {
    //split both nodes
    dual_tree_count_all_bins(node1->left, node2->left, num);
    dual_tree_count_all_bins(node1->left, node2->right, num);
    dual_tree_count_all_bins(node1->right, node2->left, num);
    dual_tree_count_all_bins(node1->right, node2->right, num);
    return;
  }
  
  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count_all_bins(node1, node2->left, num);
    dual_tree_count_all_bins(node1, node2->right, num);
    return;
  }
  
  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count_all_bins(node1->left, node2, num); 
    dual_tree_count_all_bins(node1->right, node2, num);
    return;
  }
  
  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.;
	for(int j=0; j<dim; j++)
	  rr+=pow(Ptree1->Vpart[Ptree1->indexes[i1]].pos[j]-Ptree2->Vpart[Ptree2->indexes[i2]].pos[j], 2);//+vecP[j], 2);
	rr=sqrt(rr);
	
	if( (rr<rmin) || (rr>=rmax) ) continue;
	//if(rr>=rmax) continue;
	int ibin;
	if(bin_type==LOG)
	  ibin=floor(log(rr/rmin)/delta_bin);
	else if(bin_type==LOG10)
	  ibin=floor(log10(rr/rmin)/delta_bin);
	else
	  ibin=floor((rr-rmin)/delta_bin);
	if( ibin<0 || ibin>=NBIN ) {
	  continue; 
	  cerr<<"exception: "<<ibin<<' '<<rr<<endl;
	  //continue; 
	}
	num[ibin]+=1.; //unnormalized pair counts
      }
    return;
  }
  else {
    for(int ibin=0; ibin<NBIN; ibin++) 
      if((maxDIST < rbinsL[ibin]) ) break;
    for(int ibin=0; ibin<NBIN; ibin++) {
      if( (minDIST > rbinsU[ibin]) ) continue; //|| (maxDIST < rbinsL[ibin]) ) continue;
      else if( (minDIST >= rbinsL[ibin]) && (maxDIST < rbinsU[ibin]) ) {
	//include both nodes in the count
	num[ibin]+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
	//return;
      }
    }
  }
  
}




vector<double> correlations::dual_tree_2p(string ss, bool VERBOSE)

{

  vector<double> ncorr(NBIN, 0.);
  double TIME=0.;
  if(FLAG_PER) {
    cerr<<"compute the dual tree two-point correlation ";
    cerr<<"function using periodic boundary conditions"<<endl;
    ss = ss + "_PB.dat";
  }
  else ss = ss + ".dat";
  //string file_log = "time_log_snap_63_samp_10_MM_NB_40.dat";
  //ofstream out1(file_log.c_str(), ios::app);
  ofstream out(ss.c_str(), ios::out);
  for(int ibin=0; ibin<NBIN; ibin++) {
    float rLO, rHI;
    rLO=rbinsL[ibin];
    rHI=rbinsU[ibin];
    if(VERBOSE)
      cerr<<"computing the dual_tree_count for "<<ibin<<' '<<rLO<<' '<<rHI<<endl;
    double num=0.;
    clock_t t1=clock();
    if(FLAG_PER) dual_tree_count_PB(Ptree1->root, Ptree2->root, rLO, rHI, vecPER, iPER, num);
    else 
      dual_tree_count(Ptree1->root, Ptree2->root, rLO, rHI, num);
    clock_t t2=clock();
    if(Npart1==Npart2) num=num/2.; //remove double counting
    out<<rbinsL[ibin]<<' '<<rbinsM[ibin]<<' '<<rbinsU[ibin]<<" "<<num<<" "<<num/weight<<endl;
    if(VERBOSE)
      cerr<<"time taken for DTC for bin "<<ibin<<": "<<static_cast<double>((t2-t1)/1.e+6)<<" seconds "<<endl;
    TIME+=static_cast<double>((t2-t1)/1.e+6);
    //out1<<rbinsL[ibin]<<' '<<rbinsM[ibin]<<' '<<rbinsU[ibin]<<' '<<static_cast<double>((t2-t1)/1.e+6)<<endl;
    //cerr<<"updating the file " + file_log<<endl;
    ncorr[ibin]+=(num/weight); //the return is unweighted, but without double counting
    if(VERBOSE)
      cerr<<" # pair counts "<<ncorr[ibin]<<endl;
    //SOME MPI STUFF HERE, NEED TO CHECK
  }
  out.close();
  //out1.close();
  if(VERBOSE) {
    cerr<<"done writing file " + ss<<endl;
    cerr<<"total time taken for all bins: "<<TIME<<endl<<endl;  
  }
  
  return ncorr;
  
}



vector<double> correlations::dual_tree_2p_all_bins(string ss)

{

  vector<double> ncorr(NBIN, 0.);
  
  if(FLAG_PER) {
    cerr<<"compute the dual tree two-point correlation ";
    cerr<<"function using periodic boundary conditions"<<endl;
    ss = ss + "_PB.dat";
  }
  else ss = ss + ".dat";
  ofstream out(ss.c_str(), ios::out);

  clock_t t1=clock();
  if(FLAG_PER) dual_tree_count_PB_all_bins(Ptree1->root, Ptree2->root, vecPER, iPER, ncorr);
  else dual_tree_count_all_bins(Ptree1->root, Ptree2->root, ncorr);
  clock_t t2=clock();

  if(Npart1==Npart2)
    for(int i=0; i<NBIN; i++) 
      ncorr[i]/=2.; //remove double counting
  
  cerr<<"time taken for DTCAB : "<<(t2-t1)/CLOCKS_PER_SEC<<" seconds "<<endl;
  for(int i=0; i<NBIN; i++)
    out<<rbinsL[i]<<' '<<rbinsM[i]<<' '<<rbinsU[i]<<' '<<ncorr[i]<<' '<<ncorr[i]/weight<<endl;
  out.close();
  cerr<<"done writing file " + ss <<endl;

  return ncorr;
  

}


vector<double> correlations::brute_force_XX(string ss)

{

  cerr<<"start brute-force calculation for the XX pair counts: "<<endl;
  //assume we need a brute force auto-pair counts of Tree1

  if(FLAG_PER) {
    cerr<<"compute the dual tree two-point correlation ";
    cerr<<"function using periodic boundary conditions"<<endl;
    ss = ss + "_PB.dat";
  }
  else ss = ss + ".dat";
  ofstream out(ss.c_str(), ios::out);
  vector<double> pairs(NBIN, 0.);
  for(register long unsigned i1=0; i1<Npart1; i1++)
    for(register long unsigned i2=i1+1; i2<Npart1; i2++) {
      float r=0.;
      for(int d=0; d<dim; d++) {
	float diff = Ptree1->Vpart[i1].pos[d]-Ptree1->Vpart[i2].pos[d];
	if(FLAG_PER) {
	  if(diff > rbox/2.) diff = diff-rbox;
	  else if(diff < -rbox/2.) diff = diff + rbox;
	}
	r+=pow(diff, 2);
      }
      r=sqrt(r);
      if(r<rmin) continue;
      if(r>=rmax) continue;
      int ibin;
      if(bin_type==LOG)
	ibin=floor(log(r/rmin)/delta_bin);
      else if(bin_type==LOG10)
	ibin=floor(log10(r/rmin)/delta_bin);
      else
	ibin=floor((r-rmin)/delta_bin);
      
      if( ibin<0 || ibin>=NBIN ) continue; //cerr<<"exception: "<<ibin<<' '<<dp<<endl;
      pairs[ibin]+=1.; //unnormalized pair counts
    }
  
  cerr<<"result for brute-force XX pair counts: "<<endl;
  for(int i=0; i<NBIN; i++) {
    //pairs[i]/=2;
    cerr<<rbinsL[i]<<' '<<rbinsU[i]<<' '<<pairs[i]<<endl;
  }
  cerr<<endl;
  
  for(int i=0; i<NBIN; i++)
    out<<rbinsL[i]<<' '<<rbinsM[i]<<' '<<rbinsU[i]<<' '<<pairs[i]<<' '<<pairs[i]/weight<<endl;
  out.close();
  cerr<<"done writing file " + ss <<endl;
  
  return pairs;

}



vector<double> correlations::brute_force_XY()

{

  vector<double> pairs(NBIN, 0.);

  for(int i1=0; i1<Npart1; i1++) 
    for(int i2=0; i2<Npart2; i2++) {
      float r=0.;
      for(int d=0; d<dim; d++) {
	float diff = Ptree1->Vpart[i1].pos[d]-Ptree1->Vpart[i2].pos[d];
	if(FLAG_PER) {
	  if(diff > rbox/2.) diff = diff-rbox;
	  else if(diff < -rbox/2.) diff = diff + rbox;
	}
	r+=pow(diff, 2);
      }
      r=sqrt(r);
      if(r<rmin) continue;
      if(r>=rmax) continue;
      int ibin;
      if(bin_type==LOG)
	ibin=floor(log(r/rmin)/delta_bin);
      else if(bin_type==LOG10)
	ibin=floor(log10(r/rmin)/delta_bin);
      else
	ibin=floor((r-rmin)/delta_bin);
      if( ibin<0 || ibin>=NBIN ) continue; //cerr<<"exception: "<<ibin<<' '<<dp<<endl;
      pairs[ibin]+=1.; //unnormalized pair counts
    }
  cerr<<"result for brute-force XY pair counts: "<<endl;
  for(int i=0; i<NBIN; i++)
    cerr<<rbinsL[i]<<' '<<rbinsU[i]<<' '<<pairs[i]<<endl;
  cerr<<endl;

  return pairs; 
  
}


vector<double> correlations::volume_bins()

{

  const double PI=4.*atan(1.);
  vector<double> rr;
  for(int i=0; i<NBIN; i++) 
    rr.push_back(4.*PI/3. * (pow(rbinsU[i], 3)-pow(rbinsL[i], 3)));
  
  return rr;
  
}

vector<double> correlations::generate_RR(bool write_file, const string file_name)

{

  if(!FLAG_PER) {
    cerr<<"the analytical random catalogue can be generated only for periodic boundary conditions."<<endl;
    exit(1);
  }

  //generate random-random pair counts, using PBC and Nran*(Nran-1)/2 pairs, where Nran is the expected random number per 
  //radial shell, i.e. Vshell * Nran = Vshell * Ntot/Vtot, with Ntot matching the number of galaxies  

  ostringstream osNB;
  osNB<<NBIN; 
  const double Vtot=pow(rbox, 3);
  const double PI=4.*atan(1.);
  vector<double> rr;
  for(int i=0; i<NBIN; i++) {
    rr.push_back(4.*PI/3. * (pow(rbinsU[i], 3)-pow(rbinsL[i], 3)));                                           
    cerr<<rbinsM[i]<<' '<<rr[i]<<endl; //shell volume                                        
  }   
  vector<double> pairsR(NBIN);
  for(int i=0; i<NBIN; i++)
    pairsR[i]= rr[i]/Vtot;

  if(write_file) {                                                                                                         
    ofstream out(file_name.c_str(), ios::out);                                                           
    for(int i=0; i<NBIN; i++)                                              
      out<<rbinsL[i]<<' '<<rbinsM[i]<<' '<<rbinsU[i]<<' '<<pairsR[i]<<endl;
    out.close(); 
    cerr<<"done writing file " + file_name<<endl;
  }
  return pairsR;

}




void correlations::dual_tree_count_PB_all_bins(tree_node *node1, tree_node *node2, 
					       vector<float> vecP, int iP, vector<double> &num)

{

 float dist, minDIST, maxDIST;
 dist=0.; //distance bet. the 2 nodes
  for(int i=0; i<dim; i++) 
    dist+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  dist=sqrt(dist);
  maxDIST = dist + node1->rad + node2->rad;
  minDIST = max(dist-node1->rad-node2->rad, amin);
  
  if( FLAG_PER && (iP == 0) && (minDIST > 0.) && (maxDIST > rbox/2.) ) {
    iP=1; //do periodic boundary conditions
    for(int i1=-1; i1<=1; i1++)
      for(int i2=-1; i2<=1; i2++)
	for(int i3=-1; i3<=1; i3++) {
	  vecP[0]=i1*rbox; vecP[1]=i2*rbox; vecP[2]=i3*rbox;
	  dist=0.;
	  for(int j=0; j<dim; j++)
	    dist+=pow(node1->cm[j]-node2->cm[j]+vecP[j], 2);
	  dist=sqrt(dist);
	  maxDIST = dist + node1->rad + node2->rad;
	  minDIST = max(dist-node1->rad-node2->rad, amin);
	  //inject function here
	  get_pairs_dual_tree_all_bins(node1, node2, vecP, iP, minDIST, maxDIST, num);
	}
  }

  else {
    get_pairs_dual_tree_all_bins(node1, node2, vecP, iP, minDIST, maxDIST, num);
    
    /*for(int ibin=0; ibin<NBIN; ibin++) 
      if((maxDIST < rbinsL[ibin]) ) continue;
      else if( (minDIST > rbinsU[ibin]) ) continue; 
      else if( (minDIST >= rbinsL[ibin]) && (maxDIST < rbinsU[ibin]) ) {
	//include both nodes in the count
	num[ibin]+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
	return;
	}*/
   
  }
}
    


void correlations::get_pairs_dual_tree_all_bins(tree_node *node1, tree_node *node2, vector<float> vecP,
						int iP, float minDIST, float maxDIST, vector<double> &num)
{
 
  //cerr<<"inside GTDTAB "<<endl;
  for(int ibin=0; ibin<NBIN; ibin++) 
    if((maxDIST < rbinsL[ibin]) ) continue; 
    else if( (minDIST > rbinsU[ibin]) ) continue; 
    else if( (minDIST >= rbinsL[ibin]) && (maxDIST < rbinsU[ibin]) ) {
      //include both nodes in the count
      num[ibin]+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
      return;    
    }



  if( (node1->left != 0) && (node2->left != 0) ) {
    //split both nodes
    dual_tree_count_PB_all_bins(node1->left, node2->left, vecP, iP, num);
    dual_tree_count_PB_all_bins(node1->left, node2->right, vecP, iP, num);
    dual_tree_count_PB_all_bins(node1->right, node2->left, vecP, iP, num);
    dual_tree_count_PB_all_bins(node1->right, node2->right, vecP, iP, num);
    return;
  }
  
  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count_PB_all_bins(node1, node2->left, vecP, iP, num);
    dual_tree_count_PB_all_bins(node1, node2->right, vecP, iP, num);
    return;
  }
  
  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count_PB_all_bins(node1->left, node2, vecP, iP, num); 
    dual_tree_count_PB_all_bins(node1->right, node2, vecP, iP, num);
    return;
  }

  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.;
	for(int j=0; j<dim; j++)
	  rr+=pow(Ptree1->Vpart[Ptree1->indexes[i1]].pos[j]-
		  Ptree2->Vpart[Ptree2->indexes[i2]].pos[j] + vecP[j], 2);//+vecP[j], 2);
	rr=sqrt(rr);
	if(rr<rmin) continue;
	if(rr>=rmax) continue;
	int ibin;
	if(bin_type==LOG)
	  ibin=floor(log(rr/rmin)/delta_bin);
	else if(bin_type==LOG10)
	  ibin=floor(log10(rr/rmin)/delta_bin);
	else
	  ibin=floor((rr-rmin)/delta_bin);
	if( ibin<0 || ibin>=NBIN ) {
	  continue; 
	  cerr<<"exception: "<<ibin<<' '<<rr<<endl;
	  //continue; 
	}
	num[ibin]+=1.; //unnormalized pair counts
      }
    return;
  }
}
  

