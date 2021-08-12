#include "tree_processing_projected.h"


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

projected_correlations::projected_correlations(build_tree *Ptree1_, build_tree *Ptree2_, projected_correlations::BinType bt_,
					       const float rmin_, const float rmax_, const int NBIN_, const float zmin_, 
					       const float zmax_, const int ZBIN_, bool flag_per_, const vector<float> VP_):
  Ptree1(Ptree1_), Ptree2(Ptree2_), rmin(rmin_), rmax(rmax_), NBIN(NBIN_), zmin(zmin_), zmax(zmax_), ZBIN(ZBIN_),
  FLAG_PER(flag_per_), vecPER(VP_)

{

  amin=0.;
  bin_type=bt_; 
  set_bins_rperp();
  set_bins_z();
  dim=Ptree1->dim;
  rbox=Ptree1->rbox;
  Npart1=Ptree1->Npart;
  Npart2=Ptree2->Npart;
  if(Npart1!=Npart2) weight = Npart1 * Npart2;
  else weight = Npart1*(Npart1-1.)/2.;

  //if pointers Ptree1 and Ptree2 point to different catalogues, must be careful here
  cerr<<"start correlation function code for tree 1 with "<<Npart1<<" particles."
      <<" and tree 2 with "<<Npart2<<endl;
  //Npart1 is used only for the single_tree calc, so 2 catalogues are no problem
  cerr<<"number of levels in tree 1 (root=level 1)"<<Ptree1->count_levels
      <<" and number of levels in tree 2: "<<Ptree2->count_levels<<endl;
  cerr<<"size count_buckets tree1:  "<<Ptree1->count_buckets.size()
      <<" and tree 2: "<<Ptree2->count_buckets.size()<<endl<<endl;
 
}


void projected_correlations::set_bins_rperp()

{

  if(bin_type==LOG) {

    double dd=log(rmax/rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(exp(log(rmin) + dd*i));
      rbinsU.push_back(exp(log(rmin) + dd*(i+1)));
      rbinsM.push_back(exp(log(rmin) + dd*(i+0.5)));
    }
    delta_bin_rperp=dd;

  }
 
  else if(bin_type==LOG10) {
  
    double dd=log10(rmax/rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(pow(10., (log10(rmin) + dd*i)));
      rbinsU.push_back(pow(10., (log10(rmin) + dd*(i+1))));
      rbinsM.push_back(pow(10., (log10(rmin) + dd*(i+0.5))));
    }
    delta_bin_rperp=dd;
  }
  else if(bin_type==LIN) {

    double dd=(rmax-rmin)/NBIN;
    for(int i=0; i<NBIN; i++) {
      rbinsL.push_back(rmin + dd*i);
      rbinsU.push_back(rmin + dd*(i+1));
      rbinsM.push_back(rmin + dd*(i+0.5));
    }
    delta_bin_rperp=dd;
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


void projected_correlations::set_bins_z()

{
  //do only linear option here
  double dd=(zmax-zmin)/ZBIN;
  for(int i=0; i<ZBIN; i++) {
      zbinsL.push_back(zmin + dd*i);
      zbinsU.push_back(zmin + dd*(i+1));
      zbinsM.push_back(zmin + dd*(i+0.5));
  }
  delta_bin_z=dd;

}


vector<float> projected_correlations::return_bins_rperp() 

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


vector<float> projected_correlations::return_bins_z() 

{

  vector<float> V;
  for(int i=0; i<ZBIN; i++)
    V.push_back(zbinsL[i]);
  for(int i=0; i<ZBIN; i++)
    V.push_back(zbinsU[i]);
  for(int i=0; i<ZBIN; i++)
    V.push_back(zbinsM[i]);

  return V;

}



void projected_correlations::dual_tree_count(tree_node *node1, tree_node *node2, float rLO, float rHI, 
					     float zLO, float zHI, double &num)
				     
{

  //get the max and min perpendicular distances between the 2 nodes
  //the assumption is the 0, 1 are x, y, and 2 is z (i.e. r_perp and r_parallel)
  dist_z=dist_rperp=0.; //distance bet. the 2 nodes

  for(int i=0; i<dim-1; i++)
    dist_rperp+=pow(node1->cm[i]-node2->cm[i], 2);
  dist_rperp = sqrt(dist_rperp);
  maxD_rperp = dist_rperp + node1->rad + node2->rad;
  minD_rperp = max(dist_rperp-node1->rad-node2->rad, amin);

  dist_z = node1->cm[dim-1]-node2->cm[dim-1];
  maxD_z = dist_z + node1->rad + node2->rad;
  minD_z = max(dist_z-node1->rad-node2->rad, amin);


  if( (minD_rperp > rHI) || (maxD_rperp < rLO) || (minD_z > zHI) || (maxD_z < zLO) ) {
    //no counts from these nodes
    return; //exits the function
  }
  
  else if( (minD_rperp >= rLO) && (maxD_rperp < rHI )
	   && (minD_z >= zLO) && (maxD_z < zHI) ) {
    //include both nodes in the count
    num+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
    return;
  }
  //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.; float rz;
	for(int j=0; j<dim-1; j++)
	  rr+=pow((Ptree1->Vpart[Ptree1->indexes[i1]].pos[j])-(Ptree2->Vpart[Ptree2->indexes[i2]].pos[j]), 2);
	rr=sqrt(rr);
	rz=Ptree1->Vpart[Ptree1->indexes[i1]].pos[dim-1]-Ptree2->Vpart[Ptree2->indexes[i2]].pos[dim-1];
	if( (rr < rHI) && (rr >= rLO) && (rz < zHI) && (rz >= zLO) ) 
	  num+=1.;
      }
    return;
  }
  
  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count(node1, node2->left, rLO, rHI, zLO, zHI, num);
    dual_tree_count(node1, node2->right, rLO, rHI, zLO, zHI, num);
    return;
  }
  
  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count(node1->left, node2, rLO, rHI, zLO, zHI, num); 
    dual_tree_count(node1->right, node2, rLO, rHI, zLO, zHI, num);
    return;
  }
  else {
    //split both nodes
    dual_tree_count(node1->left, node2->left, rLO, rHI, zLO, zHI, num);
    dual_tree_count(node1->left, node2->right, rLO, rHI, zLO, zHI, num);
    dual_tree_count(node1->right, node2->left, rLO, rHI, zLO, zHI, num);
    dual_tree_count(node1->right, node2->right, rLO, rHI, zLO, zHI, num);
    return;
  }
}




/*void projected_correlations::dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI, 
						float zLO, float zHI, vector<float> vecP, int iP, double &num)
				    
{

  //get the max and min distances between the 2 nodes
  dist_rperp=dist_z=0.; //distance bet. the 2 nodes

  float dist_tot = 0.;
  for(int i=0; i<dim; i++) 
    dist_tot+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  float maxD_tot = dist_tot + node1->rad + node2->rad; //to use in the first if

  for(int i=0; i<dim-1; i++) 
    dist_rperp+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  dist_rperp=sqrt(dist_rperp);
  maxD_rperp = dist_rperp + node1->rad + node2->rad;
  minD_rperp = max(dist_rperp-node1->rad-node2->rad, amin);

  dist_z = node1->cm[dim-1]-node2->cm[dim-1] + vecP[dim-1];
  maxD_z = dist_z + node1->rad + node2->rad;
  minD_z = max(dist_z-node1->rad-node2->rad, amin);

    //if( FLAG_PER && (iP == 0) && (minD_rperp > rHI) && (maxD_rperp > rbox/2.) ) {
  if( FLAG_PER && (iP == 0) && (maxD_tot > rbox/2.) ) {
    iP=1; //do periodic boundary conditions
    for(int i1=-1; i1<=1; i1++)
      for(int i2=-1; i2<=1; i2++)
	for(int i3=-1; i3<=1; i3++) {
	  vecP[0]=i1*rbox; vecP[1]=i2*rbox; vecP[2]=i3*rbox;
	  dist_rperp=dist_z=0.;
	  for(int j=0; j<dim-1; j++)
	    dist_rperp+=pow(node1->cm[j]-node2->cm[j]+vecP[j], 2);
	  dist_rperp=sqrt(dist_rperp);
	  maxD_rperp = dist_rperp + node1->rad + node2->rad;
	  minD_rperp = max(dist_rperp-node1->rad-node2->rad, amin);
	  dist_z=node1->cm[dim-1]-node2->cm[dim-1]+vecP[dim-1];
	  maxD_z = dist_z + node1->rad + node2->rad;
	  minD_z = max(dist_z-node1->rad-node2->rad, amin);
	  //inject function here
	  if( (minD_rperp > rHI) || (maxD_rperp < rLO) || 
	      (minD_z > zHI) || (maxD_z < zLO) ) {
	    //no counts from these nodes
	    continue; //exits the function
	  }

	  else if( (minD_rperp >= rLO) && (maxD_rperp < rHI) && 
		   (minD_z >= zLO) && (maxD_z < zHI) ) {
	    //include both nodes in the count
	    num+=((node1->u) - (node1->l) + 1.)*((node2->u) - (node2->l) + 1.);
	    continue;
	  }
	  else if( (node1->left == 0) && (node2->left == 0) ) {
	    //leaf nodes
	    for(register int i1=node1->l; i1<=node1->u; i1++)
	      for(register int i2=node2->l; i2<=node2->u; i2++) {
		float rr=0.; float rz;
		for(int j=0; j<dim-1; j++) 
	  	  rr+=pow(Ptree1->Vpart[Ptree1->indexes[i1]].pos[j]-
			  Ptree2->Vpart[Ptree2->indexes[i2]].pos[j]+vecP[j], 2);
		rr=sqrt(rr);
		rz=Ptree1->Vpart[Ptree1->indexes[i1]].pos[dim-1]-
		  Ptree2->Vpart[Ptree2->indexes[i2]].pos[dim-1]+vecP[dim-1];
		if( (rr < rHI) && (rr >= rLO) && (rz < zHI) && (rz >= zLO) ) 
		  num+=1.;
	      }
	    continue;
	  }
	  else if( (node1->left != 0) && (node2->left == 0) ) {
	    //split node 1 only
	    dual_tree_count_PB(node1->left, node2, rLO, rHI, zLO, zHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2, rLO, rHI, zLO, zHI, vecP, iP, num);
	    continue;
	  }
	  else if( (node1->left == 0) && (node2->left != 0) ) {
	    //split node 2 only
	    dual_tree_count_PB(node1, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
	    dual_tree_count_PB(node1, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
	    continue;
	  }
	  else if( (node1->left != 0) && (node2->left != 0) )  {
	    //split both nodes
	    dual_tree_count_PB(node1->left, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
	    dual_tree_count_PB(node1->left, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
	    dual_tree_count_PB(node1->right, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
	    continue;
	  }
	  
	  //get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minD_rperp, maxD_rperp, num) ;
	}
  }
  

  else {
    //get_pairs_dual_tree(node1, node2, rLO, rHI, vecP, iP, minD_rperp, maxD_rperp, num) ;
    if( (minD_rperp > rHI) || (maxD_rperp < rLO) || 
	(minD_z > zHI) || (maxD_z < zLO) ) {
      //no counts from these nodes
      return; //exits the function
    }
    
    else if( (minD_rperp >= rLO) && (maxD_rperp < rHI) && 
	     (minD_z >= zLO) && (maxD_z < zHI) ) {
      //include both nodes in the count
      num+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
      return;
    }
    //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
    else if( (node1->left == 0) && (node2->left == 0) ) {
      //leaf nodes
      for(register int i1=node1->l; i1<=node1->u; i1++)
	for(register int i2=node2->l; i2<=node2->u; i2++) {
	  float rr=0.; float rz;
	  for(int j=0; j<dim-1; j++) 
	    rr+=pow((Ptree1->Vpart[Ptree1->indexes[i1]].pos[j])-
		    (Ptree2->Vpart[Ptree2->indexes[i2]].pos[j])+vecP[j], 2);
	  rr=sqrt(rr);
	  rz=(Ptree1->Vpart[Ptree1->indexes[i1]].pos[dim-1])-
	    (Ptree2->Vpart[Ptree2->indexes[i2]].pos[dim-1])+vecP[dim-1];
	  if( (rr < rHI) && (rr >= rLO) && (rz < zHI) && (rz >= zLO) ) 
	    num+=1.;
	}
      return;
    }
    else if( (node1->left == 0) && (node2->left != 0) ) {
      //split node 2 only
      dual_tree_count_PB(node1, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
      dual_tree_count_PB(node1, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
      return;
    }
    else if( (node1->left != 0) && (node2->left == 0) ) {
      //split node 1 only
      dual_tree_count_PB(node1->left, node2, rLO, rHI, zLO, zHI, vecP, iP, num); 
      dual_tree_count_PB(node1->right, node2, rLO, rHI, zLO, zHI, vecP, iP, num);
      return;
    }
    else if( (node1->left != 0) && (node2->left != 0) )  {
      //split both nodes
      dual_tree_count_PB(node1->left, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
      dual_tree_count_PB(node1->left, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
      dual_tree_count_PB(node1->right, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
      dual_tree_count_PB(node1->right, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
      return;
    }
       
  }
}
*/ 

void projected_correlations::dual_tree_count_PB(tree_node *node1, tree_node *node2, float rLO, float rHI,
                                                float zLO, float zHI, vector<float> vecP, int iP, double &num)

{
  
  //get the max and min distances between the 2 nodes
  dist_rperp=dist_z=0.; //distance bet. the 2 nodes
  
  float dist_tot = 0.;
  for(int i=0; i<dim; i++)
    dist_tot+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  dist_tot=sqrt(dist_tot); //WHY IS THIS COMMENTED OUT?
  float maxD_tot = dist_tot + node1->rad + node2->rad; //to use in the first if
  
  for(int i=0; i<dim-1; i++)
    dist_rperp+=pow(node1->cm[i]-node2->cm[i] + vecP[i], 2);
  dist_rperp=sqrt(dist_rperp);
  maxD_rperp = dist_rperp + node1->rad + node2->rad;
  minD_rperp = max(dist_rperp-node1->rad-node2->rad, amin);
  
  dist_z = node1->cm[dim-1]-node2->cm[dim-1] + vecP[dim-1];
  maxD_z = dist_z + node1->rad + node2->rad;
  minD_z = max(dist_z-node1->rad-node2->rad, amin);

  //if( FLAG_PER && (iP == 0) && (minD_rperp > rHI) && (maxD_rperp > rbox/2.) ) {
  if( FLAG_PER && (iP == 0) && (maxD_tot > rbox/2.) ) {
    iP=1; //do periodic boundary conditions
    for(int i1=-1; i1<=1; i1++)
      for(int i2=-1; i2<=1; i2++)
        for(int i3=-1; i3<=1; i3++) {
          vecP[0]=i1*rbox; vecP[1]=i2*rbox; vecP[2]=i3*rbox;
          dist_rperp=dist_z=0.;
          for(int j=0; j<dim-1; j++)
            dist_rperp+=pow(node1->cm[j]-node2->cm[j]+vecP[j], 2);
          dist_rperp=sqrt(dist_rperp);
          maxD_rperp = dist_rperp + node1->rad + node2->rad;
          minD_rperp = max(dist_rperp-node1->rad-node2->rad, amin);
          dist_z=node1->cm[dim-1]-node2->cm[dim-1]+vecP[dim-1];
          maxD_z = dist_z + node1->rad + node2->rad;
          minD_z = max(dist_z-node1->rad-node2->rad, amin);
          //inject function here
	  get_pairs_dual_tree(node1, node2, rLO, rHI, zLO, zHI, vecP, iP, 
			      minD_rperp, maxD_rperp, minD_z, maxD_z, num);

	}
  }
  
  else 
    get_pairs_dual_tree(node1, node2, rLO, rHI, zLO, zHI, vecP, iP,
			minD_rperp, maxD_rperp, minD_z, maxD_z, num);

}


void projected_correlations::get_pairs_dual_tree(tree_node *node1, tree_node *node2, float rLO, float rHI,
						 float zLO, float zHI, vector<float> vecP, int iP, float minD_rperp, 
						 float maxD_rperp, float minD_z, float maxD_z, double &num)

{

  if( (minD_rperp > rHI) || (maxD_rperp < rLO) ||
      (minD_z > zHI) || (maxD_z < zLO) ) {
    //no counts from these nodes
    return; //exits the function
  }

  else if( (minD_rperp >= rLO) && (maxD_rperp < rHI) &&
	   (minD_z >= zLO) && (maxD_z < zHI) ) {
    //include both nodes in the count
    num+=((node1->u) - (node1->l) + 1.)*( (node2->u) - (node2->l) + 1.); //+=
    return;
  }
  //the assumption is that the left node is emblematic for the right one too, i.e. if left=0, then right=0
  else if( (node1->left == 0) && (node2->left == 0) ) {
    //leaf nodes
    for(register int i1=node1->l; i1<=node1->u; i1++)
      for(register int i2=node2->l; i2<=node2->u; i2++) {
	float rr=0.; float rz;
	for(int j=0; j<dim-1; j++)
	  rr+=pow((Ptree1->Vpart[Ptree1->indexes[i1]].pos[j])-
		  (Ptree2->Vpart[Ptree2->indexes[i2]].pos[j])+vecP[j], 2);
	rr=sqrt(rr);
	rz=(Ptree1->Vpart[Ptree1->indexes[i1]].pos[dim-1])-
	  (Ptree2->Vpart[Ptree2->indexes[i2]].pos[dim-1])+vecP[dim-1];
	if( (rr < rHI) && (rr >= rLO) && (rz < zHI) && (rz >= zLO) )
	  num+=1.;
      }
    return;
  }

  else if( (node1->left == 0) && (node2->left != 0) ) {
    //split node 2 only
    dual_tree_count_PB(node1, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
    dual_tree_count_PB(node1, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
    return;
  }
  else if( (node1->left != 0) && (node2->left == 0) ) {
    //split node 1 only
    dual_tree_count_PB(node1->left, node2, rLO, rHI, zLO, zHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2, rLO, rHI, zLO, zHI, vecP, iP, num);
    return;
  }
  else if( (node1->left != 0) && (node2->left != 0) )  {
    //split both nodes
    dual_tree_count_PB(node1->left, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
    dual_tree_count_PB(node1->left, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2->left, rLO, rHI, zLO, zHI, vecP, iP, num);
    dual_tree_count_PB(node1->right, node2->right, rLO, rHI, zLO, zHI, vecP, iP, num);
    return;
  }

}


vector<double> projected_correlations::dual_tree_2p(string ss)

{

  vector<double> ncorr(NBIN*ZBIN, 0.);
  double TIME=0.;
  if(FLAG_PER) {
    cerr<<"compute the dual tree two-point correlation ";
    cerr<<"function using periodic boundary conditions"<<endl;
    ss = ss + "_PB.dat";
  }
  else ss = ss + ".dat";

  //ofstream out(ss.c_str(), ios::out);
  for(int ibin=0; ibin<NBIN; ibin++)
    for(int ibiz=0; ibiz<ZBIN; ibiz++) {
      float rLO, rHI, zLO, zHI;
      rLO=rbinsL[ibin]; zLO=zbinsL[ibiz];
      rHI=rbinsU[ibin]; zHI=zbinsU[ibiz];

      double num=0.;
      clock_t t1=clock();
      if(FLAG_PER) dual_tree_count_PB(Ptree1->root, Ptree2->root, rLO, rHI, zLO, zHI, vecPER, iPER, num);
      else dual_tree_count(Ptree1->root, Ptree2->root, rLO, rHI, zLO, zHI, num);
      clock_t t2=clock();
      if(Npart1==Npart2) num=num/2.; //remove double counting
      //out<<rbinsL[ibin]<<' '<<rbinsM[ibin]<<' '<<rbinsU[ibin]<<' '<<zbinsL[ibiz]<<' '
      //<<zbinsM[ibiz]<<' '<<zbinsU[ibiz]<<' '<<num<<' '<<num/weight<<endl;
      //cerr<<"time taken for DTC for binR "<<ibin<<" for binZ: "<<ibiz<<' '<<static_cast<double>((t2-t1)/1.e+6)<<" seconds "<<endl;
      //TIME+=static_cast<double>((t2-t1)/1.e+6);
      ncorr[ibin*ZBIN+ibiz]+=(num/weight); //the return is unweighted, but without double counting

    }
  const double PI=4.*atan(1.);
  double Norm=0.;
  vector<double> rr(NBIN*ZBIN, 0.);
  for(int i=0; i<NBIN; i++)
    for(int j=0; j<ZBIN; j++) 
      rr[i*ZBIN+j]=PI * (pow(rbinsU[i], 2)-pow(rbinsL[i], 2)) * (zbinsU[j]-zbinsL[j]);
  for(int ir=0; ir<NBIN; ir++)
    for(int iz=0; iz<ZBIN; iz++)
      Norm+=(ncorr[ir*ZBIN+iz]*rr[ir*ZBIN+iz]);
  IntegralConstraint=Norm;
  
  //out.close();

  return ncorr;
  
}

vector<double> projected_correlations::volume_bins()

{
  const double PI=4.*atan(1.);
  vector<double> rr(NBIN*ZBIN, 0.);
  for(int i=0; i<NBIN; i++)
    for(int j=0; j<ZBIN; j++) 
      rr[i*ZBIN+j]=PI * (pow(rbinsU[i], 2)-pow(rbinsL[i], 2)) * (zbinsU[j]-zbinsL[j]);
  
  return rr;

}


vector<double> projected_correlations::generate_RR(bool write_file)

{

  if(!FLAG_PER) {
    cerr<<"the analytical random catalogue can be generated only for periodic boundary conditions."<<endl;
    exit(1);
  }

  //generate random-random pair counts, using PBC and Nran*(Nran-1)/2 pairs, where Nran is the expected random number per radial shell, i.e. Vshell * Nran = Vshell * Ntot/Vtot, with Ntot matching the number of galaxies
  
  ostringstream oZMAX, osNB, osNZ;
  oZMAX<<zbinsU[ZBIN-1]; osNB<<NBIN; osNZ<<ZBIN;
  const double Vtot=pow(rbox, 3);
  const double PI=4.*atan(1.);
  vector<double> rr(NBIN*ZBIN, 0.), pairsR(NBIN*ZBIN);
  for(int i=0; i<NBIN; i++) 
    for(int j=0; j<ZBIN; j++) {
      rr[i*ZBIN+j]=PI * (pow(rbinsU[i], 2)-pow(rbinsL[i], 2)) * (zbinsU[j]-zbinsL[j]);
      pairsR[i*ZBIN+j]= rr[i*ZBIN+j]/Vtot;
    }

  /*if(write_file) {
    const string path_out2 = "/u/lmarian/LauraTreeCode/Zurich/RESULTS/PROJECTED_CORR_FUNC/RANDOM/";
    const string file_rr = path_out2 + "analytical_projected_random_pair_counts_rob_binsR_" + osNB.str()  
    + "_binsZ_" + osNZ.str() + "_zmax_" + oZMAX.str() + ".dat";
    ofstream out(file_rr.c_str(), ios::out);
    for(int i=0; i<NBIN; i++) 
      for(int j=0; j<ZBIN; j++)
	out<<rbinsL[i]<<' '<<rbinsM[i]<<' '<<rbinsU[i]<<' '<<zbinsL[j]<<' '<<zbinsM[j]<<' '<<zbinsU[j]<<' '<<rr[i*ZBIN+j]/Vtot<<endl;
    out.close();
    cerr<<"done writing file " + file_rr<<endl;
    }*/

  return pairsR;

}
