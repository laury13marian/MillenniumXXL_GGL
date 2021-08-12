#include "jack_knife.h"

using namespace std;

//this is for the XXL volume, using the results from all subcubes
jack_knife::jack_knife(const string path_gen_, vector<double> mags_, int block_size_, int TASK):
  path_gen(path_gen_), mags(mags_), NM(mags_.size()-1), block_size(block_size_)
  //TASK=1000 if no separate tasks are analyzed, otherwise TASK starts at 0
{

  if(Nsub % block_size !=0) {
    cerr<<"block_size must exactly divide the number of subcubes "<<Nsub<<endl;
    exit(1);
  }
  else 
    NJK=Nsub/block_size;
  cerr<<"the number of jack-knife samplings is : "<<NJK<<endl;

  Vgg.resize(Nsub*NM*NR); SHgm.resize(Nsub*NM*NR);
  VJKgg.resize(NJK*NM*NR); VJKshgm.resize(NJK*NM*NR);
  binsR.resize(NR);

  osSNAP<<0; osSNAP<<snap; osNZ<<NZ; osNR<<NR; osNM<<NM; //CAREFUL HERE!!!!               
  osNSAMP<<nsampM; osT<<ntasks; osNRAN<<NRAN; osZMAX<<zmax;

  cerr<<"READING THE SUBCUBE DATA for the jack-knife analysis "<<endl;
  for(int is=0; is<Nsub; is++)
    read_data(is, TASK);

}

//estimate the statistics of interest (wp or shear) with one block of data removed
//the blocks are made of one or more subcubes, and they are removed sequentially
void jack_knife::jack_knife_data(int ijack) //also for the XXL volume
//ijack starts at 0, i.e. first block to be removed, then second, etc.
{
  if(ijack>=NJK) {
    cerr<<"ijack must be smaller than the number of jack-knife samplings "<<NJK<<endl;
    exit(1);
  }
  vector<int> Vcount(Nsub);
  for(int j=0; j<Nsub; j++)
    Vcount[j]=j;
  vector<int>::iterator p;
  p=Vcount.begin();
  p+=ijack*block_size;
  Vcount.erase(p, p+block_size);
  cerr<<"new size of Vcount: "<<Vcount.size()<<endl<<endl;

  //now compute averages from the remaining subcubes
  //Vgg and SHgm only, since they generate the covariances of interest

  vector<double> VTgg(NM*NR, 0.), VTshgm(NM*NR, 0.);
  for(int im=0; im<NM; im++) 
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Vcount.size(); is++) { 
	VTgg[im*NR+ir]+=Vgg[(Vcount[is]*NM+im)*NR+ir];
	VTshgm[im*NR+ir]+=SHgm[(Vcount[is]*NM+im)*NR+ir];
      }
      VTgg[im*NR+ir]/=Vcount.size();
      VTshgm[im*NR+ir]/=Vcount.size();
      VJKgg[(ijack*NM+im)*NR+ir]=VTgg[im*NR+ir];
      VJKshgm[(ijack*NM+im)*NR+ir]=VTshgm[im*NR+ir];
    }

}


void jack_knife::jack_knife_data() //also for the XXL volume

{

  for(int ijack=0; ijack<NJK; ijack++)
    jack_knife_data(ijack);
  //check to see if the new data vectors have all elements non-trivial
  for(int j=0; j<NJK*NM*NR; j++) {
    if(VJKgg[j]==0.) 
      cerr<<"JK vector for gg is 0. at position: "<<j<<endl;
    if(VJKshgm[j]==0.)
      cerr<<"JK vector for shear gm is 0. at position: "<<j<<endl;
  }

}


void jack_knife::signal_average()

{

  if(Agg.size() != 0 && Ashgm.size() != 0) {
    cerr<<"average jack-knifed functions already computed. "<<endl;
    return;
  }
  
  jack_knife_data();

  vector<double> Tave;
  //galaxy-galaxy average and standard deviation

  Tave.resize(NR*NM);
  for(int i=0; i<NR*NM; i++) 
    Tave[i]=0.; 

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<NJK; is++)
	Tave[im*NR+ir]+=VJKgg[(is*NM+im)*NR+ir];
      Tave[im*NR+ir]/=NJK;
      Agg.push_back(Tave[im*NR+ir]);
    }

  Tave.clear(); 

 //do the SHgm average

  Tave.resize(NR*NM);
  for(int i=0; i<NR*NM; i++) 
    Tave[i]=0.; 

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<NJK; is++)
        Tave[im*NR+ir]+=VJKshgm[(is*NM+im)*NR+ir];
      Tave[im*NR+ir]/=NJK;
      Ashgm.push_back(Tave[im*NR+ir]);
    }

  Tave.clear();
  
  cerr<<
    "done computing the jack-knife average for the projected gg correlation function and the gm shear."
      <<endl;
 
  cerr<<"jack-knife averages: "<<endl;
  for(int ir=0; ir<NR; ir++) {
    cerr<<binsR[ir]<<' ';
    for(int im=0; im<NM; im++)
      cerr<<Agg[im*NR+ir]<<' ';
    cerr<<endl;
  }
  cerr<<endl;
  for(int ir=0; ir<NR; ir++) {
    cerr<<binsR[ir]<<' ';
    for(int im=0; im<NM; im++)
      cerr<<Ashgm[im*NR+ir]<<' ';
    cerr<<endl;
  }

}


vector<SqDMatrix> jack_knife::signal_covariance(string type)

{

  signal_average();
  vector<SqDMatrix> VecCov;
  vector<double> V1, V2, A1, A2;
 
  if(type=="gggg") {
    V1=V2=VJKgg; A1=A2=Agg; 
  }
  else if(type=="sgmsgm") {
    V1=V2=VJKshgm; A1=A2=Ashgm;
  }
  else if(type=="ggsgm") {
    V1=VJKgg; V2=VJKshgm; A1=Agg; A2=Ashgm;
  }
  else {
    cerr<<"unknown covariance type. Allowed types are: "<<endl;
    cerr<<"gggg, sgmsgm, ggsgm."<<endl;
    exit(1);
  } 

  for(int im=0; im<NM2; im++) {
    SqDMatrix Cov(NR, 0.);
    for(int ir1=0; ir1<NR; ir1++)
      for(int ir2=0; ir2<NR; ir2++) {
	for(int is=0; is<NJK; is++) 
	  Cov(ir1, ir2)+=(V1[(is*NM+im)*NR+ir1]-A1[im*NR+ir1]) 
	    * (V2[(is*NM+im)*NR+ir2]-A2[im*NR+ir2]);
	Cov(ir1, ir2)*=((NJK-1.)/NJK); //this probably gives the covariance on the mean
      }
    VecCov.push_back(Cov);
  }
  
  return VecCov;
  
}

void jack_knife::write_covariance_matrix(vector<string> file, string type)

{

  int AL=2; //eliminate last 2 radial bins, they are probably not safe for Lproj=100 Mpc/h
  int NR2=NR-AL;
  vector<SqDMatrix> VecC;
  VecC=signal_covariance(type);

  if(VecC.size() != NM2 || NM2 != file.size()) {
    cerr<<"in write_covariance: wrong size for the covariance matrix or output file vectors"<<endl;
    exit(1);
  } 

  for(int im=0; im<NM2; im++) {
    ofstream out(file[im].c_str(), ios::out);
    out<<NR2<<endl;
    for(int ir=0; ir<NR2; ir++)
      out<<binsR[ir]<<endl;
    for(int ir1=0; ir1<NR2; ir1++) {
      for(int ir2=0; ir2<NR2; ir2++)
	out<<VecC[im](ir1, ir2)<<' ';
      out<<endl;
    }
    out.close();
    cerr<<"done writing file " + file[im]<<endl;
  }
  /*cerr<<"covariance matrix for gg, im=0:"<<endl;
  for(int ir1=0; ir1<NR2; ir1++) {
    for(int ir2=0; ir2<NR2; ir2++)
      cerr<<VecC[0](ir1, ir2)<<' ';
    cerr<<endl;
  }
  cerr<<endl;*/

}


void jack_knife::write_correlation_matrix(vector<string> file, string type)

{

  int AL=2; //eliminate last 2 radial bins, they are probably not safe for Lproj=100 Mpc/h
  int NR2=NR-AL;
  vector<SqDMatrix> VecC;

  VecC=signal_covariance(type);
  if(VecC.size() != NM2 || NM2 != file.size()) {
    cerr<<"in write_covariance: wrong size for the covariance matrix or output file vectors"<<endl;
    exit(1);
  } 

  if(type=="gggg" || type=="sgmsgm") {
    for(int im=0; im<NM2; im++) {
      ofstream out(file[im].c_str(), ios::out);
      out<<NR2<<endl;
      for(int ir=0; ir<NR2; ir++)
	out<<binsR[ir]<<endl;
      for(int ir1=0; ir1<NR2; ir1++) {
	for(int ir2=0; ir2<NR2; ir2++)
	  out<<VecC[im](ir1, ir2)/sqrt(VecC[im](ir1, ir1)*VecC[im](ir2, ir2))<<' ';
	out<<endl;
      }
      out.close();
      cerr<<"done writing file " + file[im]<<endl;
    }
  }

  else if(type=="ggsgm") {
    vector<SqDMatrix> MatGG, MatSGM;
    MatGG=signal_covariance("gggg");
    MatSGM=signal_covariance("sgmsgm");
    for(int im=0; im<NM2; im++) {
      ofstream out(file[im].c_str(), ios::out);
      out<<NR2<<endl;
      for(int ir=0; ir<NR2; ir++)
	out<<binsR[ir]<<endl;
      for(int ir1=0; ir1<NR2; ir1++) {
	for(int ir2=0; ir2<NR2; ir2++)
	  out<<VecC[im](ir1, ir2)/sqrt(MatGG[im](ir1, ir1))/sqrt(MatSGM[im](ir2, ir2))<<' ';
	out<<endl;
      }
      out.close();
      cerr<<"done writing file " + file[im]<<endl;
    }
  }

}


void jack_knife::read_data(int sub, int TASK)
//start at sub=0
{
  
  ostringstream osSC, osTASK;
  osSC<<sub;  osTASK<<TASK;

  string file_in_GG, file_in_SHgm;
  float f;
  double d1, d2, d3, b1, b2;
  
 if(TASK!=1000) { //that means that we are really using the tasks

   file_in_GG = path_gen + "CORR_GG/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
     + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
     + osNRAN.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + "_one_task_" + osTASK.str() + ".dat";
   
   file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/ONE_TASK/V2/delta_sigma_gm_no_rhob_sub_" 
     + osSC.str() + "_one_task_" + osTASK.str() + ".dat";
 }

 else {
   
   file_in_GG = path_gen + "CORR_GG/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
     + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
     + osNRAN.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + ".dat";
   
   file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/V2/delta_sigma_gm_no_rhob_sub_" + osSC.str() + ".dat";
 }
 
 ifstream inG(file_in_GG.c_str(), ios::in);
 check_stream(inG, file_in_GG);
 for(int ir=0; ir<NR; ir++) {
    inG>>f;
    if(sub==0) 
      binsR[ir]=f;
    for(int im=0; im<NM; im++) {
      inG>>b1>>b2>>d1>>d2>>d3;
      Vgg[(sub*NM+im)*NR+ir]=d3;
    }
  }
 inG.close();
 
  //read the delta-sigma files

  ifstream inSHgm(file_in_SHgm.c_str(), ios::in);
  check_stream(inSHgm, file_in_SHgm);
  for(int ir=0; ir<NR; ir++) {
    inSHgm>>f;
    for(int im=0; im<NM; im++) {
      inSHgm>>b1>>b2>>d3;
      SHgm[(sub*NM+im)*NR+ir]=d3;
    }
  }
  inSHgm.close();

}

