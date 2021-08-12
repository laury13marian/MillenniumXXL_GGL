#include "JK_statistics.h"
//run it with TASK=1000 if you already work with the average of the tasks

using namespace std;

//this is for a subcube of the XXL 
JK_statistics::JK_statistics(const string path_gen_, vector<double> mags_, const int sub_, int TASK):
  path_gen(path_gen_), mags(mags_), sub(sub_)
  //TASK=1000 if no separate tasks are analyzed, otherwise TASK starts at 0
{

  Vgg.resize(NJK*NM*NR); SHgm.resize(NJK*NM*NR);
  binsR.resize(NR);

  osSNAP<<0; osSNAP<<snap; osNR<<NR; osNM<<NM;           
  osT<<ntasks; osSC<<sub;

  cerr<<"READING THE JACK KNIFE DATA for the jack-knife analysis "<<endl;
  for(int is=0; is<NJK; is++)
    read_data(is, TASK);

}


void JK_statistics::read_data(int iJK, int TASK) //iJK from 0 to NJK-1
//start at JKsample=0
{
  
  ostringstream osJK, osTASK;
  osJK<<iJK;  osTASK<<TASK;

  string file_in_GG, file_in_SHgm;
  float f;
  double d, b1, b2;
  
  if(TASK!=1000) { //that means that we are really using the tasks
    
    file_in_GG = path_gen + "CORR_GG/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_NT_" + osT.str() + "_one_task_" + osTASK.str() + ".dat";
    
    file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/ONE_TASK/V2/delta_sigma_gm_no_rhob_sub_" 
      + osSC.str() + "_one_task_" + osTASK.str() + ".dat";
  }
  
  else {
    
    file_in_GG = path_gen + "CORR_GG/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str() + "_JKsample_" 
      + osJK.str() + "_binsM_" + osNM.str() + "_binsR_" + osNR.str()  + "_NT_" + osT.str() + ".dat";
    
    file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/delta_sigma_gm_no_rhob_sub_" + osSC.str() 
      + "_JKsample_" + osJK.str() + ".dat";
  }
  
  ifstream inG(file_in_GG.c_str(), ios::in);
  check_stream(inG, file_in_GG);
  for(int ir=0; ir<NR; ir++) {
    inG>>f;
    binsR[ir]=f;
    for(int im=0; im<NM; im++) {
      inG>>b1>>b2>>d;
      Vgg[(iJK*NM+im)*NR+ir]=d;
    }
  }
  inG.close();
  
  //read the delta-sigma files
  
  ifstream inSHgm(file_in_SHgm.c_str(), ios::in);
  check_stream(inSHgm, file_in_SHgm);
  for(int ir=0; ir<NR; ir++) {
    inSHgm>>f;
    for(int im=0; im<NM; im++) {
      inSHgm>>b1>>b2>>d;
      SHgm[(iJK*NM+im)*NR+ir]=d;
    }
  }
  inSHgm.close();
  
}


void JK_statistics::signal_average()
  
{

  if(Agg.size() != 0 && Ashgm.size() != 0) {
    cerr<<"average jack-knifed functions already computed. "<<endl;
    return;
  }
  
  vector<double> Tave;
  //galaxy-galaxy average and standard deviation

  Tave.resize(NR*NM);
  for(int i=0; i<NR*NM; i++) 
    Tave[i]=0.; 

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<NJK; is++)
	Tave[im*NR+ir]+=Vgg[(is*NM+im)*NR+ir];
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
        Tave[im*NR+ir]+=SHgm[(is*NM+im)*NR+ir];
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
      cerr<<Agg[im*NR+ir]<<' '<<Ashgm[im*NR+ir]<<' ';
    cerr<<endl;
  }
}


vector<SqDMatrix> JK_statistics::signal_covariance(string type)

{

  signal_average();
  vector<SqDMatrix> VecCov;
  vector<double> V1, V2, A1, A2;
 
  if(type=="gggg") {
    V1=V2=Vgg; A1=A2=Agg; 
  }
  else if(type=="sgmsgm") {
    V1=V2=SHgm; A1=A2=Ashgm;
  }
  else if(type=="ggsgm") {
    V1=Vgg; V2=SHgm; A1=Agg; A2=Ashgm;
  }
  else {
    cerr<<"unknown covariance type. Allowed types are: "<<endl;
    cerr<<"gggg, sgmsgm, ggsgm."<<endl;
    exit(1);
  } 

  for(int im=0; im<NM; im++) {
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

void JK_statistics::write_covariance_matrix(vector<string> file, string type)

{

  int AL=2; //eliminate last 2 radial bins, they are probably not safe for Lproj=100 Mpc/h
  int NR2=NR-AL;
  vector<SqDMatrix> VecC;
  VecC=signal_covariance(type);

  if(VecC.size() != NM || NM != file.size()) {
    cerr<<"in write_covariance: wrong size for the covariance matrix or output file vectors"<<endl;
    exit(1);
  } 

  for(int im=0; im<NM; im++) {
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


void JK_statistics::write_correlation_matrix(vector<string> file, string type)

{

  int AL=2; //eliminate last 2 radial bins, they are probably not safe for Lproj=100 Mpc/h
  int NR2=NR-AL;
  vector<SqDMatrix> VecC;

  VecC=signal_covariance(type);
  if(VecC.size() != NM || NM != file.size()) {
    cerr<<"in write_covariance: wrong size for the covariance matrix or output file vectors"<<endl;
    exit(1);
  } 

  if(type=="gggg" || type=="sgmsgm") {
    for(int im=0; im<NM; im++) {
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
    for(int im=0; im<NM; im++) {
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

