#include "statistics.h"

using namespace std;


statistics::statistics(const string path_gen_, const int ic_, vector<double> mags_, int TASK):
  path_gen(path_gen_), FLAG_IC(ic_), mags(mags_), NM(mags_.size()-1)
  //TASK=1000 if no separate tasks are analyzed, otherwise TASK starts at 0
{

  Vmm.resize(Nsub*NR); Vgg.resize(Nsub*NM*NR); Vgm.resize(Nsub*NM*NR);
  SHmm.resize(Nsub*NR); SHgm.resize(Nsub*NM*NR);
  binsR.resize(NR);

  osSNAP<<0; osSNAP<<snap; osNZ<<NZ; osNR<<NR; osNM<<NM; //CAREFUL HERE!!!!               
  osNSAMP<<nsampM; osT<<ntasks; osNRAN<<NRAN; osZMAX<<zmax;

  cerr<<"READING THE SUBCUBE DATA: "<<endl;
  for(int is=0; is<Nsub; is++)
    read_data(is, TASK);

}



void statistics::read_data(int sub, int TASK) //TASK=1000 if no separate tasks are analyzed
//start at sub=0
{
  
  ostringstream osSC, osTASK;
  osSC<<sub;
  osTASK<<TASK;
  string file_in_MM, file_in_GG, file_in_GM, file_in_SHmm, file_in_SHgm;
  
  if(TASK!=1000) { //that means that we are really using the tasks

    file_in_MM = path_gen + "CORR_MM/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_"
      + osNSAMP.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + "_one_task_" + osTASK.str() + ".dat";
    
    file_in_GG = path_gen + "CORR_GG/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
      + osNRAN.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + "_one_task_" + osTASK.str() + ".dat";
    
    file_in_GM = path_gen + "CORR_GM/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
      + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + "_one_task_" 
      + osTASK.str() + ".dat";
    
    file_in_SHmm = path_gen + "CORR_MM/DELTA_SIGMA/ONE_TASK/delta_sigma_mm_without_rhob_sub_"
      + osSC.str() + "_one_task_" + osTASK.str() + ".dat"; 
    file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/ONE_TASK/V2/delta_sigma_gm_no_rhob_sub_" 
      + osSC.str() + "_one_task_" + osTASK.str() + ".dat";
    
  }
  
  else {
    file_in_MM = path_gen + "CORR_MM/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_"
      + osNSAMP.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + ".dat";
    
    file_in_GG = path_gen + "CORR_GG/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
      + osNRAN.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + ".dat";
  
    file_in_GM = path_gen + "CORR_GM/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + osSC.str()
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
      + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + osZMAX.str() + "_NT_" + osT.str() + ".dat";

    file_in_SHmm = path_gen + "CORR_MM/DELTA_SIGMA/delta_sigma_mm_without_rhob_sub_" + osSC.str() + ".dat";
    file_in_SHgm = path_gen + "CORR_GM/DELTA_SIGMA/V2/delta_sigma_gm_no_rhob_sub_" + osSC.str() + ".dat";

  }
    
  ifstream inM(file_in_MM.c_str(), ios::in);
  check_stream(inM, file_in_MM);

  float f;
  double d1, d2, d3, b1, b2;
  for(int ir=0; ir<NR; ir++) { 
    inM>>binsR[ir]>>d1>>d2>>d3;
    if(FLAG_IC==0) 
      Vmm[sub*NR+ir]=d1;
    else if(FLAG_IC==1)
      Vmm[sub*NR+ir]=d2;
    else  
      Vmm[sub*NR+ir]=d3;
  }
  inM.close();

  ifstream inG(file_in_GG.c_str(), ios::in);
  check_stream(inG, file_in_GG);
  for(int ir=0; ir<NR; ir++) {
    inG>>f;
    for(int im=0; im<NM; im++) {
      inG>>b1>>b2>>d1>>d2>>d3;
      if(FLAG_IC==0)
	Vgg[(sub*NM+im)*NR+ir]=d1;
      else if(FLAG_IC==1)
	Vgg[(sub*NM+im)*NR+ir]=d2;
      else
	Vgg[(sub*NM+im)*NR+ir]=d3;
    }
  }
  inG.close();

  ifstream inGM(file_in_GM.c_str(), ios::in);
  check_stream(inGM, file_in_GM);
  for(int ir=0; ir<NR; ir++) {
    inGM>>f;
    for(int im=0; im<NM; im++) {
      inGM>>b1>>b2>>d1>>d2>>d3;
      if(FLAG_IC==0)
        Vgm[(sub*NM+im)*NR+ir]=d1;
      else if(FLAG_IC==1)
        Vgm[(sub*NM+im)*NR+ir]=d2;
      else
        Vgm[(sub*NM+im)*NR+ir]=d3;
    }
  }
  inGM.close();

  //read the delta-sigma files
  ifstream inSHmm(file_in_SHmm.c_str(), ios::in);
  check_stream(inSHmm, file_in_SHmm);
  for(int ir=0; ir<NR; ir++) {
    inSHmm>>f>>d3; //these have only the IC_W correction
    if(FLAG_IC==2)
      SHmm[sub*NR+ir]=d3;
  }
  inSHmm.close();

  ifstream inSHgm(file_in_SHgm.c_str(), ios::in);
  check_stream(inSHgm, file_in_SHgm);
  for(int ir=0; ir<NR; ir++) {
    inSHgm>>f;
    for(int im=0; im<NM; im++) {
      inSHgm>>b1>>b2>>d3;
      if(FLAG_IC==2)
        SHgm[(sub*NM+im)*NR+ir]=d3;
    }
  }
  inSHgm.close();

}


void statistics::average_functions()

{

  if(Amm.size()!=0 || Agm.size() !=0 || Agg.size() !=0) {
    cerr<<"average projected functions already computed. "<<endl;
    return;
  }
  //matter-matter average and standard deviation

  vector<double> Tave(NR, 0.), Tvar(NR, 0.);

  for(int ir=0; ir<NR; ir++) {
    for(int is=0; is<Nsub; is++)
      Tave[ir]+=Vmm[is*NR+ir];
    Tave[ir]/=Nsub;
    for(int is=0; is<Nsub; is++)
      Tvar[ir]+=pow(Vmm[is*NR+ir]-Tave[ir], 2);
    Tvar[ir]=sqrt(Tvar[ir]/(Nsub-1)/Nsub);
  }
  for(int ir=0; ir<NR; ir++)
    Amm.push_back(Tave[ir]);
  for(int ir=0; ir<NR; ir++)
    Amm.push_back(Tvar[ir]);

  Tave.clear(); Tvar.clear();

  //galaxy-galaxy average and standard deviation

  Tave.resize(NR*NM); Tvar.resize(NR*NM); //check if new names make a difference here
  for(int i=0; i<NR*NM; i++) {
    Tave[i]=0.; Tvar[i]=0.;
  }

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++)
	Tave[im*NR+ir]+=Vgg[(is*NM+im)*NR+ir];
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++)
	Tvar[im*NR+ir]+=pow(Vgg[(is*NM+im)*NR+ir]-Tave[im*NR+ir], 2);
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Agg.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Agg.push_back(Tvar[im*NR+ir]);

  Tave.clear(); Tvar.clear();

  //galaxy-matter average and standard deviation

  Tave.resize(NR*NM); Tvar.resize(NR*NM);
  for(int i=0; i<NR*NM; i++) {
    Tave[i]=0.; Tvar[i]=0.;
  }

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++)
	Tave[im*NR+ir]+=Vgm[(is*NM+im)*NR+ir];
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++)
        Tvar[im*NR+ir]+=pow(Vgm[(is*NM+im)*NR+ir]-Tave[im*NR+ir], 2);
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Agm.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Agm.push_back(Tvar[im*NR+ir]);

  Tave.clear(); Tvar.clear();

  //do the SHmm average
  Tave.resize(NR); Tvar.resize(NR);
  for(int i=0; i<NR; i++) {
    Tave[i]=0.; Tvar[i]=0.;
  }

  for(int ir=0; ir<NR; ir++) {
    for(int is=0; is<Nsub; is++)
      Tave[ir]+=SHmm[is*NR+ir];
    Tave[ir]/=Nsub;
    for(int is=0; is<Nsub; is++)
      Tvar[ir]+=pow(SHmm[is*NR+ir]-Tave[ir], 2);
    Tvar[ir]=sqrt(Tvar[ir]/(Nsub-1)/Nsub);
  }
  for(int ir=0; ir<NR; ir++)
    ASHmm.push_back(Tave[ir]);
  for(int ir=0; ir<NR; ir++)
    ASHmm.push_back(Tvar[ir]);

  Tave.clear(); Tvar.clear();

  //do the SHgm average

  Tave.resize(NR*NM); Tvar.resize(NR*NM);
  for(int i=0; i<NR*NM; i++) {
    Tave[i]=0.; Tvar[i]=0.;
  }
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++)
        Tave[im*NR+ir]+=SHgm[(is*NM+im)*NR+ir];
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++)
        Tvar[im*NR+ir]+=pow(SHgm[(is*NM+im)*NR+ir]-Tave[im*NR+ir], 2);
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      ASHgm.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      ASHgm.push_back(Tvar[im*NR+ir]);

  Tave.clear(); Tvar.clear();
  
  cerr<<"done computing the average and standard deviations for the mm, gg, gm projected correlation functions."<<endl;
  cerr<<"done computing the average and standard deviations for the shear mm and gm."<<endl;

}


vector<double> statistics::average_bias_gg()

{

  cerr<<"COMPUTING THE GG BIAS: "<<endl;
  vector<double> Tave(NM*NR, 0.), Tvar(NM*NR, 0.), Vres;

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++) {
	/*if(Vmm[is*NR+ir]==0.) {
	  cerr<<"in bias_gg: zero value for w_mm at sub "<<is<<" rbin"<<ir<<endl;
	  continue;
	  }*/
	Tave[im*NR+ir]+=sqrt(abs(Vgg[(is*NM+im)*NR+ir]/Vmm[is*NR+ir]));
      }
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++) {
	//if(Vmm[is*NR+ir]==0.) 
	//continue;
	Tvar[im*NR+ir]+=pow((sqrt(abs(Vgg[(is*NM+im)*NR+ir]/Vmm[is*NR+ir]))-Tave[im*NR+ir]), 2);
      }
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }
  
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tvar[im*NR+ir]);

  Tave.clear(); Tvar.clear();

  return Vres;

}

vector<double> statistics::average_bias_gm()
  
{
  cerr<<"COMPUTING THE BIAS FROM GM: "<<endl;
  vector<double> Tave(NM*NR, 0.), Tvar(NM*NR, 0.), Vres;

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++) {
	/*if(Vmm[is*NR+ir]==0.) {
          cerr<<"in bias_gm: zero value for w_mm at sub "<<is<<" rbin"<<ir<<endl;
          continue;
	  }*/
        Tave[im*NR+ir]+=(abs(Vgm[(is*NM+im)*NR+ir]/Vmm[is*NR+ir]));
      }
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++) {
	//if(Vmm[is*NR+ir]==0.) 
	//continue;
        Tvar[im*NR+ir]+=pow((abs(Vgm[(is*NM+im)*NR+ir]/Vmm[is*NR+ir])-Tave[im*NR+ir]), 2);
      }
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tvar[im*NR+ir]);

  Tave.clear(); Tvar.clear();

  return Vres;

}


vector<double> statistics::correlation_coefficient()

{

  cerr<<"COMPUTING THE CORRELATION COEFFICIENT"<<endl;
  vector<double> Tave(NM*NR, 0.), Tvar(NM*NR, 0.), Vres;

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) {
      for(int is=0; is<Nsub; is++) {
	/*if(Vmm[is*NR+ir]==0. || Vgg[(is*NM+im)*NR+ir]==0 ) {
	  cerr<<"in corr_coeff: zero value for w_mm or w_gg for sub "<<is<<" rbin "<<ir<<" mag: "<<im<<endl;
	  continue;
	  }*/
	Tave[im*NR+ir]+=(Vgm[(is*NM+im)*NR+ir]/sqrt(abs(Vmm[is*NR+ir]*Vgg[(is*NM+im)*NR+ir])));
      }
      Tave[im*NR+ir]/=Nsub;
      for(int is=0; is<Nsub; is++) {
	//if(Vmm[is*NR+ir]==0. || Vgg[(is*NM+im)*NR+ir]==0 ) 
	//continue;
	Tvar[im*NR+ir]+=pow((Vgm[(is*NM+im)*NR+ir]/sqrt(abs(Vmm[is*NR+ir]*Vgg[(is*NM+im)*NR+ir]))-Tave[im*NR+ir]), 2);
      }
      Tvar[im*NR+ir]=sqrt(Tvar[im*NR+ir]/(Nsub-1)/Nsub);
    }

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tave[im*NR+ir]);
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      Vres.push_back(Tvar[im*NR+ir]);
  
  Tave.clear(); Tvar.clear();
  
  return Vres;

}



vector<SqDMatrix> statistics::cov_corr_function(string type, string type_return)

{

  //if type_return covariance, then return covariance matrix
  //if type_return correlation, then return correlation matrix
  if(type_return != "covariance" && type_return != "correlation") {
    cerr<<"only 2 possibilities for type-return: covariance and correlation"<<endl;
    exit(1);
  }
  average_functions();
  vector<SqDMatrix> VecCov, VecCorr;
  vector<double> V1, V2, A1, A2;
 
  if(type=="gggg") {
    V1=V2=Vgg; A1=A2=Agg; 
  }
  else if(type=="gmgm") {
    V1=V2=Vgm; A1=A2=Agm; 
  }
  else if(type=="gggm") {
    V1=Vgg; V2=Vgm; A1=Agg; A2=Agm; 
  }
  else if(type=="ggmm") {
    V1=Vgg; V2=Vmm; A1=Agg; A2=Amm;
  }
  else if(type=="gmmm") {
    V1=Vgm; V2=Vmm; A1=Agm; A2=Amm;
  }
  else if(type=="mmmm") {
    V1=V2=Vmm; A1=A2=Amm;
  }
  else if(type=="sgmsgm") {
    V1=V2=SHgm; A1=A2=ASHgm;
  }
  else if(type=="ggsgm") {
    V1=Vgg; V2=SHgm; A1=Agg; A2=ASHgm;
  }
  else if(type=="smmsmm") {
    V1=V2=SHmm; A1=A2=ASHmm;
  }
  else if(type=="mmsmm") {
    V1=Vmm; V2=SHmm; A1=Amm; A2=ASHmm;
  }
  else {
    cerr<<"unknown covariance type. Allowed types are: "<<endl;
    cerr<<"gggg, gmgm, gggm, ggmm, gmmm, mmmm, sgmsgm, ggsgm, smmsmm, mmsmm."<<endl;
    exit(1);
  } 
  //NOTE: mmmm, smmsmm, mmsmm, ggmm, gmmm are basically for testing
  if(type=="gggg" || type=="gmgm" || type=="gggm" || type=="sgmsgm" || type=="ggsgm" ) {
    for(int im=0; im<NM; im++) {
      SqDMatrix Cov(NR, 0.); SqDMatrix Corr(NR);
      for(int ir1=0; ir1<NR; ir1++)
	for(int ir2=0; ir2<NR; ir2++) {
	  for(int is=0; is<Nsub; is++) 
	    Cov(ir1, ir2)+=(V1[(is*NM+im)*NR+ir1]-A1[im*NR+ir1]) 
	      * (V2[(is*NM+im)*NR+ir2]-A2[im*NR+ir2]);
	  Cov(ir1, ir2)/=(Nsub-1); //unbiased estimator
	  Cov(ir1, ir2)/=Nsub; //covariance on the mean
	  Corr(ir1, ir2)=Cov(ir1, ir2)/A1[NM*NR+im*NR+ir1]/A2[NM*NR+im*NR+ir2];
	}
      VecCov.push_back(Cov);
      VecCorr.push_back(Corr);
      if(im==0)
	for(int ir1=0; ir1<NR; ir1++)
	  for(int ir2=0; ir2<NR; ir2++) 
	    cerr<<ir1<<' '<<ir2<<' '<<Cov(ir1, ir2)<<' '<<Corr(ir1, ir2)<<endl;
    }
  }

  else if(type=="ggmm" || type=="gmmm") {
    for(int im=0; im<NM; im++) {
      SqDMatrix Cov(NR, 0.); SqDMatrix Corr(NR, 0.);
      for(int ir1=0; ir1<NR; ir1++)
	for(int ir2=0; ir2<NR; ir2++) {
	  for(int is=0; is<Nsub; is++) 
	    Cov(ir1, ir2)+=(V1[(is*NM+im)*NR+ir1]-A1[im*NR+ir1]) 
	      * (V2[is*NR+ir2]-A2[ir2]);
	  Cov(ir1, ir2)/=(Nsub-1); //unbiased estimator
	  Cov(ir1, ir2)/=Nsub; //covariance on the mean
	  Corr(ir1, ir2)=Cov(ir1, ir2)/A1[NM*NR+im*NR+ir1]/A2[NR+ir2];
	}
      VecCov.push_back(Cov);
      VecCorr.push_back(Corr);
    }
  }
  
  else if(type=="mmmm" || type=="smmsmm" || type=="mmsmm") {
    SqDMatrix Cov(NR, 0.); SqDMatrix Corr(NR, 0.);
    for(int ir1=0; ir1<NR; ir1++)
      for(int ir2=0; ir2<NR; ir2++) {
	for(int is=0; is<Nsub; is++) 
	  Cov(ir1, ir2)+=(V1[is*NR+ir1]-A1[ir1]) * (V2[is*NR+ir2]-A2[ir2]);
	Cov(ir1, ir2)/=(Nsub-1); //unbiased estimator
	Cov(ir1, ir2)/=Nsub; //covariance on the mean
	Corr(ir1, ir2)=Cov(ir1, ir2)/A1[NR+ir1]/A2[NR+ir2];
      }
    VecCov.push_back(Cov);
    VecCorr.push_back(Corr);
  }


  //return the covariance or correlation matrix

  if(type_return=="covariance")
    return VecCov;
  else if(type_return=="correlation") 
    return VecCorr;
  
}


 
void statistics::write_average_functions(string file, string type)

{

  if(Amm.size()==0) 
    average_functions();
  vector<double> A;
  if(type=="mm") A=Amm;
  else if(type=="gg") A=Agg;
  else if(type=="gm") A=Agm;
  else if(type=="smm") A=ASHmm;
  else if(type=="sgm") A=ASHgm;
  else {
    cerr<<"Wrong type. Allowed types are: "<<endl;
    cerr<<"mm, gg, gm, smm, sgm."<<endl;
    exit(1);
  }
  
  ofstream out(file.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    if(type=="mm" || type=="smm") 
      out<<binsR[ir]<<' '<<A[ir]<<' '<<A[NR+ir]<<endl;
    else if(type=="gg" || type=="gm" || type=="sgm") {
      out<<binsR[ir];
      for(int im=0; im<NM; im++)
        out<<' '<<A[im*NR+ir]<<' '<<A[NR*NM+im*NR+ir];
      out<<endl;
    }
  }

  out.close();
  cerr<<"done writing output file " + file<<endl;

}

void statistics::write_bias(string file[2])

{

  vector<double> V1=average_bias_gg();
  ofstream out0(file[0].c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    out0<<binsR[ir];
    for(int im=0; im<NM; im++)
      out0<<' '<<V1[im*NR+ir]<<' '<<V1[NR*NM+im*NR+ir];
    out0<<endl;
  }
  out0.close();

  vector<double> V2=average_bias_gm();
  ofstream out1(file[1].c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    out1<<binsR[ir];
    for(int im=0; im<NM; im++)
      out1<<' '<<V2[im*NR+ir]<<' '<<V2[NR*NM+im*NR+ir];
    out1<<endl;
  }
  out1.close();

}


void statistics::write_correlation_coefficient(string file)

{

  vector<double> V=correlation_coefficient();
  ofstream out(file.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    out<<binsR[ir];
    for(int im=0; im<NM; im++)
      out<<' '<<V[im*NR+ir]<<' '<<V[NR*NM+im*NR+ir];
    out<<endl;
  }
  out.close();

}



void statistics::write_cov_corr(vector<string> file, string type, string type_return)

{

  int AL=2; //eliminate last 2 radial bins, they are probably not safe for Lproj=100 Mpc/h
  int NR2=NR-AL;
  vector<SqDMatrix> VecC;
  VecC=cov_corr_function(type, type_return);

  if(VecC.size() != file.size()) {
    cerr<<"In write_cov_corr: mismatch between size of output files and size of covariance matrix vector"<<endl;
    exit(1);
  } 

  for(int im=0; im<file.size(); im++) {
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

}

