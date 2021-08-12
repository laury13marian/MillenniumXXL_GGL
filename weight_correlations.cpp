#include "general.h"
#include <vector>
#include <cmath>

int main(int argc, char *argv[])

{

  if(argc!=4) { //3 here usually
    cerr<<"introduce the subcube number and zmax, and task number"<<endl;
    exit(1);
  }

  const int sub=atoi(argv[1]);
  const double zMAX=atof(argv[2]);
  const int TASK=atof(argv[3]);
  const int snap=54;
  const int NZ=10;
  const int NR=30; //CAREFUL HERE!!!!!!!!!!!!!!!!!!!!
  //const int NR=NRini-2;//CAREFUL HERE!!!!!!!!!!!!!!!!!!!! 
  const int nsampM=4000; //3000;
  const int ntasks=32; //64;
  const long long NRAN=1000000;

  const int Nmags=10;
  const int NM=Nmags-1; //number of magnitude/mass bins        
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;

  ostringstream osSC, osSNAP, osNZ, osNR, osNSAMP, osT, osNRAN, osNM;
  osSC<<sub; osSNAP<<0; osSNAP<<snap; osNZ<<NZ; osNR<<NR; osNM<<NM; 
  osNSAMP<<nsampM; osT<<ntasks; osNRAN<<NRAN;

  /*const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/TESTS_NEW_BINS/CORR_FUNC/"; 
  string file_in_MM = path_gen + "MM_proj_corr_snap_" + osSNAP.str() + "_NRAN_" + osNRAN.str() + "_binsR_" + osNB.str() 
  + "_binsZ_" + osNZ.str() + "_nsamp_" + osM.str() + "_ntasks_" + osT.str() + "_zmax_" + argv[2] + ".dat";*/

  const string path_gen = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/";

  const string file_in_MM = path_gen + "CORR_MM/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] 
    + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" 
    + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" + argv[3] + ".dat";

  const string file_out_MM = path_gen + "CORR_MM/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_"
    + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" + argv[3] + ".dat";

  const string file_in_GG = path_gen + "CORR_GG/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" 
    + osNRAN.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" + argv[3] + ".dat";

  const string file_out_GG = path_gen + "CORR_GG/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
    + osNRAN.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" + argv[3] + ".dat";
 
  const string file_in_GM = path_gen + "CORR_GM/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" 
    + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" 
    + argv[3] + ".dat";

  const string file_out_GM = path_gen + "CORR_GM/ONE_TASK/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
    + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + "_one_task_" 
    + argv[3] + ".dat";


  /*const string file_in_GG = path_gen + "CORR_GG/CF_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" 
    + osNRAN.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + ".dat";
  
  const string file_out_GG = path_gen + "CORR_GG/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
    + osNRAN.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + ".dat";

  const string file_in_GM = path_gen + "CORR_GM/CF_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" 
    + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + ".dat";
  
  const string file_out_GM = path_gen + "CORR_GM/WEIGHTED/CFw_snap_" + osSNAP.str() + "_sub_" + argv[1]
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_"
    + osNRAN.str() + "_nsamp_" + osNSAMP.str() + "_zmax_" + argv[2] + "_NT_" + osT.str() + ".dat";*/

  vector<float> binsR(NR), binsZ(NZ);
  vector<double> WP(NZ, 1.);
  double Tmm[NR][NZ], TmmD[NR][NZ], TmmW[NR][NZ];
  vector<double> Cmm(NR, 0.), CmmD(NR, 0.), CmmW(NR, 0.);
  double Tgg[NM][NR][NZ], TggD[NM][NR][NZ], TggW[NM][NR][NZ];
  vector<double> Cgg(NM*NR, 0.), CggD(NM*NR, 0.), CggW(NM*NR, 0.);
  double Tgm[NM][NR][NZ], TgmD[NM][NR][NZ], TgmW[NM][NR][NZ];
  vector<double> Cgm(NM*NR, 0.), CgmD(NM*NR, 0.), CgmW(NM*NR, 0.);

  ifstream inM(file_in_MM.c_str(), ios::in);
  check_stream(inM, file_in_MM);
 
  for(int ir=0; ir<NR; ir++)   //CHECK FILE STRUCTURE HERE
    for(int jz=0; jz<NZ; jz++) 
      inM>>binsR[ir]>>binsZ[jz]>>Tmm[ir][jz]>>TmmD[ir][jz]>>TmmW[ir][jz];
  inM.close();
  cerr<<"closing file " + file_in_MM<<endl;
  cerr<<"rperp bins: "<<endl;
  for(int ir=0; ir<NR; ir++)
    cerr<<ir<<' '<<binsR[ir]<<endl;
  cerr<<endl;
  cerr<<"z bins: "<<endl;
  for(int iz=0; iz<NZ; iz++)
    cerr<<iz<<' '<<binsZ[iz]<<endl;
  cerr<<endl;
  
  //assumes dz is constant for all bins
  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) {  
      Cmm[ir]+=(Tmm[ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.); //equidistant bins
      CmmD[ir]+=(TmmD[ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.);
      CmmW[ir]+=(TmmW[ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.);
    }

  ifstream inG(file_in_GG.c_str(), ios::in);
  check_stream(inG, file_in_GG);
  float f1, f2; 
  double d1, d2;
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      inG>>f1>>f2;
      for(int im=0; im<NM; im++)
	inG>>d1>>d2>>Tgg[im][ir][jz]>>TggD[im][ir][jz]>>TggW[im][ir][jz];
    }
  inG.close();
  cerr<<"closing file " + file_in_GG<<endl;

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      for(int jz=0; jz<NZ; jz++) {
	Cgg[im*NR+ir]+=(Tgg[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.); //equidistant bins 
	CggD[im*NR+ir]+=(TggD[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.); 
	CggW[im*NR+ir]+=(TggW[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.); 
      }

  ifstream inGM(file_in_GM.c_str(), ios::in);
  check_stream(inGM, file_in_GM);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      inGM>>f1>>f2;
      for(int im=0; im<NM; im++)
        inGM>>d1>>d2>>Tgm[im][ir][jz]>>TgmD[im][ir][jz]>>TgmW[im][ir][jz];
    }
  inGM.close();
  cerr<<"closing file " + file_in_GM<<endl;

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++)
      for(int jz=0; jz<NZ; jz++) {
        Cgm[im*NR+ir]+=(Tgm[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.); //equidistant bins   
	CgmD[im*NR+ir]+=(TgmD[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.);
	CgmW[im*NR+ir]+=(TgmW[im][ir][jz] * WP[jz] * (binsZ[1]-binsZ[0]) * 2.);
      } 
  //assumes dz is constant for all bins
    
  /*ofstream outM(file_out_MM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    outM<<binsR[ir]<<' '<<Cmm[ir]<<' '<<CmmD[ir]<<' '<<CmmW[ir]<<endl;
  outM.close();
  cerr<<"done writing weighted file " + file_out_MM<<endl<<endl;*/
  
  ofstream outG(file_out_GG.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    outG<<binsR[ir];
    for(int im=0; im<NM; im++)
      outG<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<Cgg[im*NR+ir]<<' '<<CggD[im*NR+ir]<<' '<<CggW[im*NR+ir];
    outG<<endl;
  }
  outG.close();
  cerr<<"done writing weighted file " + file_out_GG<<endl<<endl;

  ofstream outGM(file_out_GM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    outGM<<binsR[ir];
    for(int im=0; im<NM; im++)
      outGM<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<Cgm[im*NR+ir]<<' '<<CgmD[im*NR+ir]<<' '<<CgmW[im*NR+ir];
    outGM<<endl;
  }
  outGM.close();
  cerr<<"done writing weighted file " + file_out_GM<<endl<<endl;
  

  return 0;

}




