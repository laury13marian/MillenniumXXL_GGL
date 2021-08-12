//enforces the integral constraint and combines the z-bins for the JK sample
//we compute just the former W case, i.e. exactly the IC according the Landy-Szalay
//each JK sample is done at a time, taken as argument through the main()
#include "general.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 3) {
    cerr<<"must input JK number and task number"<<endl;
    exit(10);
  }

  const int JKsample = atoi(argv[1]);
  const int TASK = atoi(argv[2]);
  const int ntasks=32; //64;
  const int sub = 100; //CHECK THIS
  const int snap=54; 
  const int nsampM=4000;
  const float zMAX=100.;
  const double rbox=500.;//3000.;
  const float rMAX=100.;//110.;//30.;                                              
  const float rMIN=0.01;
  const float zMIN=0.;
  const int NR=30; //20; //10;
  const int NZ=10;
  const int NJK = 64;
  float rb[3][NR]; //radial bin vector
  float zb[3][NZ]; //z-bin vector 
  const size_t NRAN=1000000; //3000000;     
  const int dim=3;
  const int Nmags=5;
  const int NM=Nmags-1; //number of magnitude/mass bins  
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; 
  const int iMAG=1;

  //random-random; 
  double npairsRR[NR][NZ];

  //for matter-matter
  double ratioMM[NR][NZ], ncorrMM[NR][NZ];
  double npairsMM[NR][NZ], npairsRM[NR][NZ];

  //for galaxy-galaxy
  double ratioGG[NM][NR][NZ], ncorrGG[NM][NR][NZ];
  double npairsGG[NM][NR][NZ], npairsRG[NM][NR][NZ];

  // galaxy-matter
  double ratioGM[NM][NR][NZ], ncorrGM[NM][NR][NZ];
  double npairsGM[NM][NR][NZ];

  ostringstream osNR, osNZ, osNM, osSNAP, osSAMPM, osNTASK, osZMAX, osNRAN, osJK, osSUB;
  osNR<<NR; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; osNRAN<<NRAN;
  osNTASK<<ntasks; osNZ<<NZ; osNM<<NM; osZMAX<<zMAX; osSUB<<sub;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/"; 

  const string file_in_RR = path_gen + 
    "RANDOM_CATALOGUE/JACK_KNIFE/pair_counts_NRAN_" + osNRAN.str() + "_JKsample_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";

  const string file_in_MM_RM = path_gen + "NEW_BINS/JACK_KNIFE/PAIR_COUNTS_MM_RM/PC_snap_" + osSNAP.str() 
    + "_sub_" + osSUB.str() + "_JKsample_" + argv[1] + "_binsR_" + osNR.str() + "_NRAN_" + osNRAN.str() 
    + "_nsamp_" + osSAMPM.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";
    
  const string file_in_GG = path_gen + "NEW_BINS/JACK_KNIFE/PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() 
    + "_sub_" + osSUB.str() + "_JKsample_" + argv[1] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str()
    + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";
  
  const string file_in_RG = path_gen + "NEW_BINS/JACK_KNIFE/PAIR_COUNTS_GG_RG/RG_snap_" + osSNAP.str() 
    + "_sub_" + osSUB.str() + "_JKsample_" + argv[1] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() 
    + "_NRAN_" + osNRAN.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";
  
  const string file_in_GM = path_gen + "NEW_BINS/JACK_KNIFE/PAIR_COUNTS_GM/PC_snap_" + osSNAP.str() 
    + "_sub_" + osSUB.str() + "_JKsample_" + argv[1] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() 
    + "_nsamp_" + osSAMPM.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";
 
  //apply first the L-S integral constraint

  for(int j=0; j<NR; j++) 
    for(int jz=0; jz<NZ; jz++) {
      ncorrMM[j][jz]=0.;
      for(int im=0; im<NM; im++) {
	ncorrGG[im][j][jz]=0.; ncorrGM[im][j][jz]=0.;
      }
    }

  //get the bins and the volume element

  vector<particle> VpartM;
  vector<float> pos1(3,100.), pos2(3,150.), pos3(5,39.);
  VpartM.push_back(particle(pos1));  VpartM.push_back(particle(pos2));  VpartM.push_back(particle(pos3));
  
 const vector<float> vecPER(dim, 0.);
  projected_correlations::BinType bt = projected_correlations::LOG10;
  bool flag_per=false; //true;
  cerr<<"box size: "<<rbox<<endl;
  build_tree OtreeM(rbox, VpartM);
  build_tree *OM1, *OM2;
  OM1=&OtreeM; OM2=&OtreeM;
  projected_correlations OcorrMM(OM1, OM2, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER);
  vector<float> bin_R=OcorrMM.return_bins_rperp();
  vector<float> bin_Z=OcorrMM.return_bins_z();
  for(int i=0; i<NR; i++) {
    rb[0][i]=bin_R[i];
    rb[1][i]=bin_R[2*NR+i];
    rb[2][i]=bin_R[NR+i];
  }
  
  for(int i=0; i<NZ; i++) {
    zb[0][i]=bin_Z[i]; //lower edge
    zb[1][i]=bin_Z[2*NZ+i]; //mean
    zb[2][i]=bin_Z[NZ+i]; //upper edge
  }
  
  vector<double> dV=OcorrMM.volume_bins();
  //read in data
  
  ifstream inRR(file_in_RR.c_str(), ios::in);
  check_stream(inRR, file_in_RR);
  float x1, x2, x3, x4, x5, x6;
  for(int ir=0; ir<NR; ir++)
    for(int iz=0; iz<NZ; iz++)
      inRR>>x1>>x2>>x3>>x4>>x5>>x6>>npairsRR[ir][iz];
  inRR.close();

  ifstream inMM(file_in_MM_RM.c_str(), ios::in);
  check_stream(inMM, file_in_MM_RM);
  for(int ir=0; ir<NR; ir++)
    for(int iz=0; iz<NZ; iz++)
      inMM>>x1>>x2>>npairsMM[ir][iz]>>npairsRM[ir][iz];
  inMM.close();

  ifstream inGG(file_in_GG.c_str(), ios::in);
  check_stream(inGG, file_in_GG);
  double d1, d2;
  for(int ir=0; ir<NR; ir++)
    for(int iz=0; iz<NZ; iz++) {
      inGG>>x1>>x2;
      for(int im=0; im<NM; im++)
	inGG>>d1>>d2>>npairsGG[im][ir][iz];
    }
  inGG.close();

  ifstream inRG(file_in_RG.c_str(), ios::in);
  check_stream(inRG, file_in_RG);
  for(int ir=0; ir<NR; ir++)
    for(int iz=0; iz<NZ; iz++) {
      inRG>>x1>>x2;
      for(int im=0; im<NM; im++)
	inRG>>d1>>d2>>npairsRG[im][ir][iz];
    }
  inRG.close();

  ifstream inGM(file_in_GM.c_str(), ios::in);
  check_stream(inGM, file_in_GM);
  for(int ir=0; ir<NR; ir++)
    for(int iz=0; iz<NZ; iz++) {
      inGM>>x1>>x2;
      for(int im=0; im<NM; im++)
	inGM>>d1>>d2>>npairsGM[im][ir][iz];
    }
  inGM.close();

  double NormR, NormM, CM; 
  NormR=0.; NormM=0.;
  vector<double> NormG(NM, 0.), CG(NM);
  vector<double> NormGM(NM, 0.), CGM(NM);

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) {
      ratioMM[ir][jz]=npairsMM[ir][jz]/npairsRR[ir][jz]-2.*npairsRM[ir][jz]/npairsRR[ir][jz]+1.;
      for(int im=0; im<NM; im++) {
	ratioGG[im][ir][jz]=npairsGG[im][ir][jz]/npairsRR[ir][jz]-2.*npairsRG[im][ir][jz]/npairsRR[ir][jz]+1.;
	ratioGM[im][ir][jz]=npairsGM[im][ir][jz]/npairsRR[ir][jz]-npairsRG[im][ir][jz]/npairsRR[ir][jz]
	  -npairsRM[ir][jz]/npairsRR[ir][jz] + 1.;
      }
    }

 //compute the integral constraint

  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      NormR+=(dV[ir*NZ+jz]*npairsRR[ir][jz]);
      NormM+=(dV[ir*NZ+jz]*ratioMM[ir][jz]*npairsRR[ir][jz]); //w_vol for MM
      }
  CM=NormM/NormR; //w_vol from the whole w_LS

  for(int im=0; im<NM; im++) {
    for(int ir=0; ir<NR; ir++)
      for(int jz=0; jz<NZ; jz++) {
	NormG[im]+=(dV[ir*NZ+jz]*ratioGG[im][ir][jz]*npairsRR[ir][jz]); //w_vol for gal-gal
	NormGM[im]+=(dV[ir*NZ+jz]*ratioGM[im][ir][jz]*npairsRR[ir][jz]); //w_vol for gal-matter
      }
    CG[im]=NormG[im]/NormR; //w_vol from the whole w_LS for gal-gal
    CGM[im]=NormGM[im]/NormR; //w_vol from the whole w_LS for gal-matter
  }
  
  //compute the Landy-Szalay estimator, with the L-S integral constraint (former W)

  //matter-matter
  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) 
   	ncorrMM[ir][jz]=(1.+ratioMM[ir][jz])/(1.+CM)-1.;

  //galaxy-galaxy and galaxy-matter
  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) 
      for(int jz=0; jz<NZ; jz++) {
	ncorrGG[im][ir][jz]=(1.+ratioGG[im][ir][jz])/(1.+CG[im])-1.;
	ncorrGM[im][ir][jz]=(1.+ratioGM[im][ir][jz])/(1.+CGM[im])-1.;
      }

  //now combine the z-bins to form projected correlation fucntions

  const string path_gen_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/JACK_KNIFE/";

  const string file_out_MM = path_gen_out + "CORR_MM/ONE_TASK/CFw_snap_" + osSNAP.str() + "_sub_" + osSUB.str() + "_JKsample_"
    + argv[1] + "_task_" + argv[2] + "_binsR_" + osNR.str() + "_nsamp_" + osSAMPM.str() + "_NT_" + osNTASK.str() + ".dat";
  
  const string file_out_GG = path_gen_out + "CORR_GG/ONE_TASK/CFw_snap_" + osSNAP.str() + "_sub_" + osSUB.str() + "_JKsample_"
    + argv[1] + "_task_" + argv[2] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str()  + "_NT_" + osNTASK.str() + ".dat";
  
  const string file_out_GM = path_gen_out + "CORR_GM/ONE_TASK/CFw_snap_" + osSNAP.str() + "_sub_" + osSUB.str() + "_JKsample_"
    + argv[1] + "_task_" + argv[2] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_nsamp_" + osSAMPM.str() 
    + "_NT_" + osNTASK.str() + ".dat";
  
  vector<double> WP(NZ, 1.);
  vector<double> wMM(NR, 0.), wGG(NM*NR, 0.), wGM(NM*NR, 0.);
  
 //assumes dz is constant for all bins

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) 
      wMM[ir]+=(ncorrMM[ir][jz] * WP[jz] * (zb[0][1]-zb[0][0]) * 2.); //equidistant bins

  for(int ir=0; ir<NR; ir++) 
    for(int im=0; im<NM; im++) 
      for(int jz=0; jz<NZ; jz++) {
	wGG[im*NR+ir]+=(ncorrGG[im][ir][jz] * WP[jz] * (zb[0][1]-zb[0][0]) * 2.); //equidistant bins
	wGM[im*NR+ir]+=(ncorrGM[im][ir][jz] * WP[jz] * (zb[0][1]-zb[0][0]) * 2.);       
      }

  //write weighted correlations

  ofstream outM(file_out_MM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    outM<<rb[1][ir]<<' '<<wMM[ir]<<endl;
  outM.close();
  cerr<<"done writing weighted file " + file_out_MM<<endl<<endl;
  
  ofstream outG(file_out_GG.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    outG<<rb[1][ir];
    for(int im=0; im<NM; im++)
      outG<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<wGG[im*NR+ir];
    outG<<endl;
  }
  outG.close();
  cerr<<"done writing weighted file " + file_out_GG<<endl<<endl;
  
  ofstream outGM(file_out_GM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++) {
    outGM<<rb[1][ir];
    for(int im=0; im<NM; im++)
      outGM<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<wGM[im*NR+ir];
    outGM<<endl;
  }
  outGM.close();
  cerr<<"done writing weighted file " + file_out_GM<<endl<<endl;
  
  return 0;
  
}








  
