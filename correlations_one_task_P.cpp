//driver to process the subcube MM, GG, GM data

#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 3) {
    cerr<<"must input subcube and task numbers."<<endl;
    exit(10);
  }

  const int sub = atoi(argv[1]);
  const int TASK = atoi(argv[2]);
  const int ntasks=32; //64;
  const int snap=54; 
  const int nsampM=4000;
  const float zMAX=100.;
  const double rbox=500.;//3000.;
  const float rMAX=100.;//110.;//30.;                                              
  const float rMIN=0.01;
  const float zMIN=0.;
  const int NR=30; //20; //10;
  const int NZ=10;
  float rb[3][NR]; //radial bin vector
  float zb[3][NZ]; //z-bin vector 
  const size_t NRAN=1000000; //3000000;     

  const int Nmags=10;
  const int NM=Nmags-1; //number of magnitude/mass bins  
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;
  const int iMAG=1;

  //for matter-matter

  double ratioMM[NR][NZ], ncorrMM[NR][NZ], ncorrMMD[NR][NZ], ncorrMMW[NR][NZ];
  double npairsMM[NR][NZ], npairsRM[NR][NZ], npairsRR[NR][NZ];

  //for galaxy-galaxy

  double ratioGG[NM][NR][NZ], ncorrGG[NM][NR][NZ], ncorrGGD[NM][NR][NZ], ncorrGGW[NM][NR][NZ];
  double npairsGG[NM][NR][NZ], npairsRG[NM][NR][NZ];

  // galaxy-matter

  double ratioGM[NM][NR][NZ], ncorrGM[NM][NR][NZ], ncorrGMD[NM][NR][NZ], ncorrGMW[NM][NR][NZ];
  double npairsGM[NM][NR][NZ];
    
  const int dim=3;
  
  ostringstream osNR, osNZ, osNM, osSNAP, osSAMPM, osNTASK, osZMAX, osNRAN;
  osNR<<NR; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; osNRAN<<NRAN;
  osNTASK<<ntasks; osNZ<<NZ; osNM<<NM; osZMAX<<zMAX;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/"; 
  
  const string file_in_RR = path_gen + "RANDOM_CATALOGUE/NEW_BINS/pair_counts_NRAN_" + osNRAN.str() + "_binsR_" + osNR.str()
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";

  const string file_in_MM_RM = path_gen + "NEW_BINS/PAIR_COUNTS_MM_RM/PC_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_task_"  + argv[2] + ".dat";
  
  string file_in_GG;
  if(sub<7)
    file_in_GG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() 
      + "_NT_64_task_" + argv[2] + ".dat";
  else
    file_in_GG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() 
      + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";

  const string file_in_RG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/RG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";

  const string file_in_GM = path_gen + "NEW_BINS/PAIR_COUNTS_GM/PC_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsM_" 
    + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";

 
  const string file_out_MM = path_gen + "NEW_BINS/CORR_MM/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_one_task_" + argv[2] + ".dat";

  const string file_out_GG = path_gen + "NEW_BINS/CORR_GG/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] 
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_one_task_" + argv[2] + ".dat";

  const string file_out_GM = path_gen + "NEW_BINS/CORR_GM/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsM_" 
    + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_one_task_" + argv[2] + ".dat";
  
 
  //get the bins and the volume element
  
  /*vector<string> path_data(4);
  path_data[0] =
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_";
  const string SPM = "PARTICLE";
  const signed long rng_seed_gen = -1381718569;
  const signed long rng_seed = rng_seed_gen;
  cerr<<"READING PARTICLE DATA"<<endl;
  //load the 0 subcube for access to various functions related to the binning
  load_SScube OBM(0, path_data, SPM, nsampM, rng_seed);
  vector<particle> VpartM=OBM.return_formatted_data();
  cerr<<"size of formatted particles: "<<VpartM.size()<<endl;*/
  //one doesn't have to load subcubes, because it's too slow
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
    zb[0][i]=bin_Z[i];
    zb[1][i]=bin_Z[2*NZ+i];
    zb[2][i]=bin_Z[NZ+i];
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
  

  double NormR(0.), NormMD(0.), NormMW(0.), CMD, CMW;
  vector<double> NormGD(NM, 0.), NormGW(NM, 0.), CGD(NM), CGW(NM);
  vector<double> NormGMD(NM, 0.), NormGMW(NM, 0.), CGMD(NM), CGMW(NM);

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) 
      ratioMM[ir][jz]=npairsMM[ir][jz]/npairsRR[ir][jz]-2.*npairsRM[ir][jz]/npairsRR[ir][jz]+1.;

  for(int im=0; im<NM; im++) 
    for(int ir=0; ir<NR; ir++) 
      for(int jz=0; jz<NZ; jz++) {
	ratioGG[im][ir][jz]=npairsGG[im][ir][jz]/npairsRR[ir][jz]-2.*npairsRG[im][ir][jz]/npairsRR[ir][jz]+1.;
	ratioGM[im][ir][jz]=npairsGM[im][ir][jz]/npairsRR[ir][jz]-npairsRG[im][ir][jz]/npairsRR[ir][jz]
	  -npairsRM[ir][jz]/npairsRR[ir][jz] + 1.;
      }

  //compute the integral constraint  

  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      NormR+=(dV[ir*NZ+jz]*npairsRR[ir][jz]);
      NormMD+=(dV[ir*NZ+jz]*npairsMM[ir][jz]);
      NormMW+=(dV[ir*NZ+jz]*ratioMM[ir][jz]*npairsRR[ir][jz]); //w_vol for MM
    }
  CMD=NormMD/NormR; //1+w_vol from DD
  CMW=NormMW/NormR; //w_vol from the whole w_LS


  //for(int np=0; np<ntasks; np++) 
  for(int im=0; im<NM; im++) 
    for(int ir=0; ir<NR; ir++)
      for(int jz=0; jz<NZ; jz++) {
	NormGD[im]+=(dV[ir*NZ+jz]*npairsGG[im][ir][jz]);
	NormGW[im]+=(dV[ir*NZ+jz]*ratioGG[im][ir][jz]*npairsRR[ir][jz]); //w_vol for gal-gal
	NormGMD[im]+=(dV[ir*NZ+jz]*npairsGM[im][ir][jz]);
	NormGMW[im]+=(dV[ir*NZ+jz]*ratioGM[im][ir][jz]*npairsRR[ir][jz]); //w_vol for gal-matter
      }
  for(int im=0; im<NM; im++) {
    CGD[im]=NormGD[im]/NormR; //1+w_vol from DD for gal-gal
    CGW[im]=NormGW[im]/NormR; //w_vol from the whole w_LS for gal-gal
    CGMD[im]=NormGMD[im]/NormR; //1+w_vol from DD for gal-matter
    CGMW[im]=NormGMW[im]/NormR; //w_vol from the whole w_LS for gal-matter
  }

  //compute the Landy-Szalay estimator, with and without the integral constraint

  //matter-matter

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) {
      ncorrMM[ir][jz]=ratioMM[ir][jz]; //without the IC
      ncorrMMD[ir][jz]=(npairsMM[ir][jz]/npairsRR[ir][jz]/CMD-2.*npairsRM[ir][jz]/npairsRR[ir][jz]+1.);
      ncorrMMW[ir][jz]=(1.+ratioMM[ir][jz])/(1.+CMW);
      ncorrMMW[ir][jz]-=1.;
    }
  
  //galaxy-galaxy

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) 
      for(int jz=0; jz<NZ; jz++) {
	//for(int np=0; np<ntasks; np++) {
	ncorrGG[im][ir][jz]=ratioGG[im][ir][jz]; //no IC
	ncorrGM[im][ir][jz]=ratioGM[im][ir][jz]; //no IC
	ncorrGGD[im][ir][jz]=npairsGG[im][ir][jz]/npairsRR[ir][jz]/CGD[im]-2.*npairsRG[im][ir][jz]/npairsRR[ir][jz]+1.;
	ncorrGMD[im][ir][jz]=npairsGM[im][ir][jz]/npairsRR[ir][jz]/CGMD[im]-npairsRG[im][ir][jz]/npairsRR[ir][jz]
	  -npairsRM[ir][jz]/npairsRR[ir][jz]+1.;
	ncorrGGW[im][ir][jz]=(1.+ratioGG[im][ir][jz])/(1.+CGW[im]);
	ncorrGMW[im][ir][jz]=(1.+ratioGM[im][ir][jz])/(1.+CGMW[im]);
	ncorrGGW[im][ir][jz]-=1.;
	ncorrGMW[im][ir][jz]-=1.;
      }  

  //write the output files
  //matter-matter
  
  /*ofstream outM(file_out_MM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) 
      outM<<rb[1][ir]<<' '<<zb[1][jz]<<' '<<ncorrMM[ir][jz]<<' '<<ncorrMMD[ir][jz]<<' '<<ncorrMMW[ir][jz]<<endl;
  outM.close();
  cerr<<"\n done writing the MM output file " + file_out_MM<<endl;*/

  //galaxy-galaxy

  ofstream outGG(file_out_GG.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      outGG<<rb[1][ir]<<' '<<zb[1][jz];
      for(int im=0; im<NM; im++)
	outGG<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<ncorrGG[im][ir][jz]<<' '<<ncorrGGD[im][ir][jz]<<' '<<ncorrGGW[im][ir][jz];
      outGG<<endl;
    }
  outGG.close();
  cerr<<"\n done writing the GG output file " + file_out_GG<<endl;

  //galaxy-matter

  ofstream outGM(file_out_GM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      outGM<<rb[1][ir]<<' '<<zb[1][jz];
      for(int im=0; im<NM; im++)
	outGM<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<ncorrGM[im][ir][jz]<<' '<<ncorrGMD[im][ir][jz]<<' '<<ncorrGMW[im][ir][jz];
      outGM<<endl;
    }
  outGM.close();
  cerr<<"\n done writing the GM output file " + file_out_GM<<endl;
  
      
  return 0;
  
}

 

