//driver to process the subcube MM, GG, GM data

#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  //subcube number 0-215; task number 0-31;  
  if(argc != 3) {
    cerr<<"must input subcube number and task number"<<endl;
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
  /*vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;
  const int iMAG=1;*/

  //for matter-matter

  double ratioMM[NR][NZ], ncorrMM[NR][NZ], ncorrMMD[NR][NZ], ncorrMMW[NR][NZ];
  double npairsMM[NR][NZ], npairsRM[NR][NZ], npairsRR[NR][NZ];
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
    + "_NT_" + osNTASK.str() + "_task_" + argv[2] + ".dat";
  
  const string file_out_MM = path_gen + "NEW_BINS/CORR_MM/ONE_TASK/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_one_task_" + argv[2] + ".dat";

  const string file_shot = path_gen + "NEW_BINS/CORR_MM/ONE_TASK/SHOT_NOISE/shot_noise_nsamp_" + osSAMPM.str() 
    + "_one_task_" + argv[2] + ".dat";

  //get the bins and the volume element
  
  vector<string> path_data(4);
  path_data[0] =
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_";
  const string SPM = "PARTICLE";
  const signed long rng_seed_gen = -1381718569;
  int id=TASK;
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  //const signed long rng_seed = rng_seed_gen;
  cerr<<"READING PARTICLE DATA"<<endl;
  //load the subcube data for access to various functions related to the binning
  load_SScube OBM(sub, path_data, SPM, nsampM, rng_seed);
  vector<particle> VpartM=OBM.return_formatted_data();
  long long NNNN=OBM.info_output();
  cerr<<"number of sampled particles for subcube: "<<sub<<" is: "<<NNNN<<endl;
  //cerr<<"size of formatted particles: "<<VpartM.size()<<endl;
  //ofstream out_shot(file_shot.c_str(), ios::app);
  //out_shot<<sub<<' '<<NNNN<<endl;
  //out_shot.close();
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
  //vector<float> bin_R, bin_Z;
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
  
  vector<double> dV=OcorrMM.volume_bins();
  //vector<double> dV;

  double NormR(0.), NormMD(0.), NormMW(0.), CMD, CMW;

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) 
      ratioMM[ir][jz]=npairsMM[ir][jz]/npairsRR[ir][jz]-2.*npairsRM[ir][jz]/npairsRR[ir][jz]+1.;

  //compute the integral constraint  
  //don't use the high radius bins
  //const int NRred=NR-2;
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) {
      NormR+=(dV[ir*NZ+jz]*npairsRR[ir][jz]);
      NormMD+=(dV[ir*NZ+jz]*npairsMM[ir][jz]);
      NormMW+=(dV[ir*NZ+jz]*ratioMM[ir][jz]*npairsRR[ir][jz]); //w_vol for MM
    }
  CMD=NormMD/NormR; //1+w_vol from DD
  CMW=NormMW/NormR; //w_vol from the whole w_LS

  //compute the Landy-Szalay estimator, with and without the integral constraint

  //matter-matter

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) {
      ncorrMM[ir][jz]=ratioMM[ir][jz]; //without the IC
      ncorrMMD[ir][jz]=(npairsMM[ir][jz]/npairsRR[ir][jz]/CMD-2.*npairsRM[ir][jz]/npairsRR[ir][jz]+1.);
      ncorrMMW[ir][jz]=((1.+ratioMM[ir][jz])/(1.+CMW));
      ncorrMMW[ir][jz]-=1.;
    }

  //write the output files
  //matter-matter

  ofstream outM(file_out_MM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) 
      outM<<rb[1][ir]<<' '<<zb[1][jz]<<' '<<ncorrMM[ir][jz]<<' '<<ncorrMMD[ir][jz]<<' '<<ncorrMMW[ir][jz]<<endl;
  outM.close();
  cerr<<"\n done writing the MM output file " + file_out_MM<<endl;

  return 0;
  
}

 

