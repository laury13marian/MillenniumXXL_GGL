//driver to process the subcube MM, GG, GM data

#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 2) {
    cerr<<"must input subcube number"<<endl;
    exit(10);
  }

  const int sub = atoi(argv[1]);
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

  double ratioMM[ntasks][NR][NZ], ncorrMM[NR][NZ], ncorrMMD[NR][NZ], ncorrMMW[NR][NZ];
  double npairsMM[ntasks][NR][NZ], npairsRM[ntasks][NR][NZ], npairsRR[ntasks][NR][NZ];

  //for galaxy-galaxy

  double ratioGG[ntasks][NM][NR][NZ], ncorrGG[NM][NR][NZ], ncorrGGD[NM][NR][NZ], ncorrGGW[NM][NR][NZ];
  double npairsGG[ntasks][NM][NR][NZ], npairsRG[ntasks][NM][NR][NZ];

  // galaxy-matter

  double ratioGM[ntasks][NM][NR][NZ], ncorrGM[NM][NR][NZ], ncorrGMD[NM][NR][NZ], ncorrGMW[NM][NR][NZ];
  double npairsGM[ntasks][NM][NR][NZ];

  string file_in_RR[ntasks], file_in_MM_RM[ntasks], file_in_GG[ntasks], file_in_RG[ntasks], file_in_GM[ntasks];
  //string file_in_RM[ntasks], file_in_MM[ntasks];
   
  const int dim=3;
  
  ostringstream osNR, osNZ, osNM, osSNAP, osSAMPM, osNTASK, osZMAX, osNRAN;
  osNR<<NR; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; osNRAN<<NRAN;
  osNTASK<<ntasks; osNZ<<NZ; osNM<<NM; osZMAX<<zMAX;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/"; 
  
  const string file_RR = path_gen + "RANDOM_CATALOGUE/NEW_BINS/pair_counts_NRAN_" + osNRAN.str() + "_binsR_" + osNR.str()
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" ; 

  const string file_MM = path_gen + "NEW_BINS/PAIR_COUNTS_MM_RM/PC_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_task_";
  
  string file_GG;
  if(sub<7)
    file_GG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() 
      + "_NT_64_task_";
  else
    file_GG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
      + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() 
      + "_NT_" + osNTASK.str() + "_task_";

  const string file_RG = path_gen + "NEW_BINS/PAIR_COUNTS_GG_RG/RG_snap_" + osSNAP.str() + "_sub_" + argv[1] 
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_";

  const string file_GM = path_gen + "NEW_BINS/PAIR_COUNTS_GM/PC_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsM_" 
    + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + "_task_";

  /*const string file_MM = path_gen + "TESTS_NEW_BINS/PAIR_COUNTS_MM/pair_counts_nsamp_" + osSAMPM.str() + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_";
  const string file_RM = path_gen + "TESTS_NEW_BINS/PAIR_COUNTS_RM/pair_counts_nsamp_" + osSAMPM.str() + "_NRAN_" + osNRAN.str()
  + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_";*/

  for(int np=0; np<ntasks; np++) {
    ostringstream oss;
    oss<<np;
    file_in_RR[np] = file_RR + oss.str() + ".dat";
    file_in_MM_RM[np] = file_MM + oss.str() + ".dat";
    file_in_GG[np] = file_GG + oss.str() + ".dat";
    file_in_RG[np] = file_RG + oss.str() + ".dat";
    file_in_GM[np] = file_GM + oss.str() + ".dat";
    //file_in_RM[np] = file_RM + oss.str() + ".dat";
    //file_in_MM[np] = file_MM + oss.str() + ".dat";
  }


  const string file_out_MM = path_gen + "NEW_BINS/CORR_MM/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsR_" 
    + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() 
    + "_NT_" + osNTASK.str() + ".dat";

  const string file_out_GG = path_gen + "NEW_BINS/CORR_GG/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] 
    + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + ".dat";

  const string file_out_GM = path_gen + "NEW_BINS/CORR_GM/CF_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_binsM_" 
    + osNM.str() + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() 
    + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + ".dat";
  
  /*const string file_out_MM = path_gen + "TESTS_NEW_BINS/CORR_FUNC/MM_proj_corr_snap_" + osSNAP.str() + "_NRAN_" + osNRAN.str() + 
    "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() + "_nsamp_" + osSAMPM.str() + "_ntasks_" + osNTASK.str() 
    + "_zmax_" + osZMAX.str() + ".dat";
    cerr<<file_out_MM<<endl;*/
  
  for(int j=0; j<NR; j++) 
    for(int jz=0; jz<NZ; jz++) {
      ncorrMM[j][jz]=0.; ncorrMMD[j][jz]=0.; ncorrMMW[j][jz]=0.; 
    }
  for(int im=0; im<NM; im++)
    for(int j=0; j<NR; j++) 
      for(int jz=0; jz<NZ; jz++) {
	ncorrGG[im][j][jz]=0.; ncorrGGD[im][j][jz]=0.; ncorrGGW[im][j][jz]=0.; 
	ncorrGM[im][j][jz]=0.; ncorrGMD[im][j][jz]=0.; ncorrGMW[im][j][jz]=0.; 
      }

  //get the bins and the volume element
  
  vector<string> path_data(4);
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
  cerr<<"size of formatted particles: "<<VpartM.size()<<endl;

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
  
  for(int np=0; np<ntasks; np++) {
    
    ifstream inRR(file_in_RR[np].c_str(), ios::in);
    check_stream(inRR, file_in_RR[np]);
    float x1, x2, x3, x4, x5, x6;
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++)
	inRR>>x1>>x2>>x3>>x4>>x5>>x6>>npairsRR[np][ir][iz];
    inRR.close();

    ifstream inMM(file_in_MM_RM[np].c_str(), ios::in);
    check_stream(inMM, file_in_MM_RM[np]);
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++)
        inMM>>x1>>x2>>npairsMM[np][ir][iz]>>npairsRM[np][ir][iz];
    inMM.close();

    ifstream inGG(file_in_GG[np].c_str(), ios::in);
    check_stream(inGG, file_in_GG[np]);
    double d1, d2;
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++) {
	inGG>>x1>>x2;
	for(int im=0; im<NM; im++)
	  inGG>>d1>>d2>>npairsGG[np][im][ir][iz];
      }
    inGG.close();

    ifstream inRG(file_in_RG[np].c_str(), ios::in);
    check_stream(inRG, file_in_RG[np]);
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++) {
	inRG>>x1>>x2;
	for(int im=0; im<NM; im++)
	  inRG>>d1>>d2>>npairsRG[np][im][ir][iz];
      }
    inRG.close();

    ifstream inGM(file_in_GM[np].c_str(), ios::in);
    check_stream(inGM, file_in_GM[np]);
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++) {
	inGM>>x1>>x2;
	for(int im=0; im<NM; im++)
	  inGM>>d1>>d2>>npairsGM[np][im][ir][iz];
      }
    inGM.close();

    /*ifstream inRM(file_in_RM[np].c_str(), ios::in);
    check_stream(inRM, file_in_RM[np]);
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++)
        inRM>>x1>>x2>>npairsRM[np][ir][iz];
    inRM.close();

    ifstream inMM(file_in_MM[np].c_str(), ios::in);
    check_stream(inMM, file_in_MM[np]);
    for(int ir=0; ir<NR; ir++)
      for(int iz=0; iz<NZ; iz++)
        inMM>>x1>>x2>>npairsMM[np][ir][iz];
	inMM.close();*/
  }
  
  vector<double> dV=OcorrMM.volume_bins();
  //vector<double> dV;

  vector<double> NormR(ntasks, 0.), NormMD(ntasks, 0.), NormMW(ntasks, 0.), CMD(ntasks), CMW(ntasks);
  vector<double> NormGD(ntasks*NM, 0.), NormGW(ntasks*NM, 0.), CGD(ntasks*NM), CGW(ntasks*NM);
  vector<double> NormGMD(ntasks*NM, 0.), NormGMW(ntasks*NM, 0.), CGMD(ntasks*NM), CGMW(ntasks*NM);

  for(int np=0; np<ntasks; np++) {
    for(int ir=0; ir<NR; ir++) 
      for(int jz=0; jz<NZ; jz++) 
	ratioMM[np][ir][jz]=npairsMM[np][ir][jz]/npairsRR[np][ir][jz]-2.*npairsRM[np][ir][jz]/npairsRR[np][ir][jz]+1.;

    for(int im=0; im<NM; im++) 
      for(int ir=0; ir<NR; ir++) 
	for(int jz=0; jz<NZ; jz++) {
	  ratioGG[np][im][ir][jz]=npairsGG[np][im][ir][jz]/npairsRR[np][ir][jz]-2.*npairsRG[np][im][ir][jz]/npairsRR[np][ir][jz]+1.;
	  ratioGM[np][im][ir][jz]=npairsGM[np][im][ir][jz]/npairsRR[np][ir][jz]-npairsRG[np][im][ir][jz]/npairsRR[np][ir][jz]
	    -npairsRM[np][ir][jz]/npairsRR[np][ir][jz] + 1.;
	}
  }

  //compute the integral constraint  

  for(int np=0; np<ntasks; np++) {
    for(int ir=0; ir<NR; ir++)
      for(int jz=0; jz<NZ; jz++) {
	NormR[np]+=(dV[ir*NZ+jz]*npairsRR[np][ir][jz]);
	NormMD[np]+=(dV[ir*NZ+jz]*npairsMM[np][ir][jz]);
	NormMW[np]+=(dV[ir*NZ+jz]*ratioMM[np][ir][jz]*npairsRR[np][ir][jz]); //w_vol for MM
      }
    CMD[np]=NormMD[np]/NormR[np]; //1+w_vol from DD
    CMW[np]=NormMW[np]/NormR[np]; //w_vol from the whole w_LS
  }


  for(int np=0; np<ntasks; np++) 
    for(int im=0; im<NM; im++) {
      for(int ir=0; ir<NR; ir++)
	for(int jz=0; jz<NZ; jz++) {
	  NormGD[np*NM+im]+=(dV[ir*NZ+jz]*npairsGG[np][im][ir][jz]);
	  NormGW[np*NM+im]+=(dV[ir*NZ+jz]*ratioGG[np][im][ir][jz]*npairsRR[np][ir][jz]); //w_vol for gal-gal
	  NormGMD[np*NM+im]+=(dV[ir*NZ+jz]*npairsGM[np][im][ir][jz]);
	  NormGMW[np*NM+im]+=(dV[ir*NZ+jz]*ratioGM[np][im][ir][jz]*npairsRR[np][ir][jz]); //w_vol for gal-matter
	}
      CGD[np*NM+im]=NormGD[np*NM+im]/NormR[np]; //1+w_vol from DD for gal-gal
      CGW[np*NM+im]=NormGW[np*NM+im]/NormR[np]; //w_vol from the whole w_LS for gal-gal
      CGMD[np*NM+im]=NormGMD[np*NM+im]/NormR[np]; //1+w_vol from DD for gal-matter
      CGMW[np*NM+im]=NormGMW[np*NM+im]/NormR[np]; //w_vol from the whole w_LS for gal-matter
    }

  //compute the Landy-Szalay estimator, with and without the integral constraint

  //matter-matter

  for(int ir=0; ir<NR; ir++) 
    for(int jz=0; jz<NZ; jz++) {
      for(int np=0; np<ntasks; np++) {
	ncorrMM[ir][jz]+=ratioMM[np][ir][jz]; //without the IC
	ncorrMMD[ir][jz]+=(npairsMM[np][ir][jz]/npairsRR[np][ir][jz]/CMD[np]
			   -2.*npairsRM[np][ir][jz]/npairsRR[np][ir][jz]+1.);
	ncorrMMW[ir][jz]+=((1.+ratioMM[np][ir][jz])/(1.+CMW[np]));
      }
      //compute average
      ncorrMM[ir][jz]/=ntasks;
      ncorrMMD[ir][jz]/=ntasks;
      ncorrMMW[ir][jz]/=ntasks;
      ncorrMMW[ir][jz]-=1.;
    }
      
  //galaxy-galaxy

  for(int im=0; im<NM; im++)
    for(int ir=0; ir<NR; ir++) 
      for(int jz=0; jz<NZ; jz++) {
	for(int np=0; np<ntasks; np++) {
	  ncorrGG[im][ir][jz]+=ratioGG[np][im][ir][jz]; //no IC
	  ncorrGM[im][ir][jz]+=ratioGM[np][im][ir][jz]; //no IC
	  ncorrGGD[im][ir][jz]+=(npairsGG[np][im][ir][jz]/npairsRR[np][ir][jz]/CGD[np*NM+im]
				 -2.*npairsRG[np][im][ir][jz]/npairsRR[np][ir][jz]+1.);
	  ncorrGMD[im][ir][jz]+=(npairsGM[np][im][ir][jz]/npairsRR[np][ir][jz]/CGMD[np*NM+im]-npairsRG[np][im][ir][jz]/npairsRR[np][ir][jz]
				 -npairsRM[np][ir][jz]/npairsRR[np][ir][jz]+1.);
	  ncorrGGW[im][ir][jz]+=((1.+ratioGG[np][im][ir][jz])/(1.+CGW[np*NM+im]));
	  ncorrGMW[im][ir][jz]+=((1.+ratioGM[np][im][ir][jz])/(1.+CGMW[np*NM+im]));

      }
	//compute average
	ncorrGG[im][ir][jz]/=ntasks;
	ncorrGGD[im][ir][jz]/=ntasks;
	ncorrGGW[im][ir][jz]/=ntasks;
	ncorrGGW[im][ir][jz]-=1.;
	ncorrGM[im][ir][jz]/=ntasks;
	ncorrGMD[im][ir][jz]/=ntasks;
	ncorrGMW[im][ir][jz]/=ntasks;
	ncorrGMW[im][ir][jz]-=1.;
      }    

  //write the output files
  //matter-matter

  ofstream outM(file_out_MM.c_str(), ios::out);
  for(int ir=0; ir<NR; ir++)
    for(int jz=0; jz<NZ; jz++) 
      outM<<rb[1][ir]<<' '<<zb[1][jz]<<' '<<ncorrMM[ir][jz]<<' '<<ncorrMMD[ir][jz]<<' '<<ncorrMMW[ir][jz]<<endl;
  outM.close();
  cerr<<"\n done writing the MM output file " + file_out_MM<<endl;

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

 

