//generate tangential shear from projected correlation functions

#include "general.h"
#include "tree_processing_projected.h"
#include "/u/lmarian/ZURICH_CODES/Cosmology/utilities/Table2.h"

using namespace std;
using namespace table2;

int main(int argc, char *argv[])

{

  if(argc != 4) {
    cerr<<"must input subcube number, task number, JKsample number."<<endl;
    exit(10);
  }

  const int sub = atoi(argv[1]);
  const int TASK = atoi(argv[2]);
  const int JKsample = atoi(argv[3]);
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
  vector<float> rbL(NR), rbM(NR), rbU(NR); //radial bin vector
  float zb[3][NZ]; //z-bin vector  
  const size_t NRAN=1000000; //3000000; 
  int dim =3;
  const int Nmags=5;
  const int NM=Nmags-1; //number of magnitude/mass bins 
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; 
  const int iMAG=1;

  ostringstream osNR, osNZ, osNM, osSNAP, osSAMPM, osNTASK, osZMAX, osNRAN;
  osNR<<NR; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; osNRAN<<NRAN;
  osNTASK<<ntasks; osNZ<<NZ; osNM<<NM; osZMAX<<zMAX;
 
  //get the radial bins

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

  for(int i=0; i<NR; i++) {
    rbL[i]=bin_R[i]; //lower
    rbM[i]=bin_R[2*NR+i]; //middle
    rbU[i]=bin_R[NR+i]; //upper
  }

  for(int ib=0; ib<NR; ib++)
    cerr<<rbL[ib]<<' '<<rbM[ib]<<' '<<rbU[ib]<<endl;

  const string path_gen =
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/JACK_KNIFE/CORR_GM/";

  //now get w_mm

  /*vector<double> Wmm(NR), BarWmm(NR), Gmm(NR);
    
  const string path_in_MM = path_gen + "CORR_MM/ONE_TASK/WEIGHTED/";
  const string file_in_MM = path_in_MM + "CFw_snap_054_sub_" + argv[1] + "_binsR_" + osNR.str() + "_binsZ_" 
    + osNZ.str() + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() + + "_NT_" 
    + osNTASK.str() + "_one_task_" + argv[2] + ".dat";
  const string path_out_MM = path_gen + "CORR_MM/DELTA_SIGMA/ONE_TASK/";                    
  const string file_out_MM = path_out_MM + "delta_sigma_mm_without_rhob_sub_" + argv[1] + "_one_task_" + argv[2] + ".dat";        

  /*const string path_in_MM = path_gen + "CORR_MM/WEIGHTED/";  
  const string file_in_MM = path_in_MM + "CFw_snap_054_sub_" + argv[1] + "_binsR_" + osNR.str() + "_binsZ_" + osNZ.str() 
    + "_NRAN_" + osNRAN.str() + "_nsamp_" + osSAMPM.str() + "_zmax_" + osZMAX.str() + + "_NT_" + osNTASK.str() + ".dat";
  const string path_out_MM = path_gen + "CORR_MM/DELTA_SIGMA/";
  const string file_out_MM = path_out_MM + "delta_sigma_mm_without_rhob_sub_" + argv[1] + ".dat";*/

  /*ifstream in(file_in_MM.c_str(), ios::in);
  check_stream(in, file_in_MM);
  float a1; 
  double a2, a3, a4;
  for(int ib=0; ib<NR; ib++) {
    in>>a1>>a2>>a3>>a4;
    Wmm[ib]=a4; //want the wpp with the W integral constraint
  }
  in.close();

  //now compute wmm bar, doing a sum, instead of an integral; 
  //compute Gmm, ignoring for now the mean density and the critical density

  vector<double> X, Y;
  for(int ib=0; ib<NR; ib++) {
    double sum=0.;
    for(int j=0; j<=ib; j++) 
      sum+=(Wmm[j] * rbM[j] * (rbU[j]-rbL[j]));
    BarWmm[ib]=2./pow(rbU[ib], 2) * sum;
  }

  for(int ib=0; ib<NR; ib++) {
    X.push_back(log(rbU[ib]));
    if(BarWmm[ib]<0.) {
      cerr<<"NEGATIVE SIGMA BAR IN SHEAR MATTER-MATTER."<<endl;
      exit(1);
    }
    else 
      Y.push_back(log(BarWmm[ib]));
  }
  Table2 OI(X, Y);
  for(int ib=1; ib<NR; ib++) 
    Gmm[ib]=exp(OI(log(rbM[ib])))-Wmm[ib];
  Gmm[0]=exp(OI(log(rbM[1])))-Wmm[0];//Gmm[1]; //APPROXIMATION HERE, GIVEN HOW SMALL THE BINS ARE IN THIS REGION


  for(int ib=0; ib<NR; ib++)
    cerr<<rbM[ib]<<' '<<Gmm[ib]<<endl;
  cerr<<endl;

  ofstream out(file_out_MM.c_str(), ios::out);
  for(int ib=0; ib<NR; ib++)
    out<<rbM[ib]<<' '<<Gmm[ib]<<endl;

  out.close();
  cerr<<"done writing output file: "<<file_out_MM<<endl<<endl;
  */    
  //now do the same for the galaxies

  double Wgm[NR][NM], BarWgm[NR][NM], Ggm[NR][NM];

  const string path_in_GM = path_gen + "ONE_TASK/";
  const string file_in_GM = path_in_GM + "CFw_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_JKsample_"
    + argv[2] + "_task_" + argv[3] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_nsamp_" + osSAMPM.str() 
    + "_NT_" + osNTASK.str() + ".dat";
  
  const string path_out_GM = path_gen + "CORR_GM/DELTA_SIGMA/ONE_TASK/";
  const string file_out_GM = path_out_GM + "delta_sigma_gm_no_rhob_sub_" + argv[1] + "_JKsample_" 
    + argv[2] + "_one_task_" + argv[3] + ".dat";

  /*const string file_in_GM = path_gen + "CFw_snap_" + osSNAP.str() + "_sub_" + argv[1] + "_JKsample_"
    + argv[3] + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_nsamp_" + osSAMPM.str() 
    + "_NT_" + osNTASK.str() + ".dat";
  const string path_out_GM = path_gen + "DELTA_SIGMA/";
  const string file_out_GM = path_out_GM + "delta_sigma_gm_no_rhob_sub_" 
  + argv[1] + "_JKsample_" + argv[3] + ".dat";*/

  ifstream in2(file_in_GM.c_str(), ios::in);
  check_stream(in2, file_in_GM);
  for(int ib=0; ib<NR; ib++) {
    float b1; //radial bin
    in2>>b1;
    for(int im=0; im<NM; im++) {
      double b2, b3, b4; //mag 1, mag 2, wgm
      in2>>b2>>b3>>b4;
      Wgm[ib][im]=b4;  
    }
  }
  in2.close();
  
  for(int im=0; im<NM; im++) {
    for(int ib=0; ib<NR; ib++) {
      double sum=0.;
      for(int j=0; j<=ib; j++)
	sum+=(Wgm[j][im] * rbM[j] * (rbU[j]-rbL[j]));
      BarWgm[ib][im]=2./pow(rbU[ib], 2) * sum;
    }
  }
  
  for(int im=0; im<NM; im++) {
    vector<double> X, Y;
    for(int ib=0; ib<NR; ib++) {
      X.push_back(log(rbU[ib]));
      if(BarWgm[ib][im]<0.) {
	cerr<<"NEGATIVE SIGMA BAR IN SHEAR GALAXY-MATTER magnitude bin "<<im<<endl;
	exit(1);
      }
      else
	Y.push_back(log(BarWgm[ib][im]));
    }
    Table2 OI(X, Y);
    for(int ib=1; ib<NR; ib++)
      Ggm[ib][im]=exp(OI(log(rbM[ib])))-Wgm[ib][im];
    Ggm[0][im]=exp(OI(log(rbM[1])))-Wgm[0][im];
    //Ggm[1]; //APPROXIMATION HERE, GIVEN HOW SMALL THE BINS ARE IN THIS REGION         
    cerr<<"im is "<<im<<' '<<Ggm[0][im]<<' '<<Ggm[1][im]<<' '<<exp(OI(log(rbM[1])))<<' '<<(OI(log(rbM[2])))
	<<' '<<Wgm[1][im]<<' '<<Wgm[2][im]<<' '<<Wgm[0][im]<<endl;
  }
  cerr<<endl;

  /*for(int ib=0; ib<NR; ib++) {
    cerr<<rbM[ib]<<' ';
    for(int im=0; im<NM; im++)
      cerr<<Wgm[ib][im]<<' '<<Ggm[ib][im]<<"  ";
    cerr<<endl;
    }*/

  ofstream out2(file_out_GM.c_str(), ios::out);
  for(int ib=0; ib<NR; ib++) {
    out2<<rbM[ib]<<' ';
    for(int im=0; im<NM; im++)
      out2<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<Ggm[ib][im];
    out2<<endl;
  }
  out2.close();
  cerr<<"done writing output file: "<<file_out_GM<<endl<<endl;
  
  return 0;

}
