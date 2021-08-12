#include "mpi.h"
#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 7) {
    cerr<<"must input snapshot number, subcube number, MM sampling rate, GG sampling rate, NRAN, zmax."<<endl;
    exit(10);
  }

  MPI::Init(argc, argv);
  int ntasks, id, master, tag;
  double wtime, InitTime, FinTime;
  id=MPI::COMM_WORLD.Get_rank(); //get this processor's id 
  ntasks=MPI::COMM_WORLD.Get_size(); //get number of processors 
  wtime=MPI::Wtime();
  InitTime=MPI::Wtime();
  MPI::Status status;
  master=0;
  
  int snap=atoi(argv[1]);
  int ns=atoi(argv[2])-1; //because the script starts subcube at 1, not 0!
  int nsampM=atoi(argv[3]);
  int nsampG=atoi(argv[4]);
  const size_t NRAN=atoi(argv[5]);
  const float zMAX=atof(argv[6]);
  if(id==master) 
    cerr<<"Run specifics: snap: "<<snap<<" subcube: "<<ns<<" nsampM: "
	<<nsampM<<" nsampG: "<<nsampG<<" NRAN: "<<NRAN<<" NTASKS: "<<ntasks<<" zmax: "<<zMAX<<endl<<endl;
  const double rbox=500.;
  const float rMAX=100.;//110.;//30.;                                              
  const float rMIN=0.01;
  const float zMIN=0.;
  //const float zMAX=150.;                                                                 
  const int NB=30; //20; //10;
  const int NZ=10;
  float rb[3][NB]; //radial bin vector
  float zb[3][NZ]; //z-bin vector      
  
  double npairsMM[ntasks][NB][NZ], npairsGM[ntasks][NB][NZ], npairsGG[ntasks][NB][NZ], npairsRM[ntasks][NB][NZ],
    npairsRG[ntasks][NB][NZ], npairsRR[ntasks][NB][NZ]; //array containing the pair counts
  double ratioMM[ntasks][NB][NZ], ratioGG[ntasks][NB][NZ], ratioGM[ntasks][NB][NZ];
  double ncorrMM[NB][NZ], ncorrGG[NB][NZ], ncorrGM[NB][NZ];
  double nvarMM[NB][NZ], nvarGG[NB][NZ], nvarGM[NB][NZ];
  string file_in_RAN[ntasks], file_in_MM[ntasks], file_in_RM[ntasks];
  const int dim=3;
  
  ostringstream osNB, osNZ, osSNAP, osSC, osSAMPM, osSAMPG, osNTASK, osNRAN, osZMAX;
  osNB<<NB; osSNAP<<0;  osSNAP<<snap; osSC<<ns; osSAMPM<<nsampM; osSAMPG<<nsampG; 
  osNTASK<<ntasks; osNRAN<<NRAN; osNZ<<NZ; osZMAX<<zMAX;
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/TESTS_NEW_BINS/";
  //to register the MM and RM results
  const string file_MM = path_out + "PAIR_COUNTS_MM/pair_counts_nsamp_" + osSAMPM.str() + "_binsR_" + osNB.str()
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" ;
  const string file_RM = path_out + "PAIR_COUNTS_RM/pair_counts_nsamp_" + osSAMPM.str() + "_NRAN_" + osNRAN.str() 
    + "_binsR_" + osNB.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" ;

  const string path_RAN = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/RANDOM_CATALOGUE/NEW_BINS/";
  const string file_RAN = path_RAN + "pair_counts_NRAN_" + osNRAN.str() + "_binsR_" + osNB.str() 
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_task_" ;
  
  /*string file_out_MM, file_out_GG, file_out_GM;
    //file_out_MM = path_out + "MM/" + "proj_corr_snap_" + osSNAP.str() + 
    file_out_MM = path_out + "MM_proj_corr_snap_" + osSNAP.str() +
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_binsZ_" 
    + osNZ.str() + "_nsamp_" + osSAMPM.str() + "_NRAN_" + osNRAN.str() 
    + "_ntasks_" + osNTASK.str() + "_zmax_" + osZMAX.str() + ".dat";
  
  //file_out_GG = path_out + "GG/" + "proj_corr_snap_" + osSNAP.str() + 
  file_out_GG = path_out + "GG_proj_corr_snap_" + osSNAP.str() +
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_binsZ_" + osNZ.str() 
    + "_nsamp_" + osSAMPG.str() + "_NRAN_" + osNRAN.str() + "_ntasks_" 
    + osNTASK.str() + "_zmax_" + osZMAX.str() + ".dat";

  //file_out_GM = path_out + "GM/" + "proj_corr_snap_" + osSNAP.str() +
  file_out_GM = path_out + "GM_proj_corr_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_binsZ_" + osNZ.str() 
    + "_nsamp_" + osSAMPM.str() + '_' + osSAMPG.str() + "_NRAN_" + osNRAN.str() 
    + "_ntasks_" + osNTASK.str() + "_zmax_" + osZMAX.str() + ".dat";
  cerr<<file_out_MM<<endl;
  cerr<<file_out_GG<<endl;
  cerr<<file_out_GM<<endl; */
  
  if(id==0) {
    cerr<<"beginning the master process for snapshot: "<<snap
	<<" subcube: "<<ns<<" and MM sampling rate 1/"<<nsampM
	<<" and GG sampling rate 1/"<<nsampG<<endl;
    for(int j=0; j<NB; j++) 
      for(int jz=0; jz<NZ; jz++) {
	ncorrMM[j][jz]=0.; ncorrGG[j][jz]=0.; ncorrGM[j][jz]=0.;
	nvarMM[j][jz]=0.; nvarGG[j][jz]=0.; nvarGM[j][jz]=0.;
      }
  }
  
  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  vector<string> path_data(4);
  path_data[0] = 
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_"; 
  
  const string SPM = "PARTICLE"; 
  const string SPG = "GALAXY"; //"GROUP";
  const signed long rng_seed_gen = -1381718569;
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '<<rng_seed<<endl; 
  
  cerr<<"READING PARTICLE DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBM(ns, path_data, SPM, nsampM, rng_seed);
  vector<particle> VpartM=OBM.return_formatted_data();
  cerr<<"size of formatted particles: "<<VpartM.size()<<endl;
  OBM.info_output();
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the particle data for slave process: "
      <<id<<" is: "<<wtime2-wtime1<<endl;
  
  vector<double> Bound=OBM.lower_bound();
  if(id==master) {
    cerr<<"lower subcube bound: ";
    for(int i=0; i<Bound.size(); i++)
      cerr<<Bound[i]<<' ';
    cerr<<endl<<endl;
  }
  
  /*cerr<<"READING GALAXY DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBG(ns, path_data, SPG, nsampG, rng_seed);
  vector<particle> VpartG=OBG.return_formatted_data();
  cerr<<"size of formatted galaxies: "<<VpartG.size()<<endl;
  OBG.info_output();
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the group data for slave process: "
  <<id<<" is: "<<wtime2-wtime1<<endl;*/
    
  const signed long rng_seed1 = rng_seed_gen + id*id;
  const signed long rng_seed2 = rng_seed_gen/3 + id*id*id*id;
  const signed long rng_seed3 = 2*rng_seed_gen + id*id;
  cerr<<"seeds for random catalogue for process "<<id<<' '
      <<rng_seed1<<' '<<rng_seed2<<' '<<rng_seed3<<endl;
  vector<particle> VpartR;
  
  random2 OBran1(rng_seed1), OBran2(rng_seed2), OBran3(rng_seed3);
  wtime1=MPI::Wtime();
  for(size_t i=0; i<NRAN; i++) {
    vector<float> temp(3);
    temp[0]=Bound[0] + rbox*OBran1(); temp[1]=Bound[1] + rbox*OBran2(); temp[2]=Bound[2] + rbox*OBran3();
    particle Otemp(temp);
    VpartR.push_back(Otemp);
  }   
  cerr<<"size of random groups: "<<VpartR.size()<<endl;
  wtime2=MPI::Wtime();
  cerr<<"time to generate random data for slave process:"<<id<<" is: "<<wtime2-wtime1<<endl;

  //now build tree

  wtime1=MPI::Wtime();
  build_tree OtreeM(rbox, VpartM);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the M tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OM1, *OM2;
  OM1=&OtreeM; OM2=&OtreeM; 
  
  /*  wtime1=MPI::Wtime();
  build_tree OtreeG(rbox, VpartG);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the G tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OG1, *OG2;
  OG1=&OtreeG; OG2=&OtreeG; */
  
  wtime1=MPI::Wtime();
  build_tree OtreeR(rbox, VpartR);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the R tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OR1, *OR2;
  OR1=&OtreeR; OR2=&OtreeR;     

  const vector<float> vecPER(dim, 0.);

    projected_correlations::BinType bt = projected_correlations::LOG10;
    bool flag_per=false; //true;
    cerr<<"box size: "<<rbox<<endl;

    projected_correlations OcorrMM(OM1, OM2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);
    //projected_correlations OcorrGG(OG1, OG2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);
    //projected_correlations OcorrGM(OM1, OG2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);
    projected_correlations OcorrRM(OM1, OR2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);
    //projected_correlations OcorrRG(OG1, OR2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);
    //projected_correlations OcorrRR(OR1, OR2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    if(id==master) {
      vector<float> bin_R=OcorrMM.return_bins_rperp();
      vector<float> bin_Z=OcorrMM.return_bins_z();
      for(int i=0; i<NB; i++) {
     	rb[0][i]=bin_R[i]; 
	rb[1][i]=bin_R[2*NB+i];
	rb[2][i]=bin_R[NB+i];
      }

      for(int i=0; i<NZ; i++) {
	zb[0][i]=bin_Z[i]; 
	zb[1][i]=bin_Z[2*NZ+i];
	zb[2][i]=bin_Z[NZ+i];
      }
    }
    
    const string dummy = "dummy";
    
    wtime1=MPI::Wtime();
    vector<double> tcorrMM=OcorrMM.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
	npairsMM[id][j][jz]=tcorrMM[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
      cerr<<"time taken to compute MM correlations:= "<<wtime2-wtime1<<endl<<endl;

    wtime1=MPI::Wtime();
    vector<double> tcorrRM=OcorrRM.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
        npairsRM[id][j][jz]=tcorrRM[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master)
      cerr<<"time taken to compute RM correlations:= "<<wtime2-wtime1<<endl<<endl;

    /*    wtime1=MPI::Wtime();
    vector<double> tcorrGG=OcorrGG.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
	npairsGG[id][j][jz]=tcorrGG[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
      cerr<<"time taken to compute GG correlations:= "<<wtime2-wtime1<<endl<<endl;
   
    wtime1=MPI::Wtime();
    vector<double> tcorrGM=OcorrGM.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
	npairsGM[id][j][jz]=tcorrGM[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
    cerr<<"time taken to compute GM correlations:= "<<wtime2-wtime1<<endl<<endl;
    
    wtime1=MPI::Wtime();
    vector<double> tcorrRG=OcorrRG.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
	npairsRG[id][j][jz]=tcorrRG[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
      cerr<<"time taken to compute RG correlations:= "<<wtime2-wtime1<<endl<<endl;*/
    
    /*wtime1=MPI::Wtime();
    vector<double> tcorrRR=OcorrRR.dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      for(int jz=0; jz<NZ; jz++)
	npairsRR[id][j][jz]=tcorrRR[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
    cerr<<"time taken to compute RR correlations:= "<<wtime2-wtime1<<endl<<endl;*/
    
   
    if(id!=master) 
      {
	MPI::COMM_WORLD.Send(&npairsMM[id], NB*NZ, MPI::DOUBLE, master, tag);    
	MPI::COMM_WORLD.Send(&npairsRM[id], NB*NZ, MPI::DOUBLE, master, tag);
	//MPI::COMM_WORLD.Send(&npairsGG[id], NB*NZ, MPI::DOUBLE, master, tag);    
	//MPI::COMM_WORLD.Send(&npairsGM[id], NB*NZ, MPI::DOUBLE, master, tag);    
	//MPI::COMM_WORLD.Send(&npairsRG[id], NB*NZ, MPI::DOUBLE, master, tag);    
	//MPI::COMM_WORLD.Send(&npairsRR[id], NB*NZ, MPI::DOUBLE, master, tag);    
      }

    else {      
      //now the master part
      for(int np=1; np<ntasks; np++) 
	{
	  tag=np;
	  MPI::COMM_WORLD.Recv(&npairsMM[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  MPI::COMM_WORLD.Recv(&npairsRM[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  //MPI::COMM_WORLD.Recv(&npairsGG[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  //MPI::COMM_WORLD.Recv(&npairsGM[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  //MPI::COMM_WORLD.Recv(&npairsRG[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  //MPI::COMM_WORLD.Recv(&npairsRR[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	}
      //output MM and RM files

      for(int np=0; np<ntasks; np++) {
        ostringstream oss;
        oss<<np;
	file_in_MM[np] = file_MM + oss.str() + ".dat";
	file_in_RM[np] = file_RM + oss.str() + ".dat";
      }

      for(int np=0; np<ntasks; np++) {
	ofstream outMM(file_in_MM[np].c_str(), ios::out);
	for(int j=0; j<NB; j++)
	  for(int iz=0; iz<NZ; iz++)
	    outMM<<rb[1][j]<<' '<<zb[1][iz]<<' '<<npairsMM[np][j][iz]<<endl;
	outMM.close();
	ofstream outRM(file_in_RM[np].c_str(), ios::out);
        for(int j=0; j<NB; j++)
          for(int iz=0; iz<NZ; iz++)
            outRM<<rb[1][j]<<' '<<zb[1][iz]<<' '<<npairsRM[np][j][iz]<<endl;
        outRM.close();
      }
    }
      
      /*for(int np=0; np<ntasks; np++) {  
	ostringstream oss;
	oss<<np;
	vector<double> Vtemp;
	file_in_RAN[np] = file_RAN + oss.str() + ".dat";
	ifstream in(file_in_RAN[np].c_str(), ios::in);
	check_stream(in, file_in_RAN[np]);
	for(int i=0; i<NB*NZ; i++) {
	  double d1;// d2;
	  float r1, r2, r3, z1, z2, z3;
	  in>>r1>>r2>>r3>>z1>>z2>>z3>>d1; //>>d2;
	  Vtemp.push_back(d1);
	}
	in.close();
	for(int ir=0; ir<NB; ir++)
	  for(int iz=0; iz<NZ; iz++)
	    npairsRR[np][ir][iz]=Vtemp[ir*NZ+iz];
      }
      for(int np=0; np<ntasks; np++) 
	for(int ib=0; ib<NB; ib++) 
	  for(int jz=0; jz<NZ; jz++) { 
	    ratioMM[np][ib][jz]=npairsMM[np][ib][jz]/npairsRR[np][ib][jz]-2.*npairsRM[np][ib][jz]/npairsRR[np][ib][jz];
	    ratioGG[np][ib][jz]=npairsGG[np][ib][jz]/npairsRR[np][ib][jz]-2.*npairsRG[np][ib][jz]/npairsRR[np][ib][jz];
	    ratioGM[np][ib][jz]=npairsGM[np][ib][jz]/npairsRR[np][ib][jz]-npairsRM[np][ib][jz]/npairsRR[np][ib][jz]
	      -npairsRG[np][ib][jz]/npairsRR[np][ib][jz];
	  }
      
      for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	  for(int np=0; np<ntasks; np++) {
	    ncorrMM[ib][jz]+=ratioMM[np][ib][jz];
	    ncorrGG[ib][jz]+=ratioGG[np][ib][jz];
	    ncorrGM[ib][jz]+=ratioGM[np][ib][jz];
	  }
	  //compute average
	  ncorrMM[ib][jz]/=ntasks;
	  ncorrGG[ib][jz]/=ntasks;
	  ncorrGM[ib][jz]/=ntasks;
	}

      //compute variance
      for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	  for(int np=0; np<ntasks; np++) {
	    nvarMM[ib][jz]+=pow(ratioMM[np][ib][jz]-ncorrMM[ib][jz], 2);
	    nvarGG[ib][jz]+=pow(ratioGG[np][ib][jz]-ncorrGG[ib][jz], 2);
	    nvarGM[ib][jz]+=pow(ratioGM[np][ib][jz]-ncorrGM[ib][jz], 2);
	  }
	  nvarMM[ib][jz]=sqrt(nvarMM[ib][jz]/(ntasks-1)/ntasks);
	  nvarGG[ib][jz]=sqrt(nvarGG[ib][jz]/(ntasks-1)/ntasks);
	  nvarGM[ib][jz]=sqrt(nvarGM[ib][jz]/(ntasks-1)/ntasks);
	}

      //formal average, i.e. add the missing 1, from RR/RR
      for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	  ncorrMM[ib][jz]+=1.;
	  ncorrGG[ib][jz]+=1.;
	  ncorrGM[ib][jz]+=1.;
	}
      
      cerr<<"the master writes the file for snap: "<<snap<<" SC: "<<ns<<" MM samp: "<<nsampM
	  <<" and GG samp: "<<nsampG<<endl;
      
      ofstream outM(file_out_MM.c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++) 
	  outM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<zb[0][jz]<<' '<<zb[1][jz]
	      <<' '<<zb[2][jz]<<' '<<ncorrMM[ib][jz]<<' '<<nvarMM[ib][jz]<<endl;
      outM.close();
      cerr<<"\n the master has written the MM output file " + file_out_MM<<endl;
      
      ofstream outG(file_out_GG.c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++) 
	  outG<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<zb[0][jz]<<' '<<zb[1][jz]
	      <<' '<<zb[2][jz]<<' '<<ncorrGG[ib][jz]<<' '<<nvarGG[ib][jz]<<endl;
      outG.close();
      cerr<<"\n the master has written the GG output file " + file_out_GG<<endl;
      
      ofstream outGM(file_out_GM.c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++) 
	  outGM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<zb[0][jz]<<' '<<zb[1][jz]
	      <<' '<<zb[2][jz]<<' '<<ncorrGM[ib][jz]<<' '<<nvarGM[ib][jz]<<endl;
      outGM.close();
      cerr<<"\n the master has written the GM output file " + file_out_GM<<endl;
      }*/
    
    MPI::COMM_WORLD.Barrier();
    FinTime=MPI::Wtime();
    cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
    MPI::Finalize();
  
    return 0;
    
}

 
