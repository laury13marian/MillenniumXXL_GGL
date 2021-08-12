#include "mpi.h"
#include "read_data.h"
#include "tree_processing.h"


//assume that the master will send the same subcube to the slaves to 
//compute multiple correlation functions; will take average afterwards
//all correlations

using namespace std;

int main(int argc, char *argv[])

{

  if(argc != 6) {
    cerr<<"must input snapshot number, subcube number, MM sampling rate, GG sampling rate, NRAN."<<endl;
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
  if(id==master) 
    cerr<<"Run specifics: snap: "<<snap<<" subcube: "<<ns<<" nsampM: "
	<<nsampM<<" nsampG: "<<nsampG<<" NRAN: "<<NRAN<<" NTASKS: "<<ntasks<<endl<<endl;
  int NB=40; //20; //10;
  float rb[3][NB]; //bin vector
  double npairsMM[ntasks][NB], npairsGM[ntasks][NB], npairsGG[ntasks][NB], npairsRM[ntasks][NB],
    npairsRG[ntasks][NB], npairsRR[ntasks][NB]; //array containing the pair counts
  double ratioMM[ntasks][NB], ratioGG[ntasks][NB], ratioGM[ntasks][NB];
  double ncorrMM[NB], ncorrGG[NB], ncorrGM[NB];
  double ncorrMM10[NB], ncorrMM11[NB], ncorrMM20[NB], ncorrMM21[NB];
  double nvarMM[NB], nvarGG[NB], nvarGM[NB];
  double npairsRRAVE[NB]; //to test the average RR theory

  const double rbox=500.;
  const int dim=3;

  ostringstream osNB, osSNAP, osSC, osSAMPM, osSAMPG, osNTASK, osNRAN;
  osNB<<NB; osSNAP<<0;  osSNAP<<snap; osSC<<ns; osSAMPM<<nsampM; osSAMPG<<nsampG; 
  osNTASK<<ntasks; osNRAN<<NRAN;
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/UNPROJECTED_CORR_FUNC/TESTS/";
  string file_out_MM, file_out_GG, file_out_GM;
  //0=MM; 1=GG; 2=GM; 
  
  file_out_MM = path_out + "MM/" + "IC_test_snap_" + osSNAP.str() +   
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
    + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";
   
  /*file_out_GG = path_out + "GG/" + "cf_AVE_RR_snap_" + osSNAP.str() +  
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[4] 
    + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";
  
  file_out_GM = path_out + "GM/" + "cf_AVE_RR_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
    + '_' + argv[4] + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";*/

  /*file_out_MM = path_out + "MM/" + "corr_func_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
    + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";
   
  file_out_GG = path_out + "GG/" + "corr_func_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[4] 
    + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";
  
  file_out_GM = path_out + "GM/" + "corr_func_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
    + '_' + argv[4] + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";*/
  //cerr<<file_out_MM<<endl;
  //cerr<<file_out_GG<<endl;
  //cerr<<file_out_GM<<endl; 
   
  if(id==0) {
    cerr<<"beginning the master process for snapshot: "<<snap
	<<" subcube: "<<ns<<" and MM sampling rate 1/"<<nsampM
	<<" and GG sampling rate 1/"<<nsampG<<endl;
    for(int j=0; j<NB; j++) {
      ncorrMM[j]=0.; ncorrGG[j]=0.; ncorrGM[j]=0.;
      nvarMM[j]=0.; nvarGG[j]=0.; nvarGM[j]=0.;
      ncorrMM10[j]=0.; ncorrMM11[j]=0.; ncorrMM20[j]=0.; ncorrMM21[j]=0.;
    }
  }

  //MPI::COMM_WORLD.Bcast(&snap, 1, MPI::INT, master);
 
  //begin slave processes

  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  vector<string> path_data(4);
  path_data[0] = 
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_"; 
  
  const string SPM = "PARTICLE"; 
  const string SPG = "GROUP";
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed_gen = -1381718569; //just for subcube 0 testing, but always for RAN;
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  //cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '<<rng_seed<<endl;
  
  cerr<<"READING PARTICLE DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBM(ns, path_data, SPM, nsampM, rng_seed);
  //OBM.read_species();
  vector<particle> VpartM=OBM.return_formatted_data();
  cerr<<"size of formatted particle vector: "<<VpartM.size()<<endl;
  cerr<<endl;  
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

  /*vector<double> BINS(2);
  BINS[0]=1.e+11; BINS[1]=1.e+18;
  bool FLAG_BINNED=true;
  vector<int> NSA(1,1);*/

  /*cerr<<"READING GROUP DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBG(ns, path_data, SPG, nsampG, rng_seed);
  // load_SScube OBG(ns, path_data, SPG, rng_seed, NSA, BINS, FLAG_BINNED); 
  vector<particle> VpartG=OBG.return_formatted_data();
  cerr<<"size of formatted groups: "<<VpartG.size()<<endl;
  OBG.info_output();
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the group data for slave process: "
      <<id<<" is: "<<wtime2-wtime1<<endl;
  
  */
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
  
  /* wtime1=MPI::Wtime();
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
  
  //logarithmic bins
  
  const float rMAX=200.;//110.;//30.;
  const float rMIN=0.2;
  correlations::BinType bt = correlations::LOG10;
  bool flag_per = false; //true;
  cerr<<"box size: "<<rbox<<endl;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  correlations OcorrMM(OM1, OM2, bt, rMIN, rMAX, NB, flag_per, vecPER);
  //correlations OcorrGG(OG1, OG2, bt, rMIN, rMAX, NB, flag_per, vecPER);
  //correlations OcorrGM(OM1, OG2, bt, rMIN, rMAX, NB, flag_per, vecPER);
  correlations OcorrRM(OM1, OR1, bt, rMIN, rMAX, NB, flag_per, vecPER);
  //correlations OcorrRG(OG1, OR1, bt, rMIN, rMAX, NB, flag_per, vecPER);
  //correlations OcorrRR(OR1, OR2, bt, rMIN, rMAX, NB, flag_per, vecPER);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(id==master) {
    vector<float> bins=OcorrMM.bins_return();
    for(int i=0; i<NB; i++) {
      rb[0][i]=bins[i]; 
      rb[1][i]=bins[2*NB+i];
      rb[2][i]=bins[NB+i];
    }
  }
  //    if(id==1) { //the first slave should send the bins to the master
  //  MPI::COMM_WORLD.Send(&rb, 3*NB, MPI::FLOAT, master, tag);
  //    }
  const string dummy = "dummy";
  
  wtime1=MPI::Wtime();
  vector<double> tcorrMM=OcorrMM.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsMM[id][j]=tcorrMM[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute MM correlations:= "<<wtime2-wtime1<<endl<<endl;

  /*wtime1=MPI::Wtime();
  vector<double> tcorrGG=OcorrGG.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsGG[id][j]=tcorrGG[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute GG correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  wtime1=MPI::Wtime();
  vector<double> tcorrGM=OcorrGM.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsGM[id][j]=tcorrGM[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute GM correlations:= "<<wtime2-wtime1<<endl<<endl;*/
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRM=OcorrRM.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRM[id][j]=tcorrRM[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute RM correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  /*wtime1=MPI::Wtime();
  vector<double> tcorrRG=OcorrRG.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRG[id][j]=tcorrRG[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute RG correlations:= "<<wtime2-wtime1<<endl<<endl;*/
  
  /*wtime1=MPI::Wtime();
  vector<double> tcorrRR=OcorrRR.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRR[id][j]=tcorrRR[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute RR correlations:= "<<wtime2-wtime1<<endl<<endl;*/
    
  if(id!=master) 
    {
      MPI::COMM_WORLD.Send(&npairsMM[id], NB, MPI::DOUBLE, master, tag);    
      //MPI::COMM_WORLD.Send(&npairsGG[id], NB, MPI::DOUBLE, master, tag);    
      //MPI::COMM_WORLD.Send(&npairsGM[id], NB, MPI::DOUBLE, master, tag);    
      MPI::COMM_WORLD.Send(&npairsRM[id], NB, MPI::DOUBLE, master, tag);    
      // MPI::COMM_WORLD.Send(&npairsRG[id], NB, MPI::DOUBLE, master, tag);    
      //MPI::COMM_WORLD.Send(&npairsRR[id], NB, MPI::DOUBLE, master, tag);    
    }
  else {      
    //now the master part
    //MPI::COMM_WORLD.Recv(&rb, 3*NB, MPI::FLOAT, 1, 1, status);
    for(int np=1; np<ntasks; np++) 
      {
	tag=np;
	MPI::COMM_WORLD.Recv(&npairsMM[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsGG[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsGM[np], NB, MPI::DOUBLE, np, tag, status);
	MPI::COMM_WORLD.Recv(&npairsRM[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsRG[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsRR[np], NB, MPI::DOUBLE, np, tag, status);
      }
    //import random files
    for(int np=0; np<ntasks; np++) {
      const string path = 
	"/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/UNPROJECTED_CORR_FUNC/RANDOM_CATALOGUE/pair_counts_NRAN_" 
	+ osNRAN.str() + "_binsR_" + osNB.str() + "_NT_64"; //+ osNTASK.str();
      ostringstream oo;
      oo<<np;
      const string file = path + "_task_" + oo.str() + ".dat";
      ifstream in(file.c_str(), ios::in);
      check_stream(in, file);
      for(int j=0; j<NB; j++) {
	float b1, b2, b3; double x;
	in>>b1>>b2>>b3>>x;
	npairsRR[np][j]=x;
      }
      in.close();
    }

    //now import the average file:
  /*
    const string path = 
      "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/UNPROJECTED_CORR_FUNC/RANDOM_CATALOGUE/ave_pair_counts_NRAN_" 
      + osNRAN.str() + "_binsR_" + osNB.str() + "_NT_" + osNTASK.str() + ".dat"; 
    ifstream in(path.c_str(), ios::in);
    check_stream(in, path);
    for(int j=0; j<NB; j++) {
      float b1, b2, b3; double x, y; 
      in>>b1>>b2>>b3>>x>>y;
      npairsRRAVE[j]=x;
    }
    in.close();

    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) {
	ratioMM[np][ib]=npairsMM[np][ib]/npairsRRAVE[ib]-2.*npairsRM[np][ib]/npairsRRAVE[ib];
	ratioGG[np][ib]=npairsGG[np][ib]/npairsRRAVE[ib]-2.*npairsRG[np][ib]/npairsRRAVE[ib];
	ratioGM[np][ib]=npairsGM[np][ib]/npairsRRAVE[ib]-npairsRM[np][ib]/npairsRRAVE[ib]
	  -npairsRG[np][ib]/npairsRRAVE[ib];
      }
      }*/
    
    vector<double> dV=OcorrMM.volume_bins();
    vector<double> NormR(ntasks, 0.), NormM(ntasks, 0.), NormM2(ntasks, 0.), C1(ntasks), C2(ntasks);

    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) {
	ratioMM[np][ib]=npairsMM[np][ib]/npairsRR[np][ib]-2.*npairsRM[np][ib]/npairsRR[np][ib]+1.; //!!!
	//ratioGG[np][ib]=npairsGG[np][ib]/npairsRR[np][ib]-2.*npairsRG[np][ib]/npairsRR[np][ib];
	//ratioGM[np][ib]=npairsGM[np][ib]/npairsRR[np][ib]-npairsRM[np][ib]/npairsRR[np][ib]
	//-npairsRG[np][ib]/npairsRR[np][ib];
      }
    }
    
    for(int np=0; np<ntasks; np++) {
      for(int ib=0; ib<NB; ib++) {
	NormM[np]+=(dV[ib]*npairsMM[np][ib]);
	NormR[np]+=(dV[ib]*npairsRR[np][ib]);
	NormM2[np]+=(dV[ib]*ratioMM[np][ib]*npairsRR[np][ib]);
	}
      C1[np]=NormM[np]/NormR[np]; //1+w_vol                                                                                            
      C2[np]=NormM2[np]/NormR[np]; //1+w_vol                                                                                           
    }

    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) {
	ncorrMM[ib]+=ratioMM[np][ib];
	ncorrMM10[ib]+=(ratioMM[np][ib]*C1[np]+C1[np]-1.);
	ncorrMM11[ib]+=((ratioMM[np][ib]+1.)/C1[np]-1.);
	ncorrMM20[ib]+=(ratioMM[np][ib]*C2[np]+C2[np]-1.);
	ncorrMM21[ib]+=((ratioMM[np][ib]+1.)/C2[np]-1.);
      }
      //compute average                                                                                                              
      ncorrMM[ib]/=ntasks;                                                                                                     
      ncorrMM10[ib]/=ntasks;
      ncorrMM11[ib]/=ntasks;
      ncorrMM20[ib]/=ntasks;
      ncorrMM21[ib]/=ntasks;
    }

    /*for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) {
	ncorrMM[ib]+=ratioMM[np][ib];
	ncorrGG[ib]+=ratioGG[np][ib];
	ncorrGM[ib]+=ratioGM[np][ib];
      }
      //compute average
      ncorrMM[ib]/=ntasks;
      ncorrGG[ib]/=ntasks;
      ncorrGM[ib]/=ntasks;
    }
    //compute variance
    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) {
	nvarMM[ib]+=pow(ratioMM[np][ib]-ncorrMM[ib], 2);
	nvarGG[ib]+=pow(ratioGG[np][ib]-ncorrGG[ib], 2);
	nvarGM[ib]+=pow(ratioGM[np][ib]-ncorrGM[ib], 2);
      }
      nvarMM[ib]=sqrt(nvarMM[ib]/(ntasks-1)/ntasks);
      nvarGG[ib]=sqrt(nvarGG[ib]/(ntasks-1)/ntasks);
      nvarGM[ib]=sqrt(nvarGM[ib]/(ntasks-1)/ntasks);
      }

    //formal average, i.e. add the missing 1, from RR/RR
    for(int ib=0; ib<NB; ib++) {
      ncorrMM[ib]+=1.;
      ncorrGG[ib]+=1.;
      ncorrGM[ib]+=1.;
      }*/
    
    cerr<<"the master writes the file for snap: "<<snap<<" SC: "<<ns<<" MM samp: "<<nsampM
	<<" and GG samp: "<<nsampG<<endl;
    
    ofstream outM(file_out_MM.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outM<<rb[1][ib]<<' '<<ncorrMM[ib]<<' '<<ncorrMM10[ib]<<' '<<ncorrMM11[ib]<<' '<<ncorrMM20[ib]<<' '<<ncorrMM21[ib]<<endl;
    //outM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrMM[ib]<<' '<<nvarMM[ib]<<endl;
    outM.close();
    cerr<<"\n the master has written the MM output file " + file_out_MM<<endl;
  }    
    /*ofstream outG(file_out_GG.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outG<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrGG[ib]<<' '<<nvarGG[ib]<<endl;
    outG.close();
    cerr<<"\n the master has written the GG output file " + file_out_GG<<endl;
    
    ofstream outGM(file_out_GM.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outGM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrGM[ib]<<' '<<nvarGM[ib]<<endl;
    outGM.close();
    cerr<<"\n the master has written the GM output file " + file_out_GM<<endl;
  }
  */
  MPI::COMM_WORLD.Barrier();

  FinTime=MPI::Wtime();  
  
  cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;

  MPI::Finalize();
  
  return 0;
  
}

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  


 
