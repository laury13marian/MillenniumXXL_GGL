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

  //adjust code to take binned groups

  const int NMASS=4;
  vector<double> MASSBINS(NMASS+1);
  MASSBINS[0]=1.e+11; MASSBINS[1]=5.e+11; MASSBINS[2]=1.e+12; MASSBINS[3]=1.e+13;
  MASSBINS[4]=1.e+18;
  vector<int> nsampB(NMASS);
  nsampB[0]=nsampG; nsampB[1]=nsampG/2; 
  nsampB[2]=1; nsampB[3]=1;
  //for(int i=1; i<NMASS; i++)
  //nsampB[i]=1.;
  bool FLAG_BINNED = true;

  if(id==master) 
    cerr<<"Run specifics: snap: "<<snap<<" subcube: "<<ns<<" nsampM: "
	<<nsampM<<" nsampG: "<<nsampG<<" NRAN: "<<NRAN<<" NTASKS: "<<ntasks<<endl<<endl;
  int NB=40; //20; //10;
  float rb[3][NB]; //bin vector
  double npairsMM[ntasks][NB], npairsGM[ntasks][NMASS][NB], npairsGG[ntasks][NMASS][NB], npairsRM[ntasks][NB],
    npairsRG[ntasks][NMASS][NB], npairsRR[ntasks][NB]; //array containing the pair counts
  double ratioMM[ntasks][NB], ratioGG[ntasks][NMASS][NB], ratioGM[ntasks][NMASS][NB];
  double ncorrMM[NB], ncorrGG[NMASS][NB], ncorrGM[NMASS][NB];
  double nvarMM[NB], nvarGG[NMASS][NB], nvarGM[NMASS][NB];

  const double rbox=500.;
  const int dim=3;

  ostringstream osNB, osSNAP, osSC, osSAMPM, osNTASK, osNRAN;
  osNB<<NB; osSNAP<<0;  osSNAP<<snap; osSC<<ns; osSAMPM<<nsampM; 
  osNTASK<<ntasks; osNRAN<<NRAN;
  const string path_out = 
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/UNPROJECTED_CORR_FUNC/BINNED/";
  string file_out_MM, file_out_GG[NMASS], file_out_GM[NMASS];
  //0=MM; 1=GG; 2=GM; 
  file_out_MM = path_out + "MM/" + "corr_func_snap_" + osSNAP.str() + 
    + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
    + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + ".dat";
  
  for(int ib=0; ib<NMASS; ib++) {
    ostringstream osM1, osM2, osSAMPG;
    osM1<<MASSBINS[ib]; osM2<<MASSBINS[ib+1]; osSAMPG<<nsampB[ib];
    
    file_out_GG[ib] = path_out + "GG/" + "corr_func_snap_" + osSNAP.str() + 
      + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + osSAMPG.str() 
      + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() + "_binM_" + osM1.str() 
      + '_' + osM2.str() + ".dat";
    
    file_out_GM[ib] = path_out + "GM/" + "corr_func_snap_" + osSNAP.str() + 
      + "_SubCube_" + osSC.str() + "_binsR_" + osNB.str() + "_nsamp_" + argv[3] 
      + '_' + osSAMPG.str() + "_NRAN_" + osNRAN.str() + "_ntasks_" + osNTASK.str() 
      + "_binM_" + osM1.str() + '_' + osM2.str()+ ".dat";
    //osM1.clear();
    //osM2.clear();    
    //cerr<<file_out_GG[ib]<<endl;
    //cerr<<file_out_GM[ib]<<endl; 
  }
  //cerr<<file_out_MM<<endl;
    
    if(id==0) {
    cerr<<"beginning the master process for snapshot: "<<snap
	<<" subcube: "<<ns<<" and MM sampling rate 1/"<<nsampM<<endl;
    for(int j=0; j<NB; j++) {
      ncorrMM[j]=0.; nvarMM[j]=0.;
      for(int ib=0; ib<NMASS; ib++) {
	ncorrGG[ib][j]=0.; ncorrGM[ib][j]=0.; nvarGG[ib][j]=0.; nvarGM[ib][j]=0.;
      }
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
  const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '<<rng_seed<<endl;
  
  cerr<<"READING PARTICLE DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBM(ns, path_data, SPM, nsampM, rng_seed);
  vector<particle> VpartM=OBM.return_formatted_data();
  cerr<<"size of formatted groups: "<<VpartM.size()<<endl;
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
  
  cerr<<"READING GROUP DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBG(ns, path_data, SPG, rng_seed, nsampB, MASSBINS, FLAG_BINNED);
  vector<particle> VpartG[NMASS];
  for(int ib=0; ib<NMASS; ib++) {
    VpartG[ib]=OBG.return_formatted_data(ib);
    cerr<<"size of formatted groups for bin "<<ib<<' '<<MASSBINS[ib]<<' '<<VpartG[ib].size()<<endl;
  }
  OBG.info_output();
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the group data for slave process: "
      <<id<<" is: "<<wtime2-wtime1<<endl;
  
 
  /*vector<double> Bound=OBG.lower_bound();
  if(id==master) {
    cerr<<"lower subcube bound: ";
    for(int i=0; i<Bound.size(); i++)
      cerr<<Bound[i]<<' ';
    cerr<<endl<<endl;
    }*/
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
  
  wtime1=MPI::Wtime();
  build_tree OtreeG[NMASS] = { build_tree(rbox, VpartG[0]), build_tree(rbox, VpartG[1]), 
			       build_tree(rbox, VpartG[2]), build_tree(rbox, VpartG[3])} ;
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the G tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OG1[NMASS], *OG2[NMASS];
  for(int ib=0; ib<NMASS; ib++) {
    OG1[ib]=&OtreeG[ib]; OG2[ib]=&OtreeG[ib]; 
  }
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
  correlations OcorrRM(OM1, OR1, bt, rMIN, rMAX, NB, flag_per, vecPER);
  correlations OcorrRR(OR1, OR2, bt, rMIN, rMAX, NB, flag_per, vecPER);

  correlations OcorrGG[NMASS] = { correlations(OG1[0], OG2[0], bt, rMIN, rMAX, NB, flag_per, vecPER), 
				  correlations(OG1[1], OG2[1], bt, rMIN, rMAX, NB, flag_per, vecPER), 
				  correlations(OG1[2], OG2[2], bt, rMIN, rMAX, NB, flag_per, vecPER), 
				  correlations(OG1[3], OG2[3], bt, rMIN, rMAX, NB, flag_per, vecPER)};

  correlations OcorrGM[NMASS] = {correlations(OM1, OG2[0], bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OM1, OG2[1], bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OM1, OG2[2], bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OM1, OG2[3], bt, rMIN, rMAX, NB, flag_per, vecPER)};

  correlations OcorrRG[NMASS] = {correlations(OG1[0], OR1, bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OG1[1], OR1, bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OG1[2], OR1, bt, rMIN, rMAX, NB, flag_per, vecPER),
				 correlations(OG1[3], OR1, bt, rMIN, rMAX, NB, flag_per, vecPER)};
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(id==master) {
    vector<float> bins=OcorrMM.bins_return();
    for(int i=0; i<NB; i++) {
      rb[0][i]=bins[i]; 
      rb[1][i]=bins[2*NB+i];
      rb[2][i]=bins[NB+i];
    }
  }
  
  const string dummy = "dummy";
  
  wtime1=MPI::Wtime();
  vector<double> tcorrMM=OcorrMM.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsMM[id][j]=tcorrMM[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute MM correlations:= "<<wtime2-wtime1<<endl<<endl;
    //const string dummy = "dummy";
  wtime1=MPI::Wtime();
  vector<double> tcorrGG[NMASS];
  for(int ib=0; ib<NMASS; ib++) {
    tcorrGG[ib]=OcorrGG[ib].dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      npairsGG[id][ib][j]=tcorrGG[ib][j];
  }
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute GG correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  wtime1=MPI::Wtime();
  vector<double> tcorrGM[NMASS];
  for(int ib=0; ib<NMASS; ib++) {    
    tcorrGM[ib]=OcorrGM[ib].dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      npairsGM[id][ib][j]=tcorrGM[ib][j];
  }
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute GM correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRM=OcorrRM.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRM[id][j]=tcorrRM[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute RM correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRG[NMASS];
  for(int ib=0; ib<NMASS; ib++) {
    tcorrRG[ib]=OcorrRG[ib].dual_tree_2p(dummy);
    for(int j=0; j<NB; j++)
      npairsRG[id][ib][j]=tcorrRG[ib][j];
  }
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute RG correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRR=OcorrRR.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRR[id][j]=tcorrRR[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
  cerr<<"time taken to compute RR correlations:= "<<wtime2-wtime1<<endl<<endl;
  
  if(id!=master) 
    {
      for(int ib=0; ib<NMASS; ib++) {
	MPI::COMM_WORLD.Send(&npairsGG[id][ib], NB, MPI::DOUBLE, master, tag);
	MPI::COMM_WORLD.Send(&npairsGM[id][ib], NB, MPI::DOUBLE, master, tag);   
	MPI::COMM_WORLD.Send(&npairsRG[id][ib], NB, MPI::DOUBLE, master, tag); 
      }
      MPI::COMM_WORLD.Send(&npairsMM[id], NB, MPI::DOUBLE, master, tag);    
      MPI::COMM_WORLD.Send(&npairsRM[id], NB, MPI::DOUBLE, master, tag);    
      MPI::COMM_WORLD.Send(&npairsRR[id], NB, MPI::DOUBLE, master, tag);   
      //MPI::COMM_WORLD.Send(&npairsGG[id], NB, MPI::DOUBLE, master, tag);    
      //MPI::COMM_WORLD.Send(&npairsGM[id], NB, MPI::DOUBLE, master, tag); 
      //MPI::COMM_WORLD.Send(&npairsRG[id], NB, MPI::DOUBLE, master, tag); 
    }
  else {      
    //now the master part
    for(int np=1; np<ntasks; np++) {
	tag=np;
	for(int ib=0; ib<NMASS; ib++) {
	  MPI::COMM_WORLD.Recv(&npairsGG[np][ib], NB, MPI::DOUBLE, np, tag, status);
	  MPI::COMM_WORLD.Recv(&npairsGM[np][ib], NB, MPI::DOUBLE, np, tag, status);
	  MPI::COMM_WORLD.Recv(&npairsRG[np][ib], NB, MPI::DOUBLE, np, tag, status);
	}

	MPI::COMM_WORLD.Recv(&npairsMM[np], NB, MPI::DOUBLE, np, tag, status);
	MPI::COMM_WORLD.Recv(&npairsRM[np], NB, MPI::DOUBLE, np, tag, status);
	MPI::COMM_WORLD.Recv(&npairsRR[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsGG[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsGM[np], NB, MPI::DOUBLE, np, tag, status);
	//MPI::COMM_WORLD.Recv(&npairsRG[np], NB, MPI::DOUBLE, np, tag, status);
      }
  
    for(int j=0; j<NB; j++) {
      for(int np=0; np<ntasks; np++) 
	ratioMM[np][j]=npairsMM[np][j]/npairsRR[np][j]-2.*npairsRM[np][j]/npairsRR[np][j];
      for(int im=0; im<NMASS; im++) 
	for(int np=0; np<ntasks; np++) { 
	  ratioGG[np][im][j]=npairsGG[np][im][j]/npairsRR[np][j]-2.*npairsRG[np][im][j]/npairsRR[np][j];
	  ratioGM[np][im][j]=npairsGM[np][im][j]/npairsRR[np][j]-npairsRM[np][j]/npairsRR[np][j]
	    -npairsRG[np][im][j]/npairsRR[np][j];
	}
      }
    for(int j=0; j<NB; j++) {
      for(int np=0; np<ntasks; np++) 
	ncorrMM[j]+=ratioMM[np][j]; 
      ncorrMM[j]/=ntasks;  //compute average
      for(int im=0; im<NMASS; im++) { 
	for(int np=0; np<ntasks; np++) {
	  ncorrGG[im][j]+=ratioGG[np][im][j];
	  ncorrGM[im][j]+=ratioGM[np][im][j];
	}
  	ncorrGG[im][j]/=ntasks; //compute average
	ncorrGM[im][j]/=ntasks;
      }
    }
    
    //compute variance
    for(int j=0; j<NB; j++) {
      for(int np=0; np<ntasks; np++) 
	nvarMM[j]+=pow(ratioMM[np][j]-ncorrMM[j], 2);
      nvarMM[j]=sqrt(nvarMM[j]/(ntasks-1)/ntasks);

      for(int im=0; im<NMASS; im++) {
	for(int np=0; np<ntasks; np++) {
	  nvarGG[im][j]+=pow(ratioGG[np][im][j]-ncorrGG[im][j], 2);
	  nvarGM[im][j]+=pow(ratioGM[np][im][j]-ncorrGM[im][j], 2);
	}
	nvarGG[im][j]=sqrt(nvarGG[im][j]/(ntasks-1)/ntasks);
	nvarGM[im][j]=sqrt(nvarGM[im][j]/(ntasks-1)/ntasks);
      }
    }
    //formal average, i.e. add the missing 1, from RR/RR
    for(int j=0; j<NB; j++) {
      ncorrMM[j]+=1.;
      for(int im=0; im<NMASS; im++) {
	ncorrGG[im][j]+=1.;
	ncorrGM[im][j]+=1.;
      }
    }
  
    cerr<<"the master writes the file for snap: "<<snap<<" SC: "<<ns<<" MM samp: "<<nsampM
        <<" and GG samp: "<<nsampG<<endl;
    
    ofstream outM(file_out_MM.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrMM[ib]<<' '<<nvarMM[ib]<<endl;
    outM.close();
    cerr<<"\n the master has written the MM output file " + file_out_MM<<endl;
  
    for(int im=0; im<NMASS; im++) {  
      ofstream outG(file_out_GG[im].c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	outG<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrGG[im][ib]<<' '<<nvarGG[im][ib]<<endl;
      outG.close();
      cerr<<"\n the master has written the GG output file " + file_out_GG[im]<<endl;
    
      ofstream outGM(file_out_GM[im].c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	outGM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrGM[im][ib]<<' '<<nvarGM[im][ib]<<endl;
      outGM.close();
      cerr<<"\n the master has written the GM output file " + file_out_GM[im]<<endl;
    }
  }
  
  MPI::COMM_WORLD.Barrier();

  FinTime=MPI::Wtime();  
  
  cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
  
  MPI::Finalize();
  
  return 0;
  
}
