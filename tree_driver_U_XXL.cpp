//driver to read the whole XXL data, not subcubes
#include "mpi.h"
#include "read_data.h"
#include "tree_processing.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 3) {
    cerr<<"must input snapshot number, MM sampling rate"<<endl;
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
  int nsampM=atoi(argv[2]);
  if(id==master) 
    cerr<<"Run specifics: snap: "<<snap<<" nsampM: "<<nsampM<<" NTASKS: "<<ntasks<<endl<<endl;
  const double rbox=3000.;
  const float rMAX=200.;//110.;//30.;                                                                                        
  const float rMIN=0.2;
  const int NB=40; 
  float rb[3][NB]; //radial bin vector
  double npairsMM[ntasks][NB], ratioMM[ntasks][NB], ncorrMM[NB], nvarMM[NB];
  const int dim=3;
  
  ostringstream osNB, osSNAP, osSAMPM, osNTASK;
  osNB<<NB; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; 
  osNTASK<<ntasks; 
  const string path_out = 
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/UNPROJECTED_CORR_FUNC/WHOLE_XXL/";
  string file_out_MM;
  file_out_MM = path_out + "MM_corr_snap_" + osSNAP.str() + "_binsR_" + osNB.str() 
    + "_nsamp_" + osSAMPM.str() + "_ntasks_" + osNTASK.str() + ".XXL.PB.dat";
  cerr<<file_out_MM<<endl;
  
  if(id==0) {
    cerr<<"beginning the master process for snapshot: "<<snap
	<<" and MM sampling rate 1/"<<nsampM<<endl;
    for(int j=0; j<NB; j++) {
      ncorrMM[j]=0.; nvarMM[j]=0.; 
    }
  }
  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  cerr<<"READING PARTICLE DATA"<<endl;
  vector<particle> VpartM;
  ostringstream osID;
  osID<<id;
  const string path_in = "/u/lmarian/LauraTreeCode/LauraMPI/XXL_DATA/COMBINED/";
  const string file_in = path_in + "particles_XXL_snap_" + osSNAP.str() + "_nsamp_" 
    + osSAMPM.str() + "_ntasks_" + osNTASK.str() + "_task_" + osID.str() + ".dat";
  ifstream in(file_in.c_str(), ios::in);
  check_stream(in, file_in);
  while(!in.eof()) {
    vector<float> temp(3);
    in>>temp[0]>>temp[1]>>temp[2];
    if(temp[0]==0. && temp[1]==0. && temp[2]==0.)
      continue;
    particle OT(temp);
    VpartM.push_back(OT);
  }
  in.close();
  cerr<<"done reading particle vector for task "<<id<<endl;
  cerr<<"size of imported particle vector: "<<VpartM.size()<<endl;
  MPI::COMM_WORLD.Barrier();
  
  //now build tree
  
  wtime1=MPI::Wtime();
  build_tree OtreeM(rbox, VpartM);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the M tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OM1, *OM2;
  OM1=&OtreeM; OM2=&OtreeM; 
  
  const vector<float> vecPER(dim, 0.);
  correlations::BinType bt = correlations::LOG10;
  bool flag_per=true; //false; //true;
  cerr<<"box size: "<<rbox<<endl;
  
  correlations OcorrMM(OM1, OM2, bt, rMIN, rMAX, NB, flag_per, vecPER);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(id==master) {
    vector<float> bin_R=OcorrMM.bins_return();
    for(int i=0; i<NB; i++) {
      rb[0][i]=bin_R[i]; 
      rb[1][i]=bin_R[2*NB+i];
      rb[2][i]=bin_R[NB+i];
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
  
  if(id!=master) 
    MPI::COMM_WORLD.Send(&npairsMM[id], NB, MPI::DOUBLE, master, tag);    
  else {      
    //now the master part
    for(int np=1; np<ntasks; np++) 
      {
	tag=np;
	MPI::COMM_WORLD.Recv(&npairsMM[np], NB, MPI::DOUBLE, np, tag, status);
      }
    vector<double> npairsRR=OcorrMM.generate_RR(false, dummy);
    
    for(int np=0; np<ntasks; np++) 
      for(int ib=0; ib<NB; ib++) 
	ratioMM[np][ib]=npairsMM[np][ib]/npairsRR[ib];
	  
    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) 
	ncorrMM[ib]+=ratioMM[np][ib];
	//compute average
      ncorrMM[ib]/=ntasks;
    }

    //compute variance
    for(int ib=0; ib<NB; ib++) {
      for(int np=0; np<ntasks; np++) 
	nvarMM[ib]+=pow(ratioMM[np][ib]-ncorrMM[ib], 2);
      nvarMM[ib]=sqrt(nvarMM[ib]/(ntasks-1)/ntasks);
    }

    //formal average, i.e. subtract the missing 1, from -2DR/RR+RR/RR, and DR=RR
    for(int ib=0; ib<NB; ib++) 
      ncorrMM[ib]-=1.;
    
    cerr<<"the master writes the file for snap: "<<snap<<" MM samp: "<<nsampM<<endl;      
    ofstream outM(file_out_MM.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outM<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<ncorrMM[ib]<<' '<<nvarMM[ib]<<endl;
    outM.close();
    cerr<<"\n the master has written the MM output file " + file_out_MM<<endl;
    
  }
  
  MPI::COMM_WORLD.Barrier();
  FinTime=MPI::Wtime();
  cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
  MPI::Finalize();
  
  return 0;
  
}

