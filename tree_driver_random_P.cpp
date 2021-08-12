#include "mpi.h"
#include "RandomNumbers.h" 
#include "tree_processing_projected.h"


//assume that the master will send the same subcube to the slaves to 
//compute multiple correlation functions; will take average afterwards

using namespace std;

int main(int argc, char *argv[])

{

  if(argc != 3) {
    cerr<<"must input the number of random particles and zmax"<<endl;
    exit(10);
  }
  
  MPI::Init(argc, argv);
  int ntasks, id, master, tag;
  double wtime;
  const size_t NRAN=atoi(argv[1]);
  const float zMAX=atof(argv[2]);
  id=MPI::COMM_WORLD.Get_rank(); //get this processor's id
  ntasks=MPI::COMM_WORLD.Get_size(); //get number of processors
  wtime=MPI::Wtime();
  MPI::Status status;
  master=0;
  const double rbox=500.;
  const float rMAX=100.;//110.;//30.;                                                                                                  
  const float rMIN=0.01;
  const float zMIN=0.;
  //const float zMAX=150.;
  const int NB=30; //20; //10; 
  const int NZ=10;
  float rb[3][NB]; //radial bin vector 
  float zb[3][NZ]; //z-bin vector    
  //double npairsRR[ntasks][NB][NZ]; //array containing the pair counts
  double npairsRR[ntasks][NB][NZ];
  //double npairsRRAVE[NB][NZ], npairsRRVAR[NB][NZ];
 
  ostringstream osNB, osNZ, osNRAN, osNTASK, osZMAX;
  osNB<<NB; osNRAN<<NRAN; osNTASK<<ntasks; osNZ<<NZ; osZMAX<<zMAX;

  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/RANDOM_CATALOGUE/NEW_BINS/";
  //const string file_out_ave = path_out + "ave_pair_counts_NRAN_" + osNRAN.str() + "_binsR_" 
  //+ osNB.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str() + "_XXL.dat"; 
  const string file_gen = path_out + "pair_counts_NRAN_" + osNRAN.str() + "_binsR_" 
    + osNB.str() + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_" + osNTASK.str(); 
  string file_out[ntasks];
  for(int np=0; np<ntasks; np++) {
    ostringstream oo;
    oo<<np;
    file_out[np] = file_gen + "_task_" + oo.str() + ".dat";
  }

  if(id==0) {
    cerr<<"beginning the master process for NRAN: "<<NRAN
	<<" and ntasks : "<<ntasks<<endl; }
  //for(int j=0; j<NB; j++) 
  //  for(int iz=0; iz<NZ; iz++) {
  //npairsRRAVE[j][iz]=0.; npairsRRVAR[j][iz]=0.;
  //  }
  //}    
  //begin slave processes

  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  const int dim=3;
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed_gen = -1381718569;
  const signed long rng_seed1 = rng_seed_gen + id*id;
  const signed long rng_seed2 = rng_seed_gen/3 + id*id*id*id;
  const signed long rng_seed3 = 2*rng_seed_gen + id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '
      <<rng_seed1<<' '<<rng_seed2<<' '<<rng_seed3<<endl;
  vector<particle> VpartR;
  
  random2 OBran1(rng_seed1), OBran2(rng_seed2), OBran3(rng_seed3);
  wtime1=MPI::Wtime();
  for(size_t i=0; i<NRAN; i++) {
    vector<float> temp(3);
    temp[0]=rbox*OBran1(); temp[1]=rbox*OBran2(); temp[2]=rbox*OBran3();
    particle Otemp(temp);
    VpartR.push_back(Otemp);
  }    
  cerr<<"size of random groups: "<<VpartR.size()<<endl;
  wtime2=MPI::Wtime();
  cerr<<"time to generate random data for slave process:"<<id<<" is: "<<wtime2-wtime1<<endl;
    
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // now build the tree
    
  cerr<<"\n BUILDING TREE FOR RANDOM PARTICLES\n";
  
  wtime1=MPI::Wtime();
  build_tree OtreeR(rbox, VpartR);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the tree:= "<<wtime2-wtime1<<endl;
    
  build_tree *OR1, *OR2;
  OR1=&OtreeR; OR2=&OtreeR;
  const vector<float> vecPER(dim, 0.);
    
  //logarithmic bins
  const string dummy = "dummy";
  projected_correlations::BinType bt = projected_correlations::LOG10;
  bool flag_per = false; //true;
  cerr<<"box size: "<<rbox<<endl;
    
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  projected_correlations OcorrRR(OR1, OR2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(id==master) { 
    vector<float> bin_R=OcorrRR.return_bins_rperp();
    vector<float> bin_Z=OcorrRR.return_bins_z();
    for(int i=0; i<NB; i++) {
      rb[0][i]=bin_R[i];
      rb[1][i]=bin_R[2*NB+i];
      rb[2][i]=bin_R[NB+i];
    }
    for(int i=0; i<NZ; i++) {
      zb[0][i]=bin_Z[i]; //lower limits       
      zb[1][i]=bin_Z[2*NZ+i]; //mean     
      zb[2][i]=bin_Z[NZ+i]; //upper limit   
    }
  }
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRR=OcorrRR.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    for(int iz=0; iz<NZ; iz++) 
      npairsRR[id][j][iz]=tcorrRR[j*NZ+iz];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute RR correlations:= "<<wtime2-wtime1<<endl<<endl;

  if(id!=master) 
    MPI::COMM_WORLD.Send(&npairsRR[id], NB*NZ, MPI::DOUBLE, master, tag); 
  //now the master part
  else {
    for(int np=1; np<ntasks; np++) {
      tag=np;
      MPI::COMM_WORLD.Recv(&npairsRR[np], NB*NZ, MPI::DOUBLE, np, tag, status);
    }
    
    //compute average
    /*for(int j=0; j<NB; j++) 
      for(int iz=0; iz<NZ; iz++) {
	for(int np=0; np<ntasks; np++)
	  npairsRRAVE[j][iz]+=npairsRR[np][j][iz];
	npairsRRAVE[j][iz]/=ntasks;
	}
    //compute variance
    for(int j=0; j<NB; j++) 
      for(int iz=0; iz<NZ; iz++) {
	for(int np=0; np<ntasks; np++)
	  npairsRRVAR[j][iz]+=pow(npairsRR[np][j][iz]-npairsRRAVE[j][iz], 2);
	npairsRRVAR[j][iz]=sqrt(npairsRRVAR[j][iz]/(ntasks-1)/ntasks);
	}*/
    //write individual files
    for(int np=0; np<ntasks; np++) {
      ofstream out(file_out[np].c_str(), ios::out);
      for(int j=0; j<NB; j++)
	for(int iz=0; iz<NZ; iz++)
	  out<<rb[0][j]<<' '<<rb[1][j]<<' '<<rb[2][j]<<' '<<zb[0][iz]<<' '<<zb[1][iz]
	     <<' '<<zb[2][iz]<<' '<<npairsRR[np][j][iz]<<endl;//<<' '<<npairsRRVAR[j][iz]<<endl;
      out.close();
      cerr<<"the master write invidual file for task: "<<np<<endl;
      cerr<<file_out[np]<<endl;
    }
    /*//write average file
    cerr<<"the master writes the average file for NRAN: "<<NRAN<<" ntasks: "<<ntasks<<endl;
    ofstream outy(file_out_ave.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      for(int iz=0; iz<NZ; iz++)
	outy<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<zb[0][iz]<<' '<<zb[1][iz]
	    <<' '<<zb[2][iz]<<' '<<npairsRRAVE[ib][iz]<<' '<<npairsRRVAR[ib][iz]<<endl;
    outy.close();
    cerr<<"\n the master has written the output file " + file_out_ave<<endl;*/
  }
  MPI::COMM_WORLD.Barrier();  
  
  MPI::Finalize();
  
  return 0;
  
}

  


 
