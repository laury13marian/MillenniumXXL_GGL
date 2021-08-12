#include "mpi.h"
#include "RandomNumbers.h" 
#include "tree_processing.h"


//assume that the master will send the same subcube to the slaves to 
//compute multiple correlation functions; will take average afterwards

using namespace std;

int main(int argc, char *argv[])

{

  if(argc != 2) {
    cerr<<"must input the number of random particles."<<endl;
    exit(10);
  }
  
  MPI::Init(argc, argv);
  int ntasks, id, master, tag;
  double wtime;
  const size_t NRAN=atoi(argv[1]);
  id=MPI::COMM_WORLD.Get_rank(); //get this processor's id
  ntasks=MPI::COMM_WORLD.Get_size(); //get number of processors
  wtime=MPI::Wtime();
  MPI::Status status;
  master=0;
  const double rbox=500.;
  int NB=40; //20; //10;
  float rb[3][NB]; //bin vector
  double npairsRR[ntasks][NB]; //array containing the pair counts
  double npairsAVE[NB], npairsVAR[NB];
 
  ostringstream osNB, osNRAN, osNTASK;
  osNB<<NB; osNRAN<<NRAN; osNTASK<<ntasks;

  const string path_out = 
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/UNPROJECTED_CORR_FUNC/RANDOM_CATALOGUE/";
  const string file_out_ave = path_out + "ave_pair_counts_NRAN_" + osNRAN.str() + "_binsR_" 
    + osNB.str() + "_NT_" + osNTASK.str() + ".dat"; 
  const string file_gen = path_out + "pair_counts_NRAN_" + osNRAN.str() + "_binsR_" 
    + osNB.str() + "_NT_" + osNTASK.str(); 
  string file_out[ntasks];
  for(int np=0; np<ntasks; np++) {
    ostringstream oo;
    oo<<np;
    file_out[np] = file_gen + "_task_" + oo.str() + ".dat";
  }

  if(id==0) {
    cerr<<"beginning the master process for NRAN: "<<NRAN
	<<" and ntasks : "<<ntasks<<endl;
    for(int j=0; j<NB; j++) {
      npairsAVE[j]=0.; npairsVAR[j]=0.;
    }
  }    
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
    
  const float rMAX=200.;//110.;//30.;
  const float rMIN=0.2;
  const string dummy = "dummy";
  correlations::BinType bt = correlations::LOG10;
  bool flag_per = false; //true;
  cerr<<"box size: "<<rbox<<endl;
    
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  correlations OcorrRR(OR1, OR2, bt, rMIN, rMAX, NB, flag_per, vecPER);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(id==master) { 
    vector<float> bins=OcorrRR.bins_return();
    for(int i=0; i<NB; i++) {
      rb[0][i]=bins[i]; 
      rb[1][i]=bins[2*NB+i];
      rb[2][i]=bins[NB+i];
    }
  }
  
  wtime1=MPI::Wtime();
  vector<double> tcorrRR=OcorrRR.dual_tree_2p(dummy);
  for(int j=0; j<NB; j++)
    npairsRR[id][j]=tcorrRR[j];
  wtime2=MPI::Wtime();
  MPI::COMM_WORLD.Barrier();
  if(id==master) 
    cerr<<"time taken to compute RR correlations:= "<<wtime2-wtime1<<endl<<endl;
 

  if(id!=master) 
    MPI::COMM_WORLD.Send(&npairsRR[id], NB, MPI::DOUBLE, master, tag); 
  //now the master part
  else {
    for(int np=1; np<ntasks; np++) {
      tag=np;
      MPI::COMM_WORLD.Recv(&npairsRR[np], NB, MPI::DOUBLE, np, tag, status);
    }
    
    //write the individual files
    for(int np=0; np<ntasks; np++) {
      ofstream out(file_out[np].c_str(), ios::out);
      for(int j=0; j<NB; j++)
	out<<rb[0][j]<<' '<<rb[1][j]<<' '<<rb[2][j]<<' '<<npairsRR[np][j]<<endl;
      out.close();
      cerr<<"the master write invidual file for task: "<<np<<endl;
      cerr<<file_out[np]<<endl;
    }
    //compute average
    for(int j=0; j<NB; j++) {
      for(int np=0; np<ntasks; np++)
	npairsAVE[j]+=npairsRR[np][j];
      npairsAVE[j]/=ntasks;
    }
    //compute variance
    for(int j=0; j<NB; j++) {
      for(int np=0; np<ntasks; np++) 
	npairsVAR[j]+=pow(npairsRR[np][j]-npairsAVE[j], 2);
      npairsVAR[j]=sqrt(npairsVAR[j]/(ntasks-1)/ntasks);
    }
    cerr<<"the master writes the average file for NRAN: "<<NRAN<<" ntasks: "<<ntasks<<endl;
    ofstream outy(file_out_ave.c_str(), ios::out);
    for(int ib=0; ib<NB; ib++)
      outy<<rb[0][ib]<<' '<<rb[1][ib]<<' '<<rb[2][ib]<<' '<<npairsAVE[ib]<<' '<<npairsVAR[ib]<<endl;
    outy.close();
    cerr<<"\n the master has written the output file " + file_out_ave<<endl;
  }
  MPI::COMM_WORLD.Barrier();  

  MPI::Finalize();
  
  return 0;
  
}

  
  //generate random-random pair counts, using PBC Nran*(Nran-1)/2, where Nran is the expected random number per radial shell, i.e. Vshell * nran = Vshell * Ntot/Vtot, with Ntot matching the number of galaxies

  /*const long unsigned Ntot=Vpart.size();
  ostringstream osTOT;
  osTOT<<Ntot;
  cerr<<"number of used galaxies: "<<Ntot<<endl;
  const double Vtot=pow(rbox, 3);
  const double PI=4.*atan(1.);
  vector<double> rr;
  for(int i=0; i<NB; i++) {
    rr.push_back(4.*PI/3. * (pow(rb[NB+i], 3)-pow(rb[i], 3)));
    cerr<<rb[i]<<' '<<rr[i]<<endl; //shell volume
  }
  cerr<<endl<<endl;
  
  //GENERATE DATA-RANDOM AND RANDOM-RANDOM COUNTS

  const string path_out2 = "/afs/mpa/home/lmarian/TREE_CODE/DATA/RANDOM_RESULTS/";
  const string file_dr = path_out + "analytical_data_random_pair_counts_NRAN_" 
    + osTOT.str() + "_croton_mini_bins_" + osNB.str() + ".dat";
  const string file_rr = path_out2 + "analytical_RR_pair_counts_NRAN_" + osTOT.str() 
    + "_bins_" + osNB.str() + ".dat";
  ofstream out1(file_dr.c_str(), ios::out);
  ofstream out2(file_rr.c_str(), ios::out);
  for(int i=0; i<NB; i++) {
    out1<<rb[i]<<' '<<rb[2*NB+i]<<' '<<rb[NB+i]<<' '<<rr[i]/Vtot<<endl;
    out2<<rb[i]<<' '<<rb[2*NB+i]<<' '<<rb[NB+i]<<' '<<rr[i]/Vtot<<endl;
  }
  out1.close();
  out2.close();
  cerr<<"done writing files " + file_dr<<endl;
  cerr<<" and " + file_rr<<endl;*/

  


 
