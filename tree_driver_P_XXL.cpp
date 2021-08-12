//driver to read the whole XXL data, not subcubes
#include "mpi.h"
#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 4) {
    cerr<<"must input snapshot number, MM sampling rate, zmax."<<endl;
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
  const float zMAX=atof(argv[3]);
  const double rbox=500.;//3000.;
  const float rMAX=150.;//110.;//30.;                                              
  const float rMIN=0.2;
  const float zMIN=0.;
  const int NB=30; //20; //10;
  const int NZ=10;
  float rb[3][NB]; //radial bin vector
  float zb[3][NZ]; //z-bin vector 
  const size_t NRAN=2000000;     
  const signed long rng_seed_gen = -1381718569;
  float SIDE=0.; //chosen subcube for replication
  if(id==master)
    cerr<<"Run specifics: snap: "<<snap<<" nsampM: "<<nsampM<<" NTASKS: "<<ntasks<<" zmax: "<<zMAX
	<<" box: "<<rbox<<" NRAN: "<<NRAN<<endl<<endl;
  double npairsMM[ntasks][NB][NZ], ratioMM[ntasks][NB][NZ], ncorrMM[NB][NZ], nvarMM[NB][NZ]; 
  double ncorrMM10[NB][NZ], ncorrMM11[NB][NZ], ncorrMM20[NB][NZ], ncorrMM21[NB][NZ];
  double npairsRM[ntasks][NB][NZ], npairsRR[ntasks][NB][NZ], DDave[NB][NZ], DRave[NB][NZ], RRave[NB][NZ];
  string file_in_RAN[ntasks]; 
  const int dim=3;
  
  ostringstream osNB, osNZ, osSNAP, osSAMPM, osNTASK, osZMAX, osNRAN;
  osNB<<NB; osSNAP<<0;  osSNAP<<snap; osSAMPM<<nsampM; osNRAN<<NRAN;
  osNTASK<<ntasks; osNZ<<NZ; osZMAX<<zMAX;

  const string path_RAN = 
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/PROJECTED_CORR_FUNC/RANDOM_CATALOGUE/"; //XXL_PB/ZMAX_" + osZMAX.str() +'/';
  const string file_RAN = path_RAN + "pair_counts_NRAN_" + osNRAN.str() + "_binsR_" + osNB.str() 
    + "_binsZ_" + osNZ.str() + "_zmax_" + osZMAX.str() + "_NT_64" + "_task_" ; //MODIFIED HERE TO FIXED NTASKS!!!!!!
  //  + osNTASK.str() +
  /*const string path_out = "/u/lmarian/LauraTreeCode/LauraMPI/XXL_DATA/V2/";
  string file_out[ntasks];
  for(int i=0; i<ntasks; i++) {
    ostringstream oss;
    oss<<i;
    file_out[i] = path_out + "particles_XXL_snap_" + osSNAP.str() + "_nsamp_" + osSAMPM.str() 
    + "_ntasks_" + osNTASK.str() + "_task_" + oss.str() + "_slices_1500_1600.dat";
    }*/
  const string path_out = 
    "/u/lmarian/LauraTreeCode/LauraMPI/RESULTS/PROJECTED_CORR_FUNC/TESTS/";
  string file_out_MM;
  file_out_MM = path_out + "MM_proj_corr_snap_" + osSNAP.str() + "_binsR_" + osNB.str() + "_binsZ_" 
    + osNZ.str() + "_nsamp_" + osSAMPM.str() + "_ntasks_" + osNTASK.str() 
    + "_zmax_" + osZMAX.str() + "_sub_0_v2.dat";
  cerr<<file_out_MM<<endl;
  
  if(id==0) {
    cerr<<"beginning the master process for snapshot: "<<snap
	<<" and MM sampling rate 1/"<<nsampM<<endl;
    for(int j=0; j<NB; j++) 
      for(int jz=0; jz<NZ; jz++) {
	ncorrMM[j][jz]=0.; nvarMM[j][jz]=0.; 
	DDave[j][jz]=0.; DRave[j][jz]=0.; RRave[j][jz]=0.;
	ncorrMM10[j][jz]=0.; ncorrMM11[j][jz]=0.; ncorrMM20[j][jz]=0.; ncorrMM21[j][jz]=0.;
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
    + osSAMPM.str() + "_ntasks_" + osNTASK.str() + "_task_" + osID.str() + ".dat"; //MODIFIED HERE TO FIXED NTASKS!!!!!!
  //osNTASK.str()
  ifstream in(file_in.c_str(), ios::in);
  check_stream(in, file_in);
  while(!in.eof()) {
    vector<float> temp(3);
    in>>temp[0]>>temp[1]>>temp[2];
    //if(temp[0]==0. && temp[1]==0. && temp[2]==0.)
    //continue;
    if(temp[0]>=SIDE && temp[0]<=(SIDE+rbox) && temp[1] >= SIDE && temp[1]<=(SIDE+rbox) 
       && temp[2]>=SIDE && temp[2]<=(SIDE+rbox))
      { 
	particle OT(temp); 
	VpartM.push_back(OT); 
      }
  }
  in.close();
  cerr<<"done reading particle vector for task "<<id<<" and subcube "<<SIDE<<' '<<SIDE<<' '<<SIDE<<endl;
  cerr<<"size of imported particle vector: "<<VpartM.size()<<endl;


  //generate random data;

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
    temp[0]=SIDE + rbox*OBran1(); temp[1]=SIDE + rbox*OBran2(); temp[2]=SIDE + rbox*OBran3();
    particle Otemp(temp);
    VpartR.push_back(Otemp);
  }   
  cerr<<"size of random groups: "<<VpartR.size()<<endl;
  wtime2=MPI::Wtime();
  cerr<<"time to generate random data for slave process:"<<id<<" is: "<<wtime2-wtime1<<endl;
  
  /*string path_in = "/galformod/data/mxxl/";
  const signed long rng_seed_gen = -1381718569;
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '<<rng_seed<<endl; 
  
  cerr<<"READING PARTICLE DATA"<<endl;
  wtime1=MPI::Wtime();
  read_data OB(path_in, snap);
  vector<particle> VpartM;
  VpartM.reserve(3000000); //CAREFUL HERE, MIGHT BE BIGGER FOR SMALLER NSAMPM
  VpartM=OB.stream_particles(rng_seed, nsampM);
  cerr<<"size of formatted particle data: "<<VpartM.size()<<endl;
  ofstream out(file_out[id].c_str(), ios::out);
  for(register long long nn=0; nn<VpartM.size(); nn++) 
    out<<VpartM[nn].pos[0]<<' '<<VpartM[nn].pos[1]<<' '<<VpartM[nn].pos[2]<<endl;
  out.close();
  cerr<<"done writing particle file for id: "<<id<<endl;
  wtime2=MPI::Wtime();
  //cerr<<"time to read and format the particle data for slave process: "
  //  <<id<<" is: "<<wtime2-wtime1<<endl;
  MPI::COMM_WORLD.Barrier();
  if(id==master)
  cerr<<"time taken to compute read and write the data := "<<wtime2-wtime1<<endl<<endl; */

   
  //now build tree
    
  wtime1=MPI::Wtime();
  build_tree OtreeM(rbox, VpartM);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the M tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OM1, *OM2;
  OM1=&OtreeM; OM2=&OtreeM; 

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
  projected_correlations OcorrRM(OM1, OR2, bt, rMIN, rMAX, NB, zMIN, zMAX, NZ, flag_per, vecPER);

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
    
    if(id!=master) {
      MPI::COMM_WORLD.Send(&npairsMM[id], NB*NZ, MPI::DOUBLE, master, tag);    
      //MPI::COMM_WORLD.Send(&NN[id], 1, MPI::DOUBLE, master, tag);
      MPI::COMM_WORLD.Send(&npairsRM[id], NB*NZ, MPI::DOUBLE, master, tag);
    }
    else {      
      //now the master part
      for(int np=1; np<ntasks; np++) 
	{
	  tag=np;
	  MPI::COMM_WORLD.Recv(&npairsMM[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	  //MPI::COMM_WORLD.Recv(&NN[np], 1, MPI::DOUBLE, np, tag, status);
	  MPI::COMM_WORLD.Recv(&npairsRM[np], NB*NZ, MPI::DOUBLE, np, tag, status);
	}
      //vector<double> npairsRR=OcorrMM.generate_RR(false);
      
      for(int np=0; np<ntasks; np++) {  
	ostringstream oss;
	oss<<np;
	vector<double> Vtemp;
	file_in_RAN[np] = file_RAN + oss.str() + ".dat"; //"_XXL_PB.dat"; //".dat";
	ifstream in(file_in_RAN[np].c_str(), ios::in);
	check_stream(in, file_in_RAN[np]);
	for(int i=0; i<NB*NZ; i++) {
	  double d1, d2;
	  float r1, r2, r3, z1, z2, z3;
	  in>>r1>>r2>>r3>>z1>>z2>>z3>>d1>>d2;
	  Vtemp.push_back(d1);
	}
	in.close();
	for(int ir=0; ir<NB; ir++)
	  for(int iz=0; iz<NZ; iz++)
	    npairsRR[np][ir][iz]=Vtemp[ir*NZ+iz];
      }

      vector<double> dV=OcorrMM.volume_bins();
      vector<double> NormR(ntasks, 0.), NormM(ntasks, 0.), NormM2(ntasks, 0.), C1(ntasks), C2(ntasks);

      for(int np=0; np<ntasks; np++) 
	for(int ib=0; ib<NB; ib++) 
	  for(int jz=0; jz<NZ; jz++) 
	    //npairsMM[np][ib][jz]/npairsRR[ib*NZ+jz]
 	    ratioMM[np][ib][jz]=npairsMM[np][ib][jz]/npairsRR[np][ib][jz]-2.*npairsRM[np][ib][jz]/npairsRR[np][ib][jz];
      
      for(int np=0; np<ntasks; np++) {
	for(int ib=0; ib<NB; ib++)
	  for(int jz=0; jz<NZ; jz++) {
	    NormM[np]+=(dV[ib*NZ+jz]*npairsMM[np][ib][jz]);
	    NormR[np]+=(dV[ib*NZ+jz]*npairsRR[np][ib][jz]);
	    NormM2[np]+=(dV[ib*NZ+jz]*ratioMM[np][ib][jz]*npairsRR[np][ib][jz]);
	  }
	C1[np]=NormM[np]/NormR[np]; //1+w_vol
	C2[np]=NormM2[np]/NormR[np]; //1+w_vol
      }
      double ncorrMM2[NB][NZ];
      for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	  for(int np=0; np<ntasks; np++) {
	    ncorrMM10[ib][jz]+=(ratioMM[np][ib][jz]*C1[np]+C1[np]-1.);
	    ncorrMM11[ib][jz]+=((ratioMM[np][ib][jz]+1.)/C1[np]-1.);
	    ncorrMM20[ib][jz]+=(ratioMM[np][ib][jz]*C2[np]+C2[np]-1.);
	    ncorrMM21[ib][jz]+=((ratioMM[np][ib][jz]+1.)/C2[np]-1.);
	    DDave[ib][jz]+=npairsMM[np][ib][jz];
	    DRave[ib][jz]+=npairsRM[np][ib][jz];
	    RRave[ib][jz]+=npairsRR[np][ib][jz];
	  }
	  //compute average
	  //ncorrMM[ib][jz]/=ntasks;
	  ncorrMM10[ib][jz]/=ntasks;
	  ncorrMM11[ib][jz]/=ntasks;
	  ncorrMM20[ib][jz]/=ntasks;
	  ncorrMM21[ib][jz]/=ntasks;
	  DDave[ib][jz]/=ntasks;
	  DRave[ib][jz]/=ntasks;
	  RRave[ib][jz]/=ntasks;
	}
      
      //compute variance
      /*for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	  for(int np=0; np<ntasks; np++) 
	    nvarMM[ib][jz]+=pow(ratioMM[np][ib][jz]-ncorrMM[ib][jz], 2);
	  nvarMM[ib][jz]=sqrt(nvarMM[ib][jz]/(ntasks-1)/ntasks);
	}
      /*for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++) {
	NormR+=(dV[ib*NZ+jz]*RRave[ib][jz]);
	NormM+=(dV[ib*NZ+jz]*DDave[ib][jz]);
	cerr<<"NormM: "<<NormM<<" NormR: "<<NormR<<" Cv: "<<Cv<<endl;
      //formal average, i.e. subtract the missing 1, from -2DR/RR+RR/RR, and DR=RR
      for(int ib=0; ib<NB; ib++) 
	for(int jz=0; jz<NZ; jz++) {
	//ncorrMM[ib][jz]-=1.;
	ncorrMM[ib][jz]+=1.; //since we input the RR, and RM!=RR
	}
      double Cv2(0.);
      for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++)
	  Cv2+=(ncorrMM[ib][jz]*RRave[ib][jz]*dV[ib*NZ+jz]);
      Cv2/=NormR;
      Cv2+=1.;*/      
      cerr<<"the master writes the file for snap: "<<snap<<" MM samp: "<<nsampM<<endl;      
      ofstream outM(file_out_MM.c_str(), ios::out);
      for(int ib=0; ib<NB; ib++)
	for(int jz=0; jz<NZ; jz++) 
	  //outM<<rb[1][ib]<<' '<<zb[1][jz]<<' '<<ncorrMM[ib][jz]<<' '<<DDave[ib][jz]<<' '
	  //<<DRave[ib][jz]<<' '<<RRave[ib][jz]<<' '<<Cv<<' '<<Cv2<<endl;
	  outM<<rb[1][ib]<<' '<<zb[1][jz]<<' '<<ncorrMM10[ib][jz]<<' '<<ncorrMM11[ib][jz]<<' '
	      <<ncorrMM20[ib][jz]<<' '<<ncorrMM21[ib][jz]<<endl;
	    //ncorrMM[ib][jz]*Cv+(Cv-1)<<' '<<(ncorrMM[ib][jz]+1)/Cv-1.<<' '<<Cv<<' '<<
      outM.close();
      cerr<<"\n the master has written the MM output file " + file_out_MM<<endl;
      
    }
    
    MPI::COMM_WORLD.Barrier();
    FinTime=MPI::Wtime();
    cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
    MPI::Finalize();
    
    return 0;
  
}

 

